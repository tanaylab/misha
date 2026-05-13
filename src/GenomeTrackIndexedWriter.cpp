/*
 * GenomeTrackIndexedWriter.cpp
 *
 * Implementation of the streaming indexed writer. Produces output that
 * is bit-for-bit identical to what gtrack_pack_per_chrom_to_indexed
 * would have produced from the equivalent per-chrom files:
 *
 *   - track.dat is the byte-for-byte concatenation, in chromid order
 *     (0..total_contigs-1), of every per-chrom file's bytes. Chroms
 *     that never called begin_chrom contribute zero bytes (matching
 *     pack's behaviour when a per-chrom file is absent).
 *   - track.idx uses the same on-disk layout (see TrackIndex.h) and
 *     the same CRC64 checksum domain (chrom_id + offset + length only,
 *     reserved excluded), via TrackIndexWriter.
 *
 * The writer accepts begin_chrom() calls in any order (the scanner
 * emits chroms strictly in ascending chromid order, but we don't rely
 * on that here - finalize() sorts the recorded entries and pads the
 * gaps with length=0 entries pointing at the current track.dat end,
 * matching pack's "missing per-chrom file" behaviour).
 */

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "GenomeTrackIndexedWriter.h"
#include "TGLException.h"
#include "TrackIndexWriter.h"

using namespace std;

namespace rdb {

GenomeTrackIndexedWriter::~GenomeTrackIndexedWriter()
{
    if (m_dat_fp) {
        // Close without renaming the .idx; caller forgot to finalize().
        // Leftover track.dat will be cleaned up by the R-level atomic
        // create wrapper (which trashes the whole track dir on error).
        fclose(m_dat_fp);
        m_dat_fp = nullptr;
    }
}

MishaTrackType GenomeTrackIndexedWriter::mtype_for(GenomeTrack::Type type)
{
    switch (type) {
        case GenomeTrack::FIXED_BIN: return MishaTrackType::DENSE;
        case GenomeTrack::SPARSE:    return MishaTrackType::SPARSE;
        case GenomeTrack::ARRAYS:    return MishaTrackType::ARRAY;
        default:
            TGLError<GenomeTrack>("Unsupported 1D track type for indexed writer");
    }
    // unreachable
    return MishaTrackType::DENSE;
}

void GenomeTrackIndexedWriter::init(const string &track_dir,
                                    GenomeTrack::Type type,
                                    uint32_t total_contigs)
{
    if (m_dat_fp)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::init called twice");

    // Validate type early (throws if unsupported).
    (void)mtype_for(type);

    m_track_dir = track_dir;
    m_type = type;
    m_total_contigs = total_contigs;
    m_current_offset = 0;
    m_current_chrom = -1;
    m_chrom_start_offset = 0;
    m_entries.clear();
    m_finalized = false;

    const string dat_path = m_track_dir + "/track.dat";
    m_dat_fp = fopen(dat_path.c_str(), "wb");
    if (!m_dat_fp)
        TGLError<GenomeTrack>("Failed to create %s: %s", dat_path.c_str(), strerror(errno));
}

void GenomeTrackIndexedWriter::begin_chrom(int chromid)
{
    if (!m_dat_fp)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::begin_chrom before init");
    if (m_current_chrom != -1)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::begin_chrom called while chrom %d is still open",
                              m_current_chrom);

    m_current_chrom = chromid;
    m_chrom_start_offset = m_current_offset;
}

void GenomeTrackIndexedWriter::write_bytes(const void *buf, size_t len)
{
    if (!m_dat_fp)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::write_bytes before init");
    if (m_current_chrom == -1)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::write_bytes with no chrom open");
    if (len == 0)
        return;

    if (fwrite(buf, 1, len, m_dat_fp) != len)
        TGLError<GenomeTrack>("Failed to write track.dat: %s", strerror(errno));
    m_current_offset += len;
}

void GenomeTrackIndexedWriter::end_chrom()
{
    if (m_current_chrom == -1)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::end_chrom with no chrom open");

    ChromEntry e;
    e.chromid = m_current_chrom;
    e.offset  = m_chrom_start_offset;
    e.length  = m_current_offset - m_chrom_start_offset;
    m_entries.push_back(e);

    m_current_chrom = -1;
}

void GenomeTrackIndexedWriter::finalize()
{
    if (m_finalized)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::finalize called twice");
    if (!m_dat_fp)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::finalize before init");
    if (m_current_chrom != -1)
        TGLError<GenomeTrack>("GenomeTrackIndexedWriter::finalize while chrom %d is still open",
                              m_current_chrom);

    // Flush + fsync + close track.dat first. Index references positions
    // inside this file, so it must be durable before track.idx is named.
    if (fflush(m_dat_fp) != 0) {
        const int err = errno;
        fclose(m_dat_fp);
        m_dat_fp = nullptr;
        TGLError<GenomeTrack>("Failed to flush track.dat: %s", strerror(err));
    }
    if (fsync(fileno(m_dat_fp)) != 0) {
        const int err = errno;
        fclose(m_dat_fp);
        m_dat_fp = nullptr;
        TGLError<GenomeTrack>("Failed to fsync track.dat: %s", strerror(err));
    }
    if (fclose(m_dat_fp) != 0) {
        m_dat_fp = nullptr;
        TGLError<GenomeTrack>("Failed to close track.dat: %s", strerror(errno));
    }
    m_dat_fp = nullptr;

    // Build the full per-genome entry table. One TrackContigEntry per
    // chromid in [0, total_contigs); missing chroms get length=0 with
    // offset = end-of-dat (matches pack's behaviour for absent files).
    vector<TrackContigEntry> full_entries;
    full_entries.reserve(m_total_contigs);

    // Index recorded entries by chromid (recorded order may equal
    // genome order in practice; this is defensive).
    vector<const ChromEntry *> by_id(m_total_contigs, nullptr);
    for (const ChromEntry &e : m_entries) {
        if (e.chromid < 0 || (uint32_t)e.chromid >= m_total_contigs)
            TGLError<GenomeTrack>("Recorded chromid %d outside [0,%u)",
                                  e.chromid, m_total_contigs);
        if (by_id[e.chromid])
            TGLError<GenomeTrack>("Chromid %d written twice in indexed writer", e.chromid);
        by_id[e.chromid] = &e;
    }

    for (uint32_t chromid = 0; chromid < m_total_contigs; ++chromid) {
        TrackContigEntry te;
        te.chrom_id = chromid;
        te.reserved = 0;
        if (by_id[chromid]) {
            te.offset = by_id[chromid]->offset;
            te.length = by_id[chromid]->length;
        } else {
            te.offset = m_current_offset; // end of dat
            te.length = 0;
        }
        full_entries.push_back(te);
    }

    // Write track.idx via tmp + rename.
    const string idx_path     = m_track_dir + "/track.idx";
    const string idx_path_tmp = idx_path + ".tmp";

    FILE *idx_fp = fopen(idx_path_tmp.c_str(), "wb");
    if (!idx_fp)
        TGLError<GenomeTrack>("Failed to create %s: %s", idx_path_tmp.c_str(), strerror(errno));

    try {
        TrackIndexWriter::write_header(idx_fp, mtype_for(m_type), m_total_contigs);
        for (const TrackContigEntry &te : full_entries)
            TrackIndexWriter::write_entry(idx_fp, te);
        TrackIndexWriter::finalize_checksum(idx_fp, full_entries);
    } catch (...) {
        fclose(idx_fp);
        unlink(idx_path_tmp.c_str());
        throw;
    }

    if (fflush(idx_fp) != 0) {
        const int err = errno;
        fclose(idx_fp);
        unlink(idx_path_tmp.c_str());
        TGLError<GenomeTrack>("Failed to flush %s: %s", idx_path_tmp.c_str(), strerror(err));
    }
    if (fsync(fileno(idx_fp)) != 0) {
        const int err = errno;
        fclose(idx_fp);
        unlink(idx_path_tmp.c_str());
        TGLError<GenomeTrack>("Failed to fsync %s: %s", idx_path_tmp.c_str(), strerror(err));
    }
    if (fclose(idx_fp) != 0) {
        const int err = errno;
        unlink(idx_path_tmp.c_str());
        TGLError<GenomeTrack>("Failed to close %s: %s", idx_path_tmp.c_str(), strerror(err));
    }

    if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0) {
        const int err = errno;
        unlink(idx_path_tmp.c_str());
        TGLError<GenomeTrack>("Failed to rename %s to %s: %s",
                              idx_path_tmp.c_str(), idx_path.c_str(), strerror(err));
    }

    m_finalized = true;
}

} // namespace rdb
