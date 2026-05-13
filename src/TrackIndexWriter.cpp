/*
 * TrackIndexWriter.cpp
 *
 * Implementation of TrackIndexWriter - the single source of truth for
 * the on-disk layout of misha's track.idx files. Bodies are lifted
 * verbatim from gtrack_pack_per_chrom_to_indexed (in
 * GenomeTrackSplitIndexed.cpp) and gtrack_convert_to_indexed_format
 * (in GenomeTrackIndexedFormat.cpp); see TrackIndex.h for the format
 * spec.
 */

#include <cstdint>
#include <cstdio>
#include <vector>

#include "TrackIndexWriter.h"
#include "TrackIndex.h"
#include "CRC64.h"
#include "GenomeTrack.h"
#include "TGLException.h"

using namespace std;

namespace rdb {

void TrackIndexWriter::write_header(FILE *idx_fp,
                                    MishaTrackType type,
                                    uint32_t num_contigs)
{
    // Layout (36 bytes total), in this exact order:
    //   magic[8] + version u32 + track_type u32 + num_contigs u32
    //     + flags u64 + checksum u64 (placeholder = 0; patched later)
    const char magic[8] = {'M','I','S','H','A','T','D','X'};
    const uint32_t version = 1;
    const uint32_t track_type_raw = static_cast<uint32_t>(type);
    const uint64_t flags = 0x01; // IS_LITTLE_ENDIAN
    const uint64_t checksum_placeholder = 0;

    bool ok =
        fwrite(magic, 1, 8, idx_fp) == 8 &&
        fwrite(&version, sizeof(version), 1, idx_fp) == 1 &&
        fwrite(&track_type_raw, sizeof(track_type_raw), 1, idx_fp) == 1 &&
        fwrite(&num_contigs, sizeof(num_contigs), 1, idx_fp) == 1 &&
        fwrite(&flags, sizeof(flags), 1, idx_fp) == 1 &&
        fwrite(&checksum_placeholder, sizeof(checksum_placeholder), 1, idx_fp) == 1;
    if (!ok) {
        TGLError<GenomeTrack>("Failed to write index header");
    }
}

void TrackIndexWriter::write_entry(FILE *idx_fp, const TrackContigEntry &e)
{
    // Field-by-field fwrite (no struct fwrite) - this preserves the
    // exact 24-byte layout (4 + 8 + 8 + 4) without any host-specific
    // struct padding leaking onto disk.
    if (fwrite(&e.chrom_id, sizeof(e.chrom_id), 1, idx_fp) != 1 ||
        fwrite(&e.offset,   sizeof(e.offset),   1, idx_fp) != 1 ||
        fwrite(&e.length,   sizeof(e.length),   1, idx_fp) != 1 ||
        fwrite(&e.reserved, sizeof(e.reserved), 1, idx_fp) != 1) {
        TGLError<GenomeTrack>("Failed to write index entry");
    }
}

void TrackIndexWriter::finalize_checksum(FILE *idx_fp,
                                         const vector<TrackContigEntry> &entries)
{
    // CRC64-ECMA over (chrom_id, offset, length) of every entry, in
    // order. The `reserved` field is intentionally excluded to match
    // the existing pack/convert writers' checksum domain (and what
    // TrackIndex::load validates).
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();
    for (const TrackContigEntry &e : entries) {
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&e.chrom_id, sizeof(e.chrom_id));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&e.offset,   sizeof(e.offset));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&e.length,   sizeof(e.length));
    }
    checksum = crc64.finalize_incremental(checksum);

    // Patch checksum field in header (offset documented in TrackIndex.h).
    if (fseek(idx_fp, IDX_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
        TGLError<GenomeTrack>("Failed to seek to checksum position in index");
    }
    if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
        TGLError<GenomeTrack>("Failed to update checksum in index");
    }
}

} // namespace rdb
