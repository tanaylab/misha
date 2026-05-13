/*
 * GenomeTrackIndexedWriter.h
 *
 * Streams 1D track data directly into track.dat + track.idx (indexed
 * format) on indexed-format DBs, avoiding the per-chrom files that
 * gtrack_pack_per_chrom_to_indexed would otherwise eliminate as a
 * second pass.
 *
 * Usage:
 *   GenomeTrackIndexedWriter w;
 *   w.init(track_dir, GenomeTrack::FIXED_BIN, total_contigs);
 *   for each chrom with data:
 *       w.begin_chrom(chromid);
 *       w.write_bytes(buf, len);     // same bytes that would go into
 *                                    // a per-chrom file for this type
 *       w.end_chrom();
 *   w.finalize();                    // writes track.idx with checksum
 *
 * total_contigs must equal the genome's chromosome count so that
 * track.idx contains one entry per chromid in genome order (entries
 * for chroms that never called begin_chrom get length=0 with offset
 * pointing at the current track.dat end).
 *
 * Throws TGLException on I/O errors. On error track.dat may be
 * partial; caller is responsible for cleanup (typically .gdb.trash
 * via the R-level atomic-create wrapper).
 */

#ifndef GENOMETRACKINDEXEDWRITER_H_
#define GENOMETRACKINDEXEDWRITER_H_

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "GenomeTrack.h"
#include "TrackIndex.h"

namespace rdb {

class GenomeTrackIndexedWriter {
public:
    GenomeTrackIndexedWriter() = default;
    ~GenomeTrackIndexedWriter();

    // Non-copyable (owns FILE*).
    GenomeTrackIndexedWriter(const GenomeTrackIndexedWriter &) = delete;
    GenomeTrackIndexedWriter &operator=(const GenomeTrackIndexedWriter &) = delete;

    void init(const std::string &track_dir,
              GenomeTrack::Type type,
              std::uint32_t total_contigs);

    void begin_chrom(int chromid);
    void write_bytes(const void *buf, std::size_t len);
    void end_chrom();

    void finalize();

private:
    struct ChromEntry {
        int chromid;
        std::uint64_t offset;
        std::uint64_t length;
    };

    static MishaTrackType mtype_for(GenomeTrack::Type type);

    std::string m_track_dir;
    GenomeTrack::Type m_type = GenomeTrack::NUM_TYPES;
    std::uint32_t m_total_contigs = 0;

    FILE *m_dat_fp = nullptr;
    std::uint64_t m_current_offset = 0;

    int m_current_chrom = -1;
    std::uint64_t m_chrom_start_offset = 0;

    // Largest chromid passed to begin_chrom() so far. begin_chrom()
    // requires strictly increasing chromids so that track.dat is laid
    // out in genome order; finalize() relies on this to compute
    // missing-chrom offsets identically to the legacy pack
    // (gtrack_pack_per_chrom_to_indexed) - see finalize() comment.
    int m_max_seen_chromid = -1;

    // Entries recorded in the order begin_chrom() was called. Sparse
    // (only chroms with data are present); finalize() pads to a full
    // per-genome table before writing track.idx.
    std::vector<ChromEntry> m_entries;
    bool m_finalized = false;
};

} // namespace rdb

#endif /* GENOMETRACKINDEXEDWRITER_H_ */
