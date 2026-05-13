/*
 * TrackIndexWriter.h
 *
 * Helper for writing track.idx files in misha's indexed format (see
 * TrackIndex.h for the spec). Encapsulates header write, entry append,
 * and post-hoc checksum patching so producers (the pack-per-chrom path
 * and the upcoming streaming direct-write path) share one definition
 * of the on-disk format.
 */

#ifndef TRACKINDEXWRITER_H_
#define TRACKINDEXWRITER_H_

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "TrackIndex.h"

namespace rdb {

class TrackIndexWriter {
public:
    // Write the (provisional) header with checksum=0 into idx_fp.
    // Caller is responsible for fseek/fopen; this just writes 36 bytes.
    // Throws TGLException on write failure.
    static void write_header(FILE *idx_fp,
                             MishaTrackType type,
                             std::uint32_t num_contigs);

    // Append one entry (24 bytes) to idx_fp at the current position.
    // Field-by-field fwrite preserves the exact byte layout used by the
    // existing pack/convert writers (no struct padding leaks on disk).
    // Throws TGLException on write failure.
    static void write_entry(FILE *idx_fp, const TrackContigEntry &e);

    // After all entries are written, compute CRC64-ECMA over the entry
    // table (NumContigs * 20 bytes: chrom_id + offset + length per entry,
    // matching the existing pack/convert checksums - reserved field is
    // intentionally NOT included) and patch it into the header's
    // checksum field at offset IDX_HEADER_SIZE_TO_CHECKSUM.
    // Throws TGLException on seek/write failure.
    static void finalize_checksum(FILE *idx_fp,
                                  const std::vector<TrackContigEntry> &entries);
};

} // namespace rdb

#endif /* TRACKINDEXWRITER_H_ */
