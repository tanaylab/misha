/*
 * GenomeTrackIndexedFormat.cpp
 *
 * Converts per-chromosome tracks to indexed format
 */

#include <cstdint>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>

#include "GenomeTrack.h"
#include "TrackIndex.h"
#include "CRC64.h"
#include "BufferedFile.h"
#include "TGLException.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// Offset to checksum field in index header
static const size_t IDX_HEADER_SIZE_TO_CHECKSUM =
    8 +                    // Magic header
    sizeof(uint32_t) +     // Version
    sizeof(uint32_t) +     // Track type
    sizeof(uint32_t) +     // Num contigs
    sizeof(uint64_t);      // Flags

// Helper function to write index header
static void write_index_header(FILE *fp, MishaTrackType track_type, uint32_t num_contigs, uint64_t checksum) {
    // Magic header
    const char magic[8] = {'M','I','S','H','A','T','D','X'};
    if (fwrite(magic, 1, 8, fp) != 8) {
        TGLError<GenomeTrack>("Failed to write index header");
    }

    // Version
    uint32_t version = 1;
    if (fwrite(&version, sizeof(version), 1, fp) != 1) {
        TGLError<GenomeTrack>("Failed to write index version");
    }

    // Track type
    uint32_t track_type_raw = static_cast<uint32_t>(track_type);
    if (fwrite(&track_type_raw, sizeof(track_type_raw), 1, fp) != 1) {
        TGLError<GenomeTrack>("Failed to write track type");
    }

    // Number of contigs
    if (fwrite(&num_contigs, sizeof(num_contigs), 1, fp) != 1) {
        TGLError<GenomeTrack>("Failed to write number of contigs");
    }

    // Flags (little-endian flag)
    uint64_t flags = 0x01; // IS_LITTLE_ENDIAN
    if (fwrite(&flags, sizeof(flags), 1, fp) != 1) {
        TGLError<GenomeTrack>("Failed to write flags");
    }

    // Checksum
    if (fwrite(&checksum, sizeof(checksum), 1, fp) != 1) {
        TGLError<GenomeTrack>("Failed to write checksum");
    }
}

// Helper function to copy file contents
static bool copy_file_contents(const string &src, FILE *dest, uint64_t &bytes_written) {
    FILE *src_fp = fopen(src.c_str(), "rb");
    if (!src_fp) {
        return false; // File doesn't exist or can't be opened
    }

    // Get file size
    if (fseek(src_fp, 0, SEEK_END) != 0) {
        fclose(src_fp);
        TGLError<GenomeTrack>("Failed to seek to end of %s", src.c_str());
    }
    uint64_t file_size = ftello(src_fp);
    if (fseek(src_fp, 0, SEEK_SET) != 0) {
        fclose(src_fp);
        TGLError<GenomeTrack>("Failed to seek to start of %s", src.c_str());
    }

    // Copy in chunks
    const size_t BUFFER_SIZE = 1024 * 1024; // 1MB buffer
    char *buffer = new char[BUFFER_SIZE];
    uint64_t total_read = 0;

    while (total_read < file_size) {
        size_t to_read = min((uint64_t)BUFFER_SIZE, file_size - total_read);
        size_t read_bytes = fread(buffer, 1, to_read, src_fp);
        if (read_bytes != to_read) {
            delete[] buffer;
            fclose(src_fp);
            TGLError<GenomeTrack>("Failed to read from %s", src.c_str());
        }

        size_t written = fwrite(buffer, 1, read_bytes, dest);
        if (written != read_bytes) {
            delete[] buffer;
            fclose(src_fp);
            TGLError<GenomeTrack>("Failed to write to track.dat");
        }

        total_read += read_bytes;
    }

    delete[] buffer;
    fclose(src_fp);
    bytes_written = file_size;
    return true;
}

extern "C" {

SEXP gtrack_convert_to_indexed_format(SEXP _track, SEXP _remove_old, SEXP _envir) {
    // Declare paths outside try block for cleanup access
    string dat_path_tmp;
    string idx_path_tmp;

    try {
        RdbInitializer rdb_init;

        // Parse arguments
        if (!Rf_isString(_track) || Rf_length(_track) != 1) {
            verror("Track name must be a string");
        }
        string track = CHAR(STRING_ELT(_track, 0));

        bool remove_old = Rf_asLogical(_remove_old);

        // Get track directory
        const char *groot = get_groot(_envir);
        if (!groot || !*groot) {
            verror("Genome root is not defined");
        }

        string track_dir_name = track;
        replace(track_dir_name.begin(), track_dir_name.end(), '.', '/');
        string track_dir = string(groot) + "/tracks/" + track_dir_name + ".track";

        // Get chromosome list
        IntervUtils iu(_envir);
        const GenomeChromKey &chromkey = iu.get_chromkey();

        // Determine track type from first chromosome file
        GenomeTrack::Type track_type_enum = GenomeTrack::get_type(track_dir.c_str(), chromkey);
        MishaTrackType track_type = MishaTrackType::DENSE;

        switch (track_type_enum) {
            case GenomeTrack::FIXED_BIN:
                track_type = MishaTrackType::DENSE;
                break;
            case GenomeTrack::SPARSE:
                track_type = MishaTrackType::SPARSE;
                break;
            case GenomeTrack::ARRAYS:
                track_type = MishaTrackType::ARRAY;
                break;
            default:
                verror("Only 1D tracks (dense, sparse, array) can be converted");
        }

        // Prepare paths
        dat_path_tmp = track_dir + "/track.dat.tmp";
        idx_path_tmp = track_dir + "/track.idx.tmp";
        string dat_path = track_dir + "/track.dat";
        string idx_path = track_dir + "/track.idx";

        // Open temporary files
        FILE *dat_fp = fopen(dat_path_tmp.c_str(), "wb");
        if (!dat_fp) {
            TGLError<GenomeTrack>("Failed to create %s: %s", dat_path_tmp.c_str(), strerror(errno));
        }

        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            TGLError<GenomeTrack>("Failed to create %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Write index header (checksum=0 for now)
        write_index_header(idx_fp, track_type, chromkey.get_num_chroms(), 0);

        // Collect entries and concatenate files
        vector<TrackContigEntry> entries;
        vector<string> chr_files_to_remove;
        uint64_t current_offset = 0;

        for (int chromid = 0; chromid < (int)chromkey.get_num_chroms(); chromid++) {
            string chrom_name = chromkey.id2chrom(chromid);

            // Try to find the chromosome file, handling chr prefix mismatch
            string chr_file;
            bool found = false;

            // Try 1: chromosome name as-is
            string candidate = track_dir + "/" + chrom_name;
            struct stat st;
            if (stat(candidate.c_str(), &st) == 0) {
                chr_file = candidate;
                found = true;
            } else {
                // Try 2: with "chr" prefix if not already present
                if (chrom_name.substr(0, 3) != "chr") {
                    candidate = track_dir + "/chr" + chrom_name;
                    if (stat(candidate.c_str(), &st) == 0) {
                        chr_file = candidate;
                        found = true;
                    }
                } else {
                    // Try 3: without "chr" prefix if present
                    candidate = track_dir + "/" + chrom_name.substr(3);
                    if (stat(candidate.c_str(), &st) == 0) {
                        chr_file = candidate;
                        found = true;
                    }
                }
            }

            TrackContigEntry entry;
            entry.chrom_id = chromid;
            entry.offset = current_offset;
            entry.length = 0;
            entry.reserved = 0;

            uint64_t bytes_written = 0;
            if (found && copy_file_contents(chr_file, dat_fp, bytes_written)) {
                entry.length = bytes_written;
                current_offset += bytes_written;
                chr_files_to_remove.push_back(chr_file);
            }

            // Write entry to index
            if (fwrite(&entry.chrom_id, sizeof(entry.chrom_id), 1, idx_fp) != 1 ||
                fwrite(&entry.offset, sizeof(entry.offset), 1, idx_fp) != 1 ||
                fwrite(&entry.length, sizeof(entry.length), 1, idx_fp) != 1 ||
                fwrite(&entry.reserved, sizeof(entry.reserved), 1, idx_fp) != 1) {
                fclose(dat_fp);
                fclose(idx_fp);
                TGLError<GenomeTrack>("Failed to write index entry for chromosome %s", chrom_name.c_str());
            }

            entries.push_back(entry);
        }

        // Compute checksum of entries
        misha::CRC64 crc64;
        uint64_t checksum = crc64.init_incremental();
        for (const auto &entry : entries) {
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.chrom_id, sizeof(entry.chrom_id));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.offset, sizeof(entry.offset));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.length, sizeof(entry.length));
        }
        checksum = crc64.finalize_incremental(checksum);

        // Update checksum in index header
        if (fseek(idx_fp, IDX_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError<GenomeTrack>("Failed to seek to checksum position in index");
        }
        if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError<GenomeTrack>("Failed to update checksum in index");
        }

        // Flush and sync both files
        fflush(dat_fp);
        fflush(idx_fp);
        fsync(fileno(dat_fp));
        fsync(fileno(idx_fp));

        fclose(dat_fp);
        fclose(idx_fp);

        // Atomic rename: track.dat first, then track.idx
        if (rename(dat_path_tmp.c_str(), dat_path.c_str()) != 0) {
            TGLError<GenomeTrack>("Failed to rename %s to %s: %s",
                dat_path_tmp.c_str(), dat_path.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0) {
            TGLError<GenomeTrack>("Failed to rename %s to %s: %s",
                idx_path_tmp.c_str(), idx_path.c_str(), strerror(errno));
        }

        // Validate conversion before removing old files
        // Check that track.dat has the expected size
        struct stat dat_stat;
        if (stat(dat_path.c_str(), &dat_stat) != 0) {
            TGLError<GenomeTrack>("Failed to stat %s after conversion", dat_path.c_str());
        }

        if ((uint64_t)dat_stat.st_size != current_offset) {
            TGLError<GenomeTrack>("track.dat size mismatch: expected %llu bytes, got %llu bytes",
                (unsigned long long)current_offset, (unsigned long long)dat_stat.st_size);
        }

        // Remove old per-chromosome files if requested
        if (remove_old) {
            for (const string &chr_file : chr_files_to_remove) {
                unlink(chr_file.c_str());
            }
        }

        return R_NilValue;

    } catch (TGLException &e) {
        // Clean up temporary files on error
        unlink(dat_path_tmp.c_str());
        unlink(idx_path_tmp.c_str());
        verror("%s", e.msg());
    } catch (const bad_alloc &e) {
        // Clean up temporary files on error
        unlink(dat_path_tmp.c_str());
        unlink(idx_path_tmp.c_str());
        verror("Out of memory");
    }

    return R_NilValue;
}

} // extern "C"
