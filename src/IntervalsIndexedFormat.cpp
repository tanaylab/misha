/*
 * IntervalsIndexedFormat.cpp
 *
 * Converts per-chromosome/per-pair intervals to indexed format
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
#include <dirent.h>

#include "IntervalsIndex1D.h"
#include "IntervalsIndex2D.h"
#include "CRC64.h"
#include "TGLException.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// Offset to checksum field in 1D index header
static const size_t IDX_1D_HEADER_SIZE_TO_CHECKSUM =
    8 +                    // Magic header
    sizeof(uint32_t) +     // Version
    sizeof(uint32_t) +     // Num entries
    sizeof(uint64_t);      // Flags

// Offset to checksum field in 2D index header
static const size_t IDX_2D_HEADER_SIZE_TO_CHECKSUM =
    8 +                    // Magic header
    sizeof(uint32_t) +     // Version
    sizeof(uint32_t) +     // Num entries
    sizeof(uint64_t);      // Flags

// Helper function to write 1D index header
static void write_1d_index_header(FILE *fp, uint32_t num_entries, uint64_t checksum) {
    // Magic header
    const char magic[8] = {'M','I','S','H','A','I','1','D'};
    if (fwrite(magic, 1, 8, fp) != 8) {
        TGLError("Failed to write index header");
    }

    // Version
    uint32_t version = 1;
    if (fwrite(&version, sizeof(version), 1, fp) != 1) {
        TGLError("Failed to write index version");
    }

    // Number of entries
    if (fwrite(&num_entries, sizeof(num_entries), 1, fp) != 1) {
        TGLError("Failed to write number of entries");
    }

    // Flags (little-endian flag)
    uint64_t flags = 0x01; // IS_LITTLE_ENDIAN
    if (fwrite(&flags, sizeof(flags), 1, fp) != 1) {
        TGLError("Failed to write flags");
    }

    // Checksum
    if (fwrite(&checksum, sizeof(checksum), 1, fp) != 1) {
        TGLError("Failed to write checksum");
    }

    // Reserved
    uint32_t reserved = 0;
    if (fwrite(&reserved, sizeof(reserved), 1, fp) != 1) {
        TGLError("Failed to write reserved field");
    }
}

// Helper function to write 2D index header
static void write_2d_index_header(FILE *fp, uint32_t num_entries, uint64_t checksum) {
    // Magic header
    const char magic[8] = {'M','I','S','H','A','I','2','D'};
    if (fwrite(magic, 1, 8, fp) != 8) {
        TGLError("Failed to write index header");
    }

    // Version
    uint32_t version = 1;
    if (fwrite(&version, sizeof(version), 1, fp) != 1) {
        TGLError("Failed to write index version");
    }

    // Number of entries
    if (fwrite(&num_entries, sizeof(num_entries), 1, fp) != 1) {
        TGLError("Failed to write number of entries");
    }

    // Flags (little-endian flag)
    uint64_t flags = 0x01; // IS_LITTLE_ENDIAN
    if (fwrite(&flags, sizeof(flags), 1, fp) != 1) {
        TGLError("Failed to write flags");
    }

    // Checksum
    if (fwrite(&checksum, sizeof(checksum), 1, fp) != 1) {
        TGLError("Failed to write checksum");
    }

    // Reserved
    uint64_t reserved = 0;
    if (fwrite(&reserved, sizeof(reserved), 1, fp) != 1) {
        TGLError("Failed to write reserved field");
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
        TGLError("Failed to seek to end of %s", src.c_str());
    }
    uint64_t file_size = ftello(src_fp);
    if (fseek(src_fp, 0, SEEK_SET) != 0) {
        fclose(src_fp);
        TGLError("Failed to seek to start of %s", src.c_str());
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
            TGLError("Failed to read from %s", src.c_str());
        }

        size_t written = fwrite(buffer, 1, read_bytes, dest);
        if (written != read_bytes) {
            delete[] buffer;
            fclose(src_fp);
            TGLError("Failed to write to data file");
        }

        total_read += read_bytes;
    }

    delete[] buffer;
    fclose(src_fp);
    bytes_written = file_size;
    return true;
}

extern "C" {

SEXP ginterv_convert(SEXP _intervset, SEXP _remove_old, SEXP _envir) {
    // Declare paths outside try block for cleanup access
    string dat_path_tmp;
    string idx_path_tmp;

    try {
        RdbInitializer rdb_init;

        // Parse arguments
        if (!Rf_isString(_intervset) || Rf_length(_intervset) != 1) {
            verror("Interval set name must be a string");
        }
        string intervset = CHAR(STRING_ELT(_intervset, 0));

        bool remove_old = Rf_asLogical(_remove_old);

        // Get interval set directory
        string intervset_dir = interv2path(_envir, intervset.c_str());

        // Get chromosome list
        IntervUtils iu(_envir);
        const GenomeChromKey &chromkey = iu.get_chromkey();

        // Prepare paths
        dat_path_tmp = intervset_dir + "/intervals.dat.tmp";
        idx_path_tmp = intervset_dir + "/intervals.idx.tmp";
        string dat_path = intervset_dir + "/intervals.dat";
        string idx_path = intervset_dir + "/intervals.idx";

        // Open temporary files
        FILE *dat_fp = fopen(dat_path_tmp.c_str(), "wb");
        if (!dat_fp) {
            TGLError("Failed to create %s: %s", dat_path_tmp.c_str(), strerror(errno));
        }

        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            TGLError("Failed to create %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Write index header (checksum=0 for now)
        write_1d_index_header(idx_fp, chromkey.get_num_chroms(), 0);

        // Collect entries and concatenate files
        vector<IntervalsContigEntry> entries;
        vector<string> chr_files_to_remove;
        uint64_t current_offset = 0;

        for (int chromid = 0; chromid < (int)chromkey.get_num_chroms(); chromid++) {
            string chrom_name = chromkey.id2chrom(chromid);
            string chr_file = intervset_dir + "/" + chrom_name;

            IntervalsContigEntry entry;
            entry.chrom_id = chromid;
            entry.offset = current_offset;
            entry.length = 0;
            entry.reserved = 0;

            uint64_t bytes_written = 0;
            if (copy_file_contents(chr_file, dat_fp, bytes_written)) {
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
                TGLError("Failed to write index entry for chromosome %s", chrom_name.c_str());
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

        // Update checksum in index header (seek to checksum position)
        if (fseek(idx_fp, IDX_1D_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError("Failed to seek to checksum position in 1D index");
        }
        if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError("Failed to update checksum in index");
        }

        // Flush and sync both files
        fflush(dat_fp);
        fflush(idx_fp);
        fsync(fileno(dat_fp));
        fsync(fileno(idx_fp));

        fclose(dat_fp);
        fclose(idx_fp);

        // Atomic rename: intervals.dat first, then intervals.idx
        if (rename(dat_path_tmp.c_str(), dat_path.c_str()) != 0) {
            TGLError("Failed to rename %s to %s: %s",
                dat_path_tmp.c_str(), dat_path.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0) {
            TGLError("Failed to rename %s to %s: %s",
                idx_path_tmp.c_str(), idx_path.c_str(), strerror(errno));
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

SEXP ginterv2d_convert(SEXP _intervset, SEXP _remove_old, SEXP _envir) {
    // Declare paths outside try block for cleanup access
    string dat_path_tmp;
    string idx_path_tmp;

    try {
        RdbInitializer rdb_init;

        // Parse arguments
        if (!Rf_isString(_intervset) || Rf_length(_intervset) != 1) {
            verror("Interval set name must be a string");
        }
        string intervset = CHAR(STRING_ELT(_intervset, 0));

        bool remove_old = Rf_asLogical(_remove_old);

        // Get interval set directory
        string intervset_dir = interv2path(_envir, intervset.c_str());

        // Get chromosome list
        IntervUtils iu(_envir);
        const GenomeChromKey &chromkey = iu.get_chromkey();

        // Enumerate existing per-pair files
        DIR *dir = opendir(intervset_dir.c_str());
        if (!dir) {
            verror("Cannot open interval set directory %s: %s",
                   intervset_dir.c_str(), strerror(errno));
        }

        vector<pair<int, int>> pairs; // (chromid1, chromid2)
        struct dirent *entry;
        while ((entry = readdir(dir)) != nullptr) {
            string filename = entry->d_name;

            // Skip . and ..
            if (filename == "." || filename == "..") continue;

            // Look for files matching pattern: chrom1-chrom2
            size_t dash_pos = filename.find('-');
            if (dash_pos != string::npos && dash_pos > 0 && dash_pos < filename.length() - 1) {
                string chrom1_name = filename.substr(0, dash_pos);
                string chrom2_name = filename.substr(dash_pos + 1);

                // Try to map to chromids
                int chromid1 = chromkey.chrom2id(chrom1_name.c_str());
                int chromid2 = chromkey.chrom2id(chrom2_name.c_str());

                if (chromid1 >= 0 && chromid2 >= 0) {
                    pairs.push_back(make_pair(chromid1, chromid2));
                }
            }
        }
        closedir(dir);

        // Sort pairs by (chromid1, chromid2) for stable ordering
        sort(pairs.begin(), pairs.end());

        // Prepare paths
        dat_path_tmp = intervset_dir + "/intervals2d.dat.tmp";
        idx_path_tmp = intervset_dir + "/intervals2d.idx.tmp";
        string dat_path = intervset_dir + "/intervals2d.dat";
        string idx_path = intervset_dir + "/intervals2d.idx";

        // Open temporary files
        FILE *dat_fp = fopen(dat_path_tmp.c_str(), "wb");
        if (!dat_fp) {
            TGLError("Failed to create %s: %s", dat_path_tmp.c_str(), strerror(errno));
        }

        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            TGLError("Failed to create %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Write index header (checksum=0 for now)
        write_2d_index_header(idx_fp, pairs.size(), 0);

        // Collect entries and concatenate files
        vector<IntervalsPairEntry> entries;
        vector<string> pair_files_to_remove;
        uint64_t current_offset = 0;

        for (const auto &pair : pairs) {
            int chromid1 = pair.first;
            int chromid2 = pair.second;
            string chrom1_name = chromkey.id2chrom(chromid1);
            string chrom2_name = chromkey.id2chrom(chromid2);
            string pair_file = intervset_dir + "/" + chrom1_name + "-" + chrom2_name;

            IntervalsPairEntry entry;
            entry.chrom1_id = chromid1;
            entry.chrom2_id = chromid2;
            entry.offset = current_offset;
            entry.length = 0;
            entry.reserved = 0;

            uint64_t bytes_written = 0;
            if (copy_file_contents(pair_file, dat_fp, bytes_written)) {
                entry.length = bytes_written;
                current_offset += bytes_written;
                pair_files_to_remove.push_back(pair_file);
            }

            // Write entry to index (28 bytes total)
            if (fwrite(&entry.chrom1_id, sizeof(entry.chrom1_id), 1, idx_fp) != 1 ||
                fwrite(&entry.chrom2_id, sizeof(entry.chrom2_id), 1, idx_fp) != 1 ||
                fwrite(&entry.offset, sizeof(entry.offset), 1, idx_fp) != 1 ||
                fwrite(&entry.length, sizeof(entry.length), 1, idx_fp) != 1 ||
                fwrite(&entry.reserved, sizeof(entry.reserved), 1, idx_fp) != 1) {
                fclose(dat_fp);
                fclose(idx_fp);
                TGLError("Failed to write index entry for pair %s-%s",
                         chrom1_name.c_str(), chrom2_name.c_str());
            }

            entries.push_back(entry);
        }

        // Compute checksum of entries
        misha::CRC64 crc64;
        uint64_t checksum = crc64.init_incremental();
        for (const auto &entry : entries) {
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.chrom1_id, sizeof(entry.chrom1_id));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.chrom2_id, sizeof(entry.chrom2_id));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.offset, sizeof(entry.offset));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.length, sizeof(entry.length));
        }
        checksum = crc64.finalize_incremental(checksum);

        // Update checksum in index header (seek to checksum position)
        if (fseek(idx_fp, IDX_2D_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError("Failed to seek to checksum position in 2D index");
        }
        if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
            fclose(dat_fp);
            fclose(idx_fp);
            TGLError("Failed to update checksum in index");
        }

        // Flush and sync both files
        fflush(dat_fp);
        fflush(idx_fp);
        fsync(fileno(dat_fp));
        fsync(fileno(idx_fp));

        fclose(dat_fp);
        fclose(idx_fp);

        // Atomic rename: intervals2d.dat first, then intervals2d.idx
        if (rename(dat_path_tmp.c_str(), dat_path.c_str()) != 0) {
            TGLError("Failed to rename %s to %s: %s",
                dat_path_tmp.c_str(), dat_path.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0) {
            TGLError("Failed to rename %s to %s: %s",
                idx_path_tmp.c_str(), idx_path.c_str(), strerror(errno));
        }

        // Remove old per-pair files if requested
        if (remove_old) {
            for (const string &pair_file : pair_files_to_remove) {
                unlink(pair_file.c_str());
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
