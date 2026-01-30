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
#include <tuple>

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

            // Try to find the chromosome file, handling chr prefix mismatch
            string chr_file;
            bool found = false;

            // Try 1: chromosome name as-is
            string candidate = intervset_dir + "/" + chrom_name;
            struct stat st;
            if (stat(candidate.c_str(), &st) == 0) {
                chr_file = candidate;
                found = true;
            } else {
                // Try 2: with "chr" prefix if not already present
                if (chrom_name.substr(0, 3) != "chr") {
                    candidate = intervset_dir + "/chr" + chrom_name;
                    if (stat(candidate.c_str(), &st) == 0) {
                        chr_file = candidate;
                        found = true;
                    }
                } else {
                    // Try 3: without "chr" prefix if present
                    candidate = intervset_dir + "/" + chrom_name.substr(3);
                    if (stat(candidate.c_str(), &st) == 0) {
                        chr_file = candidate;
                        found = true;
                    }
                }
            }

            IntervalsContigEntry entry;
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

        // Validate conversion before removing old files
        // Check that intervals.dat has the expected size
        struct stat dat_stat;
        if (stat(dat_path.c_str(), &dat_stat) != 0) {
            TGLError("Failed to stat %s after conversion", dat_path.c_str());
        }

        if ((uint64_t)dat_stat.st_size != current_offset) {
            TGLError("intervals.dat size mismatch: expected %llu bytes, got %llu bytes",
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

        // Store both chromids and original filenames
        vector<tuple<int, int, string>> pair_files; // (chromid1, chromid2, filename)
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

                // Try to map to chromids (handle chr prefix mismatch)
                int chromid1 = chromkey.chrom2id(chrom1_name.c_str());
                int chromid2 = chromkey.chrom2id(chrom2_name.c_str());

                // If not found, try with/without chr prefix
                if (chromid1 < 0) {
                    if (chrom1_name.substr(0, 3) == "chr") {
                        chromid1 = chromkey.chrom2id(chrom1_name.substr(3).c_str());
                    } else {
                        chromid1 = chromkey.chrom2id(("chr" + chrom1_name).c_str());
                    }
                }
                if (chromid2 < 0) {
                    if (chrom2_name.substr(0, 3) == "chr") {
                        chromid2 = chromkey.chrom2id(chrom2_name.substr(3).c_str());
                    } else {
                        chromid2 = chromkey.chrom2id(("chr" + chrom2_name).c_str());
                    }
                }

                if (chromid1 >= 0 && chromid2 >= 0) {
                    pair_files.push_back(make_tuple(chromid1, chromid2, filename));
                }
            }
        }
        closedir(dir);

        // Sort pairs by (chromid1, chromid2) for stable ordering
        sort(pair_files.begin(), pair_files.end(),
             [](const tuple<int, int, string> &a, const tuple<int, int, string> &b) {
                 if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
                 return get<1>(a) < get<1>(b);
             });

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
        write_2d_index_header(idx_fp, pair_files.size(), 0);

        // Collect entries and concatenate files
        vector<IntervalsPairEntry> entries;
        vector<string> pair_files_to_remove;
        uint64_t current_offset = 0;

        for (const auto &pair_tuple : pair_files) {
            int chromid1 = get<0>(pair_tuple);
            int chromid2 = get<1>(pair_tuple);
            string filename = get<2>(pair_tuple);
            string pair_file = intervset_dir + "/" + filename;

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
                TGLError("Failed to write index entry for pair %s", filename.c_str());
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

        // Validate conversion before removing old files
        // Check that intervals2d.dat has the expected size
        struct stat dat_stat;
        if (stat(dat_path.c_str(), &dat_stat) != 0) {
            TGLError("Failed to stat %s after conversion", dat_path.c_str());
        }

        if ((uint64_t)dat_stat.st_size != current_offset) {
            TGLError("intervals2d.dat size mismatch: expected %llu bytes, got %llu bytes",
                (unsigned long long)current_offset, (unsigned long long)dat_stat.st_size);
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

// Streaming write support for direct indexed format creation
// These functions allow writing indexed format incrementally without creating per-chromosome files first

// Create/open indexed format files for streaming writes (1D intervals)
// Returns an external pointer containing the file handles and state
SEXP gbigintervs_indexed_create(SEXP _intervset_path, SEXP _num_chroms, SEXP _envir) {
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_intervset_path) || Rf_length(_intervset_path) != 1)
            verror("Interval set path must be a string");

        if (!Rf_isInteger(_num_chroms) && !Rf_isReal(_num_chroms))
            verror("Number of chromosomes must be numeric");

        string intervset_path = CHAR(STRING_ELT(_intervset_path, 0));
        int num_chroms = Rf_asInteger(_num_chroms);

        // Create directory if it doesn't exist
        mkdir(intervset_path.c_str(), 0777);

        string dat_path = intervset_path + "/intervals.dat.tmp";
        string idx_path = intervset_path + "/intervals.idx.tmp";

        FILE *dat_fp = fopen(dat_path.c_str(), "wb");
        if (!dat_fp) {
            verror("Cannot create %s: %s", dat_path.c_str(), strerror(errno));
        }

        FILE *idx_fp = fopen(idx_path.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            verror("Cannot create %s: %s", idx_path.c_str(), strerror(errno));
        }

        // Write index header with placeholder checksum
        write_1d_index_header(idx_fp, num_chroms, 0);

        // Return a list with file paths and current offset
        // We'll track state in R and pass it back to subsequent calls
        SEXP result;
        PROTECT(result = Rf_allocVector(VECSXP, 4));

        SEXP names;
        PROTECT(names = Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("dat_path"));
        SET_STRING_ELT(names, 1, Rf_mkChar("idx_path"));
        SET_STRING_ELT(names, 2, Rf_mkChar("intervset_path"));
        SET_STRING_ELT(names, 3, Rf_mkChar("num_chroms"));
        Rf_setAttrib(result, R_NamesSymbol, names);

        SET_VECTOR_ELT(result, 0, Rf_mkString(dat_path.c_str()));
        SET_VECTOR_ELT(result, 1, Rf_mkString(idx_path.c_str()));
        SET_VECTOR_ELT(result, 2, Rf_mkString(intervset_path.c_str()));
        SET_VECTOR_ELT(result, 3, Rf_ScalarInteger(num_chroms));

        fclose(dat_fp);
        fclose(idx_fp);

        UNPROTECT(2);
        return result;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

// Write a chromosome's intervals to the indexed format (1D)
// Returns offset and length for this chromosome
SEXP gbigintervs_indexed_write_chrom(SEXP _dat_path, SEXP _intervals, SEXP _envir) {
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_dat_path) || Rf_length(_dat_path) != 1)
            verror("dat_path must be a string");

        string dat_path = CHAR(STRING_ELT(_dat_path, 0));

        // Open dat file in append mode
        FILE *dat_fp = fopen(dat_path.c_str(), "ab");
        if (!dat_fp) {
            verror("Cannot open %s for appending: %s", dat_path.c_str(), strerror(errno));
        }

        // Get current position (will be the offset for this chromosome)
        long offset = ftell(dat_fp);
        if (offset < 0) {
            fclose(dat_fp);
            verror("Failed to get file position: %s", strerror(errno));
        }

        // Serialize the intervals dataframe
        RSaneSerialize(_intervals, dat_fp);

        // Get new position to calculate length
        long new_pos = ftell(dat_fp);
        if (new_pos < 0) {
            fclose(dat_fp);
            verror("Failed to get file position after write: %s", strerror(errno));
        }

        fclose(dat_fp);

        uint64_t length = new_pos - offset;

        // Return offset and length
        SEXP result;
        PROTECT(result = Rf_allocVector(REALSXP, 2));
        REAL(result)[0] = (double)offset;
        REAL(result)[1] = (double)length;

        SEXP names;
        PROTECT(names = Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(names, 0, Rf_mkChar("offset"));
        SET_STRING_ELT(names, 1, Rf_mkChar("length"));
        Rf_setAttrib(result, R_NamesSymbol, names);

        UNPROTECT(2);
        return result;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

// Finalize the indexed format (1D) - write index entries and compute checksum
SEXP gbigintervs_indexed_finalize(SEXP _idx_path, SEXP _dat_path, SEXP _intervset_path,
                                   SEXP _entries, SEXP _envir) {
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_idx_path) || Rf_length(_idx_path) != 1)
            verror("idx_path must be a string");
        if (!Rf_isString(_dat_path) || Rf_length(_dat_path) != 1)
            verror("dat_path must be a string");
        if (!Rf_isString(_intervset_path) || Rf_length(_intervset_path) != 1)
            verror("intervset_path must be a string");

        string idx_path_tmp = CHAR(STRING_ELT(_idx_path, 0));
        string dat_path_tmp = CHAR(STRING_ELT(_dat_path, 0));
        string intervset_path = CHAR(STRING_ELT(_intervset_path, 0));

        // _entries is a data.frame with columns: chrom_id, offset, length
        if (!Rf_isNewList(_entries))
            verror("entries must be a data.frame");

        SEXP chrom_ids = VECTOR_ELT(_entries, 0);
        SEXP offsets = VECTOR_ELT(_entries, 1);
        SEXP lengths = VECTOR_ELT(_entries, 2);

        int n_entries = Rf_length(chrom_ids);

        // Open idx file for writing entries
        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "r+b");
        if (!idx_fp) {
            verror("Cannot open %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Seek past header to write entries
        // Header size: 8 (magic) + 4 (version) + 4 (num_entries) + 8 (flags) + 8 (checksum) + 4 (reserved) = 36 bytes
        if (fseek(idx_fp, 36, SEEK_SET) != 0) {
            fclose(idx_fp);
            verror("Failed to seek in index file: %s", strerror(errno));
        }

        // Collect entries for checksum computation and write them
        misha::CRC64 crc64;
        uint64_t checksum = crc64.init_incremental();

        for (int i = 0; i < n_entries; i++) {
            IntervalsContigEntry entry;
            entry.chrom_id = INTEGER(chrom_ids)[i];
            entry.offset = (uint64_t)REAL(offsets)[i];
            entry.length = (uint64_t)REAL(lengths)[i];
            entry.reserved = 0;

            // Write entry
            if (fwrite(&entry.chrom_id, sizeof(entry.chrom_id), 1, idx_fp) != 1 ||
                fwrite(&entry.offset, sizeof(entry.offset), 1, idx_fp) != 1 ||
                fwrite(&entry.length, sizeof(entry.length), 1, idx_fp) != 1 ||
                fwrite(&entry.reserved, sizeof(entry.reserved), 1, idx_fp) != 1) {
                fclose(idx_fp);
                verror("Failed to write index entry");
            }

            // Update checksum
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.chrom_id, sizeof(entry.chrom_id));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.offset, sizeof(entry.offset));
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)&entry.length, sizeof(entry.length));
        }

        checksum = crc64.finalize_incremental(checksum);

        // Write checksum at correct position in header
        if (fseek(idx_fp, IDX_1D_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
            fclose(idx_fp);
            verror("Failed to seek to checksum position");
        }
        if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
            fclose(idx_fp);
            verror("Failed to write checksum");
        }

        fflush(idx_fp);
        fsync(fileno(idx_fp));
        fclose(idx_fp);

        // Atomic rename
        string dat_path_final = intervset_path + "/intervals.dat";
        string idx_path_final = intervset_path + "/intervals.idx";

        if (rename(dat_path_tmp.c_str(), dat_path_final.c_str()) != 0) {
            verror("Failed to rename %s to %s: %s",
                   dat_path_tmp.c_str(), dat_path_final.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path_final.c_str()) != 0) {
            verror("Failed to rename %s to %s: %s",
                   idx_path_tmp.c_str(), idx_path_final.c_str(), strerror(errno));
        }

        return R_NilValue;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

// Similar functions for 2D intervals
SEXP gbigintervs_2d_indexed_create(SEXP _intervset_path, SEXP _num_pairs, SEXP _envir) {
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_intervset_path) || Rf_length(_intervset_path) != 1)
            verror("Interval set path must be a string");

        if (!Rf_isInteger(_num_pairs) && !Rf_isReal(_num_pairs))
            verror("Number of pairs must be numeric");

        string intervset_path = CHAR(STRING_ELT(_intervset_path, 0));
        int num_pairs = Rf_asInteger(_num_pairs);

        // Create directory if it doesn't exist
        mkdir(intervset_path.c_str(), 0777);

        string dat_path = intervset_path + "/intervals2d.dat.tmp";
        string idx_path = intervset_path + "/intervals2d.idx.tmp";

        FILE *dat_fp = fopen(dat_path.c_str(), "wb");
        if (!dat_fp) {
            verror("Cannot create %s: %s", dat_path.c_str(), strerror(errno));
        }

        FILE *idx_fp = fopen(idx_path.c_str(), "wb");
        if (!idx_fp) {
            fclose(dat_fp);
            verror("Cannot create %s: %s", idx_path.c_str(), strerror(errno));
        }

        // Write index header with placeholder checksum
        write_2d_index_header(idx_fp, num_pairs, 0);

        fclose(dat_fp);
        fclose(idx_fp);

        // Return paths
        SEXP result;
        PROTECT(result = Rf_allocVector(VECSXP, 3));

        SEXP names;
        PROTECT(names = Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(names, 0, Rf_mkChar("dat_path"));
        SET_STRING_ELT(names, 1, Rf_mkChar("idx_path"));
        SET_STRING_ELT(names, 2, Rf_mkChar("intervset_path"));
        Rf_setAttrib(result, R_NamesSymbol, names);

        SET_VECTOR_ELT(result, 0, Rf_mkString(dat_path.c_str()));
        SET_VECTOR_ELT(result, 1, Rf_mkString(idx_path.c_str()));
        SET_VECTOR_ELT(result, 2, Rf_mkString(intervset_path.c_str()));

        UNPROTECT(2);
        return result;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

// Finalize 2D indexed format
SEXP gbigintervs_2d_indexed_finalize(SEXP _idx_path, SEXP _dat_path, SEXP _intervset_path,
                                      SEXP _entries, SEXP _envir) {
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_idx_path) || Rf_length(_idx_path) != 1)
            verror("idx_path must be a string");
        if (!Rf_isString(_dat_path) || Rf_length(_dat_path) != 1)
            verror("dat_path must be a string");
        if (!Rf_isString(_intervset_path) || Rf_length(_intervset_path) != 1)
            verror("intervset_path must be a string");

        string idx_path_tmp = CHAR(STRING_ELT(_idx_path, 0));
        string dat_path_tmp = CHAR(STRING_ELT(_dat_path, 0));
        string intervset_path = CHAR(STRING_ELT(_intervset_path, 0));

        // _entries is a data.frame with columns: chrom_id1, chrom_id2, offset, length
        if (!Rf_isNewList(_entries))
            verror("entries must be a data.frame");

        SEXP chrom_ids1 = VECTOR_ELT(_entries, 0);
        SEXP chrom_ids2 = VECTOR_ELT(_entries, 1);
        SEXP offsets = VECTOR_ELT(_entries, 2);
        SEXP lengths = VECTOR_ELT(_entries, 3);

        int n_entries = Rf_length(chrom_ids1);

        // Open idx file for writing entries
        FILE *idx_fp = fopen(idx_path_tmp.c_str(), "r+b");
        if (!idx_fp) {
            verror("Cannot open %s: %s", idx_path_tmp.c_str(), strerror(errno));
        }

        // Seek past header
        if (fseek(idx_fp, 36, SEEK_SET) != 0) {
            fclose(idx_fp);
            verror("Failed to seek in index file: %s", strerror(errno));
        }

        // Collect entries for checksum computation
        misha::CRC64 crc64;
        uint64_t checksum = crc64.init_incremental();

        for (int i = 0; i < n_entries; i++) {
            IntervalsPairEntry entry;
            entry.chrom1_id = INTEGER(chrom_ids1)[i];
            entry.chrom2_id = INTEGER(chrom_ids2)[i];
            entry.offset = (uint64_t)REAL(offsets)[i];
            entry.length = (uint64_t)REAL(lengths)[i];

            // Write entry
            if (fwrite(&entry.chrom1_id, sizeof(entry.chrom1_id), 1, idx_fp) != 1 ||
                fwrite(&entry.chrom2_id, sizeof(entry.chrom2_id), 1, idx_fp) != 1 ||
                fwrite(&entry.offset, sizeof(entry.offset), 1, idx_fp) != 1 ||
                fwrite(&entry.length, sizeof(entry.length), 1, idx_fp) != 1) {
                fclose(idx_fp);
                verror("Failed to write index entry");
            }

            // Update checksum
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

        // Write checksum
        if (fseek(idx_fp, IDX_2D_HEADER_SIZE_TO_CHECKSUM, SEEK_SET) != 0) {
            fclose(idx_fp);
            verror("Failed to seek to checksum position");
        }
        if (fwrite(&checksum, sizeof(checksum), 1, idx_fp) != 1) {
            fclose(idx_fp);
            verror("Failed to write checksum");
        }

        fflush(idx_fp);
        fsync(fileno(idx_fp));
        fclose(idx_fp);

        // Atomic rename
        string dat_path_final = intervset_path + "/intervals2d.dat";
        string idx_path_final = intervset_path + "/intervals2d.idx";

        if (rename(dat_path_tmp.c_str(), dat_path_final.c_str()) != 0) {
            verror("Failed to rename %s to %s: %s",
                   dat_path_tmp.c_str(), dat_path_final.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path_final.c_str()) != 0) {
            verror("Failed to rename %s to %s: %s",
                   idx_path_tmp.c_str(), idx_path_final.c_str(), strerror(errno));
        }

        return R_NilValue;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
