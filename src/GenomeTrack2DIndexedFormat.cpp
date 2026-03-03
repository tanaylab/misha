/*
 * GenomeTrack2DIndexedFormat.cpp
 *
 * Converts per-chromosome-pair 2D tracks (RECTS/POINTS) to indexed format
 * (track.dat + track.idx). Uses raw byte copy of per-pair files.
 */

#include <cstdint>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>

#include "GenomeTrack.h"
#include "TrackIndex2D.h"
#include "CRC64.h"
#include "BufferedFile.h"
#include "TGLException.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// Helper function to copy file contents (raw byte copy)
static bool copy_file_contents_2d(const string &src, FILE *dest, uint64_t &bytes_written) {
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

SEXP gtrack2d_convert_to_indexed(SEXP _track, SEXP _remove_old, SEXP _envir) {
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

        // Get chromosome key for resolving chrom names
        IntervUtils iu(_envir);
        const GenomeChromKey &chromkey = iu.get_chromkey();

        // Determine track type (must be RECTS or POINTS)
        GenomeTrack::Type track_type_enum = GenomeTrack::get_type(track_dir.c_str(), chromkey);
        MishaTrack2DType track_type_2d = MishaTrack2DType::RECTS; // initialized to suppress warning

        switch (track_type_enum) {
            case GenomeTrack::RECTS:
                track_type_2d = MishaTrack2DType::RECTS;
                break;
            case GenomeTrack::POINTS:
                track_type_2d = MishaTrack2DType::POINTS;
                break;
            default:
                verror("Only 2D tracks (rectangles, points) can be converted with this function");
        }

        // Check if already in indexed format
        string idx_path = track_dir + "/track.idx";
        string dat_path = track_dir + "/track.dat";
        struct stat st;
        if (stat(idx_path.c_str(), &st) == 0) {
            // Already indexed. Before removing, check if per-pair files exist.
            // If not (removed by previous remove.old=TRUE), extract them from track.dat.
            try {
                auto old_idx = TrackIndex2D::get_track_index_2d(track_dir);
                if (old_idx && old_idx->is_loaded()) {
                    // Check if any per-pair files exist
                    bool have_pair_files = false;
                    for (const auto &entry : old_idx->entries()) {
                        string pair_file = track_dir + "/" + GenomeTrack::get_2d_filename(chromkey, entry.chrom1_id, entry.chrom2_id);
                        struct stat pst;
                        if (stat(pair_file.c_str(), &pst) == 0) {
                            have_pair_files = true;
                            break;
                        }
                    }

                    if (!have_pair_files) {
                        // Extract per-pair data from track.dat back to per-pair files
                        FILE *dat_fp = fopen(dat_path.c_str(), "rb");
                        if (dat_fp) {
                            for (const auto &entry : old_idx->entries()) {
                                string pair_file = track_dir + "/" + GenomeTrack::get_2d_filename(chromkey, entry.chrom1_id, entry.chrom2_id);
                                vector<char> buf(entry.length);
                                if (fseek(dat_fp, entry.offset, SEEK_SET) == 0 &&
                                    fread(buf.data(), 1, entry.length, dat_fp) == entry.length) {
                                    FILE *pf = fopen(pair_file.c_str(), "wb");
                                    if (pf) {
                                        fwrite(buf.data(), 1, entry.length, pf);
                                        fclose(pf);
                                    }
                                }
                            }
                            fclose(dat_fp);
                        }
                    }
                }
            } catch (...) {
                // If index is corrupted, just proceed - per-pair files might still exist
            }

            // Now remove the old index files
            unlink(idx_path.c_str());
            unlink(dat_path.c_str());
            // Clear the index cache so stale entries are not used
            TrackIndex2D::clear_cache();
        }

        // Enumerate existing per-pair files in the track directory
        DIR *dir = opendir(track_dir.c_str());
        if (!dir) {
            verror("Cannot open track directory %s: %s",
                   track_dir.c_str(), strerror(errno));
        }

        // Collect (chromid1, chromid2, filename) tuples for all valid pair files
        vector<tuple<int, int, string>> pair_files;
        struct dirent *dentry;
        while ((dentry = readdir(dir)) != nullptr) {
            string filename = dentry->d_name;

            // Skip . and .. and hidden files
            if (filename.empty() || filename[0] == '.') continue;

            // Skip track.dat, track.idx, and any .tmp files
            if (filename == "track.dat" || filename == "track.idx") continue;
            if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".tmp") continue;

            // Look for files matching the chrom1-chrom2 pattern
            size_t dash_pos = filename.find('-');
            if (dash_pos == string::npos || dash_pos == 0 || dash_pos >= filename.length() - 1) {
                continue;
            }

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
        closedir(dir);

        if (pair_files.empty()) {
            verror("No valid chromosome pair files found in track directory %s", track_dir.c_str());
        }

        // Sort pairs by (chromid1, chromid2) for deterministic output
        sort(pair_files.begin(), pair_files.end(),
             [](const tuple<int, int, string> &a, const tuple<int, int, string> &b) {
                 if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
                 return get<1>(a) < get<1>(b);
             });

        // Prepare temporary file paths
        dat_path_tmp = track_dir + "/track.dat.tmp";
        idx_path_tmp = track_dir + "/track.idx.tmp";

        // Open temporary data file
        FILE *dat_fp = fopen(dat_path_tmp.c_str(), "wb");
        if (!dat_fp) {
            TGLError<GenomeTrack>("Failed to create %s: %s", dat_path_tmp.c_str(), strerror(errno));
        }

        // Collect index entries and concatenate pair file contents
        vector<Track2DPairEntry> entries;
        vector<string> pair_files_to_remove;
        uint64_t current_offset = 0;

        for (const auto &pair_tuple : pair_files) {
            int chromid1 = get<0>(pair_tuple);
            int chromid2 = get<1>(pair_tuple);
            string filename = get<2>(pair_tuple);
            string pair_file = track_dir + "/" + filename;

            Track2DPairEntry entry(chromid1, chromid2, current_offset, 0);

            uint64_t bytes_written = 0;
            if (copy_file_contents_2d(pair_file, dat_fp, bytes_written)) {
                entry.length = bytes_written;
                current_offset += bytes_written;
                pair_files_to_remove.push_back(pair_file);
            }

            entries.push_back(entry);
        }

        // Flush and sync data file
        fflush(dat_fp);
        fsync(fileno(dat_fp));
        fclose(dat_fp);

        // Write index file using TrackIndex2D::write_index()
        TrackIndex2D::write_index(idx_path_tmp, track_type_2d, entries);

        // Atomic rename: track.dat first, then track.idx
        if (rename(dat_path_tmp.c_str(), dat_path.c_str()) != 0) {
            TGLError<GenomeTrack>("Failed to rename %s to %s: %s",
                dat_path_tmp.c_str(), dat_path.c_str(), strerror(errno));
        }

        if (rename(idx_path_tmp.c_str(), idx_path.c_str()) != 0) {
            TGLError<GenomeTrack>("Failed to rename %s to %s: %s",
                idx_path_tmp.c_str(), idx_path.c_str(), strerror(errno));
        }

        // Validate conversion: check that track.dat has the expected size
        struct stat dat_stat;
        if (stat(dat_path.c_str(), &dat_stat) != 0) {
            TGLError<GenomeTrack>("Failed to stat %s after conversion", dat_path.c_str());
        }

        if ((uint64_t)dat_stat.st_size != current_offset) {
            TGLError<GenomeTrack>("track.dat size mismatch: expected %llu bytes, got %llu bytes",
                (unsigned long long)current_offset, (unsigned long long)dat_stat.st_size);
        }

        // Validate by loading the index and checking entries
        TrackIndex2D::clear_cache();
        auto loaded_idx = TrackIndex2D::get_track_index_2d(track_dir);
        if (!loaded_idx || !loaded_idx->is_loaded()) {
            TGLError<GenomeTrack>("Failed to validate 2D track index after conversion for %s", track_dir.c_str());
        }

        if (loaded_idx->num_entries() != entries.size()) {
            TGLError<GenomeTrack>("2D track index entry count mismatch: expected %u, got %u",
                (unsigned)entries.size(), (unsigned)loaded_idx->num_entries());
        }

        // Verify each entry matches
        for (const auto &entry : entries) {
            const Track2DPairEntry *loaded_entry = loaded_idx->get_entry(entry.chrom1_id, entry.chrom2_id);
            if (!loaded_entry) {
                TGLError<GenomeTrack>("2D track index missing pair (%u,%u) after conversion",
                    entry.chrom1_id, entry.chrom2_id);
            }
            if (loaded_entry->offset != entry.offset || loaded_entry->length != entry.length) {
                TGLError<GenomeTrack>("2D track index offset/length mismatch for pair (%u,%u): "
                    "expected offset=%llu length=%llu, got offset=%llu length=%llu",
                    entry.chrom1_id, entry.chrom2_id,
                    (unsigned long long)entry.offset, (unsigned long long)entry.length,
                    (unsigned long long)loaded_entry->offset, (unsigned long long)loaded_entry->length);
            }
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
        if (!dat_path_tmp.empty()) unlink(dat_path_tmp.c_str());
        if (!idx_path_tmp.empty()) unlink(idx_path_tmp.c_str());
        verror("%s", e.msg());
    } catch (const bad_alloc &e) {
        // Clean up temporary files on error
        if (!dat_path_tmp.empty()) unlink(dat_path_tmp.c_str());
        if (!idx_path_tmp.empty()) unlink(idx_path_tmp.c_str());
        verror("Out of memory");
    }

    return R_NilValue;
}

} // extern "C"
