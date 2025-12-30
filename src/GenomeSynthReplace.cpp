/*
 * GenomeSynthReplace.cpp
 *
 * C++ implementation for iteratively replacing k-mers in genome sequences.
 */

#include <algorithm>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

/**
 * Helper: Convert character to uppercase
 */
static inline char to_upper(char c) {
    return (c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c;
}

/**
 * Helper: Iteratively replace target k-mer in a sequence
 *
 * Scans through the sequence and replaces all occurrences of target with replacement.
 * If a replacement creates a new target instance, it will be replaced in subsequent passes.
 * Continues until no more target instances are found.
 *
 * @param seq The sequence to modify (in-place)
 * @param target The k-mer to search for
 * @param replacement The replacement sequence
 * @return Number of replacements made
 */
static size_t replace_kmer_iterative(vector<char>& seq, const string& target, const string& replacement) {
    if (target.length() != replacement.length()) {
        // Lengths must match for in-place replacement
        return 0;
    }

    size_t kmer_len = target.length();
    size_t total_replacements = 0;
    bool found_any = true;

    // Convert target and replacement to uppercase for comparison
    string target_upper = target;
    string replacement_upper = replacement;
    for (size_t i = 0; i < target_upper.length(); ++i) {
        target_upper[i] = to_upper(target_upper[i]);
        replacement_upper[i] = to_upper(replacement_upper[i]);
    }

    // Iteratively replace until no more instances are found
    while (found_any) {
        found_any = false;

        // Scan through sequence
        for (size_t i = 0; i + kmer_len <= seq.size(); ) {
            // Check if we have a match
            bool match = true;
            for (size_t j = 0; j < kmer_len; ++j) {
                if (to_upper(seq[i + j]) != target_upper[j]) {
                    match = false;
                    break;
                }
            }

            if (match) {
                // Replace this k-mer
                for (size_t j = 0; j < kmer_len; ++j) {
                    seq[i + j] = replacement_upper[j];
                }
                found_any = true;
                total_replacements++;

                // Move forward by 1 to check if we created a new instance
                // at position i+1 (bubble sort behavior)
                i++;
            } else {
                // No match, move to next position
                i++;
            }
        }
    }

    return total_replacements;
}

/**
 * Helper: Write sequence to FASTA format.
 */
static void write_fasta(ofstream& ofs, const string& chrom_name,
                        const vector<char>& seq, int line_width = 60) {
    ofs << ">" << chrom_name << "\n";
    for (size_t i = 0; i < seq.size(); i += line_width) {
        size_t len = min(static_cast<size_t>(line_width), seq.size() - i);
        ofs.write(&seq[i], len);
        ofs << "\n";
    }
}

/**
 * Helper: Write sequence to .seq binary format.
 */
static void write_seq(BufferedFile& bfile, const vector<char>& seq) {
    bfile.write(&seq[0], seq.size());
}

extern "C" {

/**
 * C_gsynth_replace_kmer: Iteratively replace a k-mer in genome sequences.
 *
 * @param _target Target k-mer to replace
 * @param _replacement Replacement sequence
 * @param _intervals R intervals object for regions to process
 * @param _output_path Output file path (ignored if output_format = 2)
 * @param _output_format Integer: 0 = misha .seq, 1 = FASTA, 2 = return vector
 * @param _envir R environment
 *
 * @return R_NilValue for file output, or character vector for output_format = 2
 */
SEXP C_gsynth_replace_kmer(SEXP _target, SEXP _replacement, SEXP _intervals,
                            SEXP _output_path, SEXP _output_format, SEXP _envir) {
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Extract parameters
        const char* target_cstr = CHAR(STRING_ELT(_target, 0));
        const char* replacement_cstr = CHAR(STRING_ELT(_replacement, 0));
        string target(target_cstr);
        string replacement(replacement_cstr);

        if (target.empty() || replacement.empty()) {
            verror("target and replacement cannot be empty");
        }

        if (target.length() != replacement.length()) {
            verror("target and replacement must have the same length");
        }

        const char* output_path = CHAR(STRING_ELT(_output_path, 0));
        int output_format = Rf_asInteger(_output_format);

        // Get chromosome key and setup
        const GenomeChromKey& chromkey = iu.get_chromkey();
        int num_chroms = chromkey.get_num_chroms();

        // Parse intervals
        vector<vector<GInterval>> intervals_per_chrom;
        intervals_per_chrom.resize(num_chroms);

        if (!Rf_isNull(_intervals)) {
            GIntervalsFetcher1D* intervals_fetcher = NULL;
            iu.convert_rintervs(_intervals, &intervals_fetcher, NULL);
            unique_ptr<GIntervalsFetcher1D> intervals_guard(intervals_fetcher);
            intervals_fetcher->sort();

            for (intervals_fetcher->begin_iter(); !intervals_fetcher->isend();
                 intervals_fetcher->next()) {
                const GInterval& iv = intervals_fetcher->cur_interval();
                if (iv.chromid >= 0 && iv.chromid < num_chroms) {
                    intervals_per_chrom[iv.chromid].push_back(iv);
                }
            }
        }

        // Vector to collect sequences when output_format = 2 (vector mode)
        vector<string> collected_seqs;

        // Compute total range for progress reporting
        uint64_t total_range = 0;
        for (int c = 0; c < num_chroms; ++c) {
            for (const auto& iv : intervals_per_chrom[c]) {
                total_range += iv.end - iv.start;
            }
        }

        Progress_reporter progress;
        progress.init(total_range, 1000000);

        // Set up sequence fetcher
        GenomeSeqFetch seqfetch;
        seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

        // Open output file (only if not vector mode)
        ofstream fasta_ofs;
        BufferedFile seq_bfile;
        if (output_format == 1) {
            // FASTA
            fasta_ofs.open(output_path);
            if (!fasta_ofs) {
                verror("Failed to open output file: %s", output_path);
            }
        } else if (output_format == 0) {
            // .seq binary
            seq_bfile.open(output_path, "wb");
            if (seq_bfile.error()) {
                verror("Failed to open output file: %s", output_path);
            }
        }
        // output_format == 2: vector mode, no file to open

        size_t total_replacements = 0;

        // Process each chromosome's intervals
        for (int chromid = 0; chromid < num_chroms; ++chromid) {
            const vector<GInterval>& intervals = intervals_per_chrom[chromid];
            if (intervals.empty()) {
                continue;
            }

            int64_t chrom_size = chromkey.get_chrom_size(chromid);
            if (chrom_size <= 0) continue;

            const string& chrom_name = chromkey.id2chrom(chromid);

            for (size_t iv_idx = 0; iv_idx < intervals.size(); ++iv_idx) {
                const GInterval& iv = intervals[iv_idx];
                int64_t interval_start = max<int64_t>(0, iv.start);
                int64_t interval_end = min<int64_t>(chrom_size, iv.end);

                if (interval_end <= interval_start) {
                    continue;
                }

                // Read original sequence
                vector<char> seq;
                GInterval read_interval(chromid, interval_start, interval_end, 0);
                seqfetch.read_interval(read_interval, chromkey, seq);

                // Perform iterative replacement
                size_t num_replaced = replace_kmer_iterative(seq, target, replacement);
                total_replacements += num_replaced;

                // Write output based on format
                if (output_format == 2) {
                    // Vector mode: collect sequence
                    string seq_str(seq.begin(), seq.end());
                    collected_seqs.push_back(seq_str);
                } else if (output_format == 1) {
                    // FASTA mode
                    string header = chrom_name;
                    if (intervals.size() > 1) {
                        header += "_" + to_string(interval_start) + "_" + to_string(interval_end);
                    }
                    write_fasta(fasta_ofs, header, seq);
                } else if (output_format == 0) {
                    // Binary .seq mode
                    write_seq(seq_bfile, seq);
                }

                progress.report(interval_end - interval_start);
            }
        }

        progress.report_last();

        // Close output files
        if (output_format == 1) {
            fasta_ofs.close();
        } else if (output_format == 0) {
            seq_bfile.close();
        }

        // Return result based on output format
        if (output_format == 2) {
            // Return character vector
            SEXP result;
            PROTECT(result = Rf_allocVector(STRSXP, collected_seqs.size()));
            for (size_t i = 0; i < collected_seqs.size(); ++i) {
                SET_STRING_ELT(result, i, Rf_mkChar(collected_seqs[i].c_str()));
            }
            UNPROTECT(1);
            return result;
        }

        return R_NilValue;

    } catch (const TGLException& e) {
        rerror("%s", e.msg());
    } catch (const exception& e) {
        rerror("Error in C_gsynth_replace_kmer: %s", e.what());
    }
    return R_NilValue;
}

}  // extern "C"
