/*
 * GenomeKmerDist.cpp
 *
 * C++ implementation for counting k-mer distribution in genomic sequences.
 */

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"
#include "GIntervalsBigSet1D.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "StratifiedMarkovModel.h"

using namespace std;
using namespace rdb;

/**
 * Helper: Check if a position is within any masked interval.
 * Uses a cursor for efficient sequential access.
 */
static bool is_masked_kmer(int64_t pos, const vector<GInterval>& mask_intervals,
                           size_t& cursor) {
    while (cursor < mask_intervals.size() &&
           mask_intervals[cursor].end <= pos) {
        ++cursor;
    }

    if (cursor < mask_intervals.size() &&
        mask_intervals[cursor].start <= pos &&
        pos < mask_intervals[cursor].end) {
        return true;
    }

    return false;
}

/**
 * Encode a k-mer string to an integer index.
 * Returns -1 if the k-mer contains invalid characters.
 */
static int encode_kmer(const char* seq, int k) {
    int idx = 0;
    for (int i = 0; i < k; ++i) {
        int base = StratifiedMarkovModel::encode_base(seq[i]);
        if (base < 0) return -1;
        idx = (idx << 2) | base;
    }
    return idx;
}

/**
 * Decode a k-mer index to a string.
 */
static string decode_kmer(int idx, int k) {
    string result(k, 'N');
    for (int i = k - 1; i >= 0; --i) {
        result[i] = StratifiedMarkovModel::decode_base(idx & 3);
        idx >>= 2;
    }
    return result;
}

extern "C" {

/**
 * C_gseq_kmer_dist: Count k-mer distribution in genomic intervals.
 *
 * @param _intervals R intervals object
 * @param _k Integer k-mer size (1-10)
 * @param _mask R intervals object for mask regions (NULL if no mask)
 * @param _envir R environment
 *
 * @return A data frame with kmer and count columns
 */
SEXP C_gseq_kmer_dist(SEXP _intervals, SEXP _k, SEXP _mask, SEXP _envir) {
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Get k-mer size
        int k = Rf_asInteger(_k);
        if (k < 1 || k > 10) {
            verror("k must be between 1 and 10");
        }

        int num_kmers = 1 << (2 * k);  // 4^k

        // Parse intervals
        GIntervalsFetcher1D* intervals = NULL;
        iu.convert_rintervs(_intervals, &intervals, NULL);
        unique_ptr<GIntervalsFetcher1D> intervals_guard(intervals);
        intervals->sort();

        // Compute total range for progress
        uint64_t total_range = 0;
        for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
            total_range += intervals->cur_interval().end -
                          intervals->cur_interval().start;
        }

        Progress_reporter progress;
        progress.init(total_range, 1000000);

        // Parse mask intervals if provided
        const GenomeChromKey& chromkey = iu.get_chromkey();
        int num_chroms = chromkey.get_num_chroms();
        vector<vector<GInterval>> mask_per_chrom(num_chroms);

        if (!Rf_isNull(_mask)) {
            GIntervalsFetcher1D* mask_intervals = NULL;
            iu.convert_rintervs(_mask, &mask_intervals, NULL);
            unique_ptr<GIntervalsFetcher1D> mask_guard(mask_intervals);
            mask_intervals->sort();

            for (mask_intervals->begin_iter(); !mask_intervals->isend();
                 mask_intervals->next()) {
                const GInterval& iv = mask_intervals->cur_interval();
                if (iv.chromid >= 0 && iv.chromid < num_chroms) {
                    mask_per_chrom[iv.chromid].push_back(iv);
                }
            }
        }

        // Set up sequence fetcher
        GenomeSeqFetch seqfetch;
        seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

        // Initialize counts
        vector<uint64_t> counts(num_kmers, 0);
        uint64_t total_valid = 0;
        uint64_t total_masked = 0;
        uint64_t total_n = 0;

        // Process each interval
        for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
            const GInterval& iv = intervals->cur_interval();
            int chromid = iv.chromid;

            // Load sequence
            vector<char> seq;
            seqfetch.read_interval(iv, chromkey, seq);

            if ((int)seq.size() < k) {
                progress.report(iv.end - iv.start);
                continue;
            }

            const vector<GInterval>& mask_ivs = mask_per_chrom[chromid];
            size_t mask_cursor = 0;

            // Skip mask intervals before this interval
            while (mask_cursor < mask_ivs.size() &&
                   mask_ivs[mask_cursor].end <= iv.start) {
                ++mask_cursor;
            }

            // Scan with k-mer sliding window
            int64_t seq_size = static_cast<int64_t>(seq.size());
            for (int64_t pos = 0; pos <= seq_size - k; ++pos) {
                int64_t genome_pos = iv.start + pos;

                // Check if masked
                if (is_masked_kmer(genome_pos, mask_ivs, mask_cursor)) {
                    ++total_masked;
                    continue;
                }

                // Encode k-mer
                int kmer_idx = encode_kmer(&seq[pos], k);
                if (kmer_idx < 0) {
                    ++total_n;
                    continue;
                }

                ++counts[kmer_idx];
                ++total_valid;
            }

            progress.report(iv.end - iv.start);
            check_interrupt();
        }

        progress.report_last();

        // Build result data frame
        SEXP answer, kmer_col, count_col, names, row_names;

        // Count non-zero k-mers for efficient storage
        int num_nonzero = 0;
        for (int i = 0; i < num_kmers; ++i) {
            if (counts[i] > 0) ++num_nonzero;
        }

        rprotect(kmer_col = Rf_allocVector(STRSXP, num_nonzero));
        rprotect(count_col = Rf_allocVector(REALSXP, num_nonzero));
        rprotect(row_names = Rf_allocVector(INTSXP, num_nonzero));

        int idx = 0;
        for (int i = 0; i < num_kmers; ++i) {
            if (counts[i] > 0) {
                SET_STRING_ELT(kmer_col, idx, Rf_mkChar(decode_kmer(i, k).c_str()));
                REAL(count_col)[idx] = static_cast<double>(counts[i]);
                INTEGER(row_names)[idx] = idx + 1;
                ++idx;
            }
        }

        // Create data frame
        rprotect(answer = Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(answer, 0, kmer_col);
        SET_VECTOR_ELT(answer, 1, count_col);

        rprotect(names = Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(names, 0, Rf_mkChar("kmer"));
        SET_STRING_ELT(names, 1, Rf_mkChar("count"));
        Rf_setAttrib(answer, R_NamesSymbol, names);

        Rf_setAttrib(answer, R_RowNamesSymbol, row_names);
        Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));

        return answer;

    } catch (TGLException& e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc& e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

}  // extern "C"
