/*
 * GenomeSynthSample.cpp
 *
 * C++ implementation for sampling a synthetic genome from a trained
 * stratified Markov-k model.
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>

#include <R_ext/Random.h>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"
#include "GIntervalsBigSet1D.h"
#include "MaskUtils.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "StratifiedMarkovModel.h"

using namespace std;
using namespace rdb;

/**
 * One entry of a FASTA index (.fai) file:
 *   <name>\t<length>\t<offset>\t<linebases>\t<linewidth>
 * where offset is the byte offset of the first sequence base (i.e., just
 * after the header line) and linewidth = linebases + 1 for LF line endings.
 */
struct FaiEntry {
    string    name;
    long long length;
    long long offset;
    int       linebases;
    int       linewidth;
};

/**
 * Helper: Write sequence to FASTA format and record a FaiEntry describing
 * the written record. The caller passes the byte offset BEFORE the header
 * write (from fasta_ofs.tellp()); we use it to compute the first-base byte
 * offset that samtools expects. line_width is the FASTA line wrap width.
 */
static void write_fasta(ofstream& ofs, const string& chrom_name,
                        const vector<char>& seq, int line_width,
                        long long header_byte_offset,
                        FaiEntry* out_entry) {
    ofs << ">" << chrom_name << "\n";
    for (size_t i = 0; i < seq.size(); i += line_width) {
        size_t len = min(static_cast<size_t>(line_width), seq.size() - i);
        ofs.write(&seq[i], len);
        ofs << "\n";
    }
    if (out_entry) {
        // First-base offset = byte position of header + ">" + name + "\n"
        out_entry->name      = chrom_name;
        out_entry->length    = static_cast<long long>(seq.size());
        out_entry->offset    = header_byte_offset + 1 +
                               static_cast<long long>(chrom_name.size()) + 1;
        out_entry->linebases = (out_entry->length > 0)
                                   ? (int)min((long long)line_width, out_entry->length)
                                   : 0;
        out_entry->linewidth = (out_entry->length > 0) ? out_entry->linebases + 1 : 0;
    }
}

// Overload kept for internal call sites that don't care about .fai.
static inline void write_fasta(ofstream& ofs, const string& chrom_name,
                                const vector<char>& seq, int line_width = 60) {
    write_fasta(ofs, chrom_name, seq, line_width, 0, nullptr);
}

/**
 * Helper: Write sequence to .seq binary format.
 */
static void write_seq(BufferedFile& bfile, const vector<char>& seq) {
    bfile.write(&seq[0], seq.size());
}

extern "C" {

/**
 * C_gsynth_sample: Sample a synthetic genome from a trained Markov model.
 *
 * @param _cdf_list List of CDF matrices (one per bin)
 * @param _breaks Numeric vector of bin boundaries
 * @param _bin_indices Integer vector of bin indices for each position
 * @param _iter_starts Integer vector of iterator interval start positions
 * @param _iter_chroms Integer vector of iterator interval chromosome IDs
 * @param _intervals R intervals object for output regions
 * @param _mask_copy R intervals object for regions to copy from original (NULL if none)
 * @param _output_path Output file path (ignored if output_format = 2)
 * @param _output_format Integer: 0 = misha .seq, 1 = FASTA, 2 = return vector
 * @param _n_samples Integer: number of samples to generate per interval
 * @param _k Integer: Markov order k
 * @param _iter_size Integer: iterator bin size in bp (e.g., model$iterator)
 * @param _envir R environment
 *
 * @return R_NilValue on success for file output, or character vector for
 *         output_format = 2
 */
SEXP C_gsynth_sample(SEXP _cdf_list, SEXP _breaks, SEXP _bin_indices,
                      SEXP _iter_starts, SEXP _iter_chroms, SEXP _intervals,
                      SEXP _mask_copy, SEXP _output_path,
                      SEXP _output_format, SEXP _n_samples, SEXP _k,
                      SEXP _iter_size, SEXP _envir) {
    try {
        struct RNGStateGuard {
            bool active = false;
            void acquire() {
                GetRNGstate();
                active = true;
            }
            void release() {
                if (active) {
                    PutRNGstate();
                    active = false;
                }
            }
            ~RNGStateGuard() {
                if (active) {
                    PutRNGstate();
                }
            }
        };

        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        // Get R's RNG state (seed is set via set.seed() in R before calling)
        RNGStateGuard rng_guard;
        rng_guard.acquire();

        // Parse Markov order k (default 5 for backward compatibility)
        int k = Rf_isNull(_k) ? 5 : INTEGER(_k)[0];
        if (k < 1 || k > StratifiedMarkovModel::MAX_K) {
            verror("k must be between 1 and %d", StratifiedMarkovModel::MAX_K);
        }
        int num_kmers = 1;
        for (int i = 0; i < k; ++i) num_kmers *= NUM_BASES;  // 4^k

        // Extract breaks
        int num_breaks = Rf_length(_breaks);
        int num_bins = num_breaks - 1;
        if (num_bins <= 0) {
            verror("breaks must have at least 2 elements");
        }
        double* breaks = REAL(_breaks);
        vector<double> breaks_vec(breaks, breaks + num_breaks);
        vector<vector<float>> cdf_data(num_bins);
        for (int b = 0; b < num_bins; ++b) {
            SEXP cdf_mat = VECTOR_ELT(_cdf_list, b);
            double* cdf_ptr = REAL(cdf_mat);
            cdf_data[b].resize(num_kmers * NUM_BASES);
            for (int ctx = 0; ctx < num_kmers; ++ctx) {
                for (int base = 0; base < NUM_BASES; ++base) {
                    // R matrices are column-major: [row + col * nrow]
                    cdf_data[b][ctx * NUM_BASES + base] =
                        static_cast<float>(cdf_ptr[ctx + base * num_kmers]);
                }
            }
        }

        // Extract bin indices from track extraction
        int num_iter_positions = Rf_length(_bin_indices);
        int* bin_indices = INTEGER(_bin_indices);
        int* iter_starts = INTEGER(_iter_starts);
        int* iter_chroms = INTEGER(_iter_chroms);

        // Iterator bin size (in bp). Passed explicitly from R (e.g. model$iterator)
        // so unaligned intervals are handled correctly (was silently inferred
        // from the first same-chrom diff — see misha#94).
        int iter_size = Rf_asInteger(_iter_size);
        if (iter_size <= 0) {
            verror("iter_size must be a positive integer; got %d", iter_size);
        }

        const GenomeChromKey& chromkey = iu.get_chromkey();
        int num_chroms = chromkey.get_num_chroms();

        // Parse sampling intervals
        vector<vector<GInterval>> sample_per_chrom;
        sample_per_chrom.resize(num_chroms);
        if (!Rf_isNull(_intervals)) {
            GIntervalsFetcher1D* sample_intervals = NULL;
            iu.convert_rintervs(_intervals, &sample_intervals, NULL);
            unique_ptr<GIntervalsFetcher1D> sample_guard(sample_intervals);
            sample_intervals->sort();

            for (sample_intervals->begin_iter(); !sample_intervals->isend();
                 sample_intervals->next()) {
                const GInterval& iv = sample_intervals->cur_interval();
                if (iv.chromid >= 0 && iv.chromid < num_chroms) {
                    sample_per_chrom[iv.chromid].push_back(iv);
                }
            }
        }

        // Parse mask_copy intervals (regions to copy from original genome)
        vector<vector<GInterval>> mask_copy_per_chrom;
        mask_copy_per_chrom.resize(num_chroms);

        if (!Rf_isNull(_mask_copy)) {
            GIntervalsFetcher1D* mask_intervals = NULL;
            iu.convert_rintervs(_mask_copy, &mask_intervals, NULL);
            unique_ptr<GIntervalsFetcher1D> mask_guard(mask_intervals);
            mask_intervals->sort();

            for (mask_intervals->begin_iter(); !mask_intervals->isend();
                 mask_intervals->next()) {
                const GInterval& iv = mask_intervals->cur_interval();
                if (iv.chromid >= 0 &&
                    iv.chromid < (int)mask_copy_per_chrom.size()) {
                    mask_copy_per_chrom[iv.chromid].push_back(iv);
                }
            }
        }

        // Get output settings
        const char* output_path = CHAR(STRING_ELT(_output_path, 0));
        int output_format = Rf_asInteger(_output_format);
        int n_samples = Rf_asInteger(_n_samples);
        if (n_samples < 1) {
            n_samples = 1;
        }

        // Vector to collect sequences when output_format = 2 (vector mode)
        vector<string> collected_seqs;

        // Compute total range for progress reporting
        uint64_t total_range = 0;
        for (int c = 0; c < num_chroms; ++c) {
            for (const auto& iv : sample_per_chrom[c]) {
                total_range += iv.end - iv.start;
            }
        }
        total_range *= n_samples;  // Multiply by number of samples

        Progress_reporter progress;
        progress.init(total_range, 1000000);

        // Set up sequence fetcher (for mask_mode = copy)
        GenomeSeqFetch seqfetch;
        seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

        // Open output file (only if not vector mode)
        ofstream fasta_ofs;
        BufferedFile seq_bfile;
        vector<FaiEntry> fai_entries;  // collected during FASTA write, flushed after
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

        // Build per-chromosome bin lookup
        vector<vector<pair<int64_t, int>>> chrom_bins(num_chroms);
        for (int i = 0; i < num_iter_positions; ++i) {
            int chromid = iter_chroms[i];
            if (chromid >= 0 && chromid < num_chroms) {
                chrom_bins[chromid].push_back({iter_starts[i], bin_indices[i]});
            }
        }
        for (int c = 0; c < num_chroms; ++c) {
            sort(chrom_bins[c].begin(), chrom_bins[c].end());
        }

        // Process each chromosome's requested intervals
        for (int chromid = 0; chromid < num_chroms; ++chromid) {
            const vector<GInterval>& sample_ivs = sample_per_chrom[chromid];
            if (sample_ivs.empty()) {
                continue;
            }
            int64_t chrom_size = chromkey.get_chrom_size(chromid);
            if (chrom_size <= 0) continue;

            const string& chrom_name = chromkey.id2chrom(chromid);
            const vector<GInterval>& mask_copy_ivs = mask_copy_per_chrom[chromid];
            const vector<pair<int64_t, int>>& bins = chrom_bins[chromid];

            for (size_t iv_idx = 0; iv_idx < sample_ivs.size(); ++iv_idx) {
                const GInterval& iv = sample_ivs[iv_idx];
                int64_t interval_start = max<int64_t>(0, iv.start);
                int64_t interval_end = min<int64_t>(chrom_size, iv.end);
                if (interval_end <= interval_start) {
                    continue;
                }

                int64_t interval_len = interval_end - interval_start;

                // Load original sequence for mask_copy regions (only once per interval)
                vector<char> original_seq;
                if (!mask_copy_ivs.empty()) {
                    GInterval interval(chromid, interval_start, interval_end, 0);
                    seqfetch.read_interval(interval, chromkey, original_seq);
                }

                // Generate n_samples samples for this interval
                for (int sample_idx = 0; sample_idx < n_samples; ++sample_idx) {
                    // Reset bin cursor for each sample. Bin queries use
                    // pos - k (the leftmost base of the (k+1)-mer
                    // context), matching the convention used by training.
                    // Initialize against the first query position the
                    // sampling loop will see.
                    size_t bin_cursor = 0;
                    int64_t first_query = interval_start;  // pos - k for pos = interval_start + k
                    if (!bins.empty()) {
                        while (bin_cursor + 1 < bins.size() &&
                               first_query >= bins[bin_cursor + 1].first) {
                            ++bin_cursor;
                        }
                    }

                    // Allocate output sequence
                    vector<char> synth_seq(interval_len);

                    size_t mask_cursor = 0;
                    while (mask_cursor < mask_copy_ivs.size() &&
                           mask_copy_ivs[mask_cursor].end <= interval_start) {
                        ++mask_cursor;
                    }

                    // Initialize first k bases, honoring mask_copy semantics.
                    int64_t init_len = min<int64_t>(k, interval_len);
                    for (int64_t i = 0; i < init_len; ++i) {
                        int64_t pos = interval_start + i;
                        // Check if this position should be copied from original
                        if (is_position_masked(pos, mask_copy_ivs, mask_cursor) &&
                            i < (int64_t)original_seq.size()) {
                            synth_seq[i] = original_seq[i];
                        } else {
                            synth_seq[i] = StratifiedMarkovModel::decode_base(
                                static_cast<int>(unif_rand() * NUM_BASES));
                        }
                    }

                    // Sample remaining bases using Markov chain
                    for (int64_t pos = interval_start + init_len; pos < interval_end; ++pos) {
                        int64_t rel_pos = pos - interval_start;

                        // Check if this position should be copied from original
                        if (is_position_masked(pos, mask_copy_ivs, mask_cursor)) {
                            if (rel_pos < (int64_t)original_seq.size()) {
                                // Copy from original
                                synth_seq[rel_pos] = original_seq[rel_pos];
                            } else {
                                // No original sequence available, sample uniformly
                                synth_seq[rel_pos] = StratifiedMarkovModel::decode_base(
                                    static_cast<int>(unif_rand() * NUM_BASES));
                            }
                            continue;
                        }

                        // Find bin for this position using a forward
                        // cursor. Use pos - k (context-leftmost) to match
                        // training; this keeps the per-bin transition
                        // tables aligned with the same iter window the
                        // training counts were attributed to.
                        int64_t bin_query_pos = pos - k;
                        int bin_idx = -1;
                        if (!bins.empty()) {
                            while (bin_cursor + 1 < bins.size() &&
                                   bin_query_pos >= bins[bin_cursor + 1].first) {
                                ++bin_cursor;
                            }
                            if (bin_query_pos >= bins[bin_cursor].first &&
                                bin_query_pos <
                                    bins[bin_cursor].first + iter_size) {
                                bin_idx = bins[bin_cursor].second;
                            }
                        }

                        // Get k-mer context from already-sampled bases
                        int context_idx = StratifiedMarkovModel::encode_kmer(
                            &synth_seq[rel_pos - k], k);

                        int next_base;
                        if (context_idx < 0 || bin_idx < 0 || bin_idx >= num_bins) {
                            // Invalid context or bin - sample uniformly
                            next_base = static_cast<int>(unif_rand() * NUM_BASES);
                        } else {
                            // Sample from CDF (flat layout: ctx * NUM_BASES + base)
                            float r = unif_rand();
                            int base_offset = context_idx * NUM_BASES;
                            next_base = NUM_BASES - 1;
                            for (int b = 0; b < NUM_BASES; ++b) {
                                if (r < cdf_data[bin_idx][base_offset + b]) {
                                    next_base = b;
                                    break;
                                }
                            }
                        }

                        synth_seq[rel_pos] =
                            StratifiedMarkovModel::decode_base(next_base);
                    }

                    // Write interval to output
                    if (output_format == 2) {
                        // Vector mode: collect sequence as string
                        collected_seqs.push_back(string(synth_seq.begin(), synth_seq.end()));
                    } else if (output_format == 1) {
                        // FASTA
                        string header = chrom_name;
                        if (!(interval_start == 0 && interval_end == chrom_size)) {
                            header = chrom_name + ":" +
                                to_string(interval_start) + "-" +
                                to_string(interval_end);
                        }
                        if (n_samples > 1) {
                            header += "_sample" + to_string(sample_idx + 1);
                        }
                        long long header_offset = static_cast<long long>(fasta_ofs.tellp());
                        FaiEntry entry;
                        write_fasta(fasta_ofs, header, synth_seq, 60,
                                    header_offset, &entry);
                        fai_entries.push_back(entry);
                    } else {
                        // .seq binary
                        write_seq(seq_bfile, synth_seq);
                    }

                    progress.report(interval_len);
                }  // end sample loop
            }

            check_interrupt();
        }

        progress.report_last();

        // Close output files
        if (output_format == 1) {
            fasta_ofs.close();
            // Write .fai alongside the FASTA so the output is a drop-in replacement
            // for samtools-indexed references. Matches the format produced by
            // `samtools faidx`: <name>\t<length>\t<offset>\t<linebases>\t<linewidth>.
            string fai_path = string(output_path) + ".fai";
            FILE* ffai = fopen(fai_path.c_str(), "w");
            if (!ffai) {
                verror("Failed to open .fai file for writing: %s", fai_path.c_str());
            }
            for (const auto& e : fai_entries) {
                fprintf(ffai, "%s\t%lld\t%lld\t%d\t%d\n",
                        e.name.c_str(), e.length, e.offset, e.linebases, e.linewidth);
            }
            fclose(ffai);
        } else if (output_format == 0) {
            seq_bfile.close();
        }

        // Save R's RNG state
        rng_guard.release();

        // Return vector of sequences if requested
        if (output_format == 2) {
            SEXP result;
            rprotect(result = Rf_allocVector(STRSXP, collected_seqs.size()));
            for (size_t i = 0; i < collected_seqs.size(); ++i) {
                SET_STRING_ELT(result, i, Rf_mkChar(collected_seqs[i].c_str()));
            }
            return result;
        }

        return R_NilValue;

    } catch (TGLException& e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc& e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

}  // extern "C"
