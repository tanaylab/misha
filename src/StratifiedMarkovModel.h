/*
 * StratifiedMarkovModel.h
 *
 * Stratified Markov model for synthetic genome generation.
 * Stores transition probabilities for k-mer -> base transitions,
 * stratified by bins of a track expression (e.g., GC content).
 * Supports configurable Markov order k (1..8).
 */

#ifndef STRATIFIED_MARKOV_MODEL_H_
#define STRATIFIED_MARKOV_MODEL_H_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Kept for backward compatibility with external callers.
// Do NOT rely on NUM_5MERS inside StratifiedMarkovModel internals.
constexpr int NUM_5MERS = 1024;

// Number of bases (A, C, G, T)
constexpr int NUM_BASES = 4;

/**
 * StratifiedMarkovModel: A k-th order Markov model stratified by track values.
 *
 * For each bin (based on a track expression like GC content), stores:
 * - Counts of (k+1)-mers (k-mer context + next base)
 * - Probability distributions P(next_base | k-mer_context)
 * - Cumulative distribution functions for sampling
 *
 * The Markov order k is configurable at init time (1..8, default 5).
 */
class StratifiedMarkovModel {
public:
    // Maximum supported Markov order
    static constexpr int MAX_K = 10;
    // Maximum number of k-mers (4^MAX_K = 1048576)
    static constexpr int MAX_KMERS = 1048576;

    StratifiedMarkovModel();

    /**
     * Initialize the model with the specified number of bins and Markov order.
     * @param num_bins Number of stratification bins
     * @param breaks Bin boundaries (size = num_bins + 1)
     * @param k Markov order (1..8, default 5 for backward compatibility)
     */
    void init(int num_bins, const std::vector<double>& breaks, int k = 5);

    /**
     * Reset all counts to zero.
     */
    void reset_counts();

    /**
     * Increment the count for a (k+1)-mer in the specified bin.
     * @param bin_idx Bin index (0 to num_bins-1)
     * @param context_kmer_idx Encoded k-mer context (0 to m_num_kmers-1)
     * @param next_base_idx Next base index (0=A, 1=C, 2=G, 3=T)
     */
    void increment_count(int bin_idx, int context_kmer_idx, int next_base_idx);

    /**
     * Apply bin mapping: merge sparse bins into target bins.
     * @param bin_map Maps source bin indices to target bin indices.
     *                bin_map[i] = j means bin i's counts go to bin j.
     *                If bin_map[i] < 0, bin i keeps its own data.
     */
    void apply_bin_mapping(const std::vector<int>& bin_map);

    /**
     * Normalize counts to probabilities and build CDFs for sampling.
     * Should be called after all counts are accumulated and bin mapping is applied.
     * @param pseudocount Value to add to all counts to avoid zero probabilities (default 1)
     */
    void normalize_and_build_cdf(double pseudocount = 1.0);

    /**
     * Get the bin index for a track value.
     * @param track_value The track value to look up
     * @return Bin index (0 to num_bins-1), or -1 if value is NA or out of range
     */
    int get_bin(double track_value) const;

    /**
     * Sample the next base given a k-mer context and bin.
     * @param bin_idx Bin index
     * @param context_kmer_idx Encoded k-mer context
     * @param random_val Random value in [0, 1)
     * @return Base index (0=A, 1=C, 2=G, 3=T)
     */
    int sample_next_base(int bin_idx, int context_kmer_idx, float random_val) const;

    /**
     * Encode a k-mer sequence to an index.
     * @param seq Pointer to k characters (A/C/G/T)
     * @param k Number of characters to encode
     * @return Index (0 to 4^k - 1), or -1 if any character is not A/C/G/T
     */
    static int encode_kmer(const char* seq, int k);

    /**
     * Encode a 5-mer sequence to an index (0-1023).
     * Backward-compatible wrapper around encode_kmer with k=5.
     * @param seq Pointer to 5 characters (A/C/G/T)
     * @return Index (0-1023), or -1 if any character is not A/C/G/T
     */
    static inline int encode_5mer(const char* seq) {
        return encode_kmer(seq, 5);
    }

    /**
     * Encode a single base to an index.
     * @param base Character (A/C/G/T, case-insensitive)
     * @return Index (0=A, 1=C, 2=G, 3=T), or -1 if not a valid base
     */
    static int encode_base(char base);

    /**
     * Decode a base index to a character.
     * @param idx Base index (0-3)
     * @return Character ('A', 'C', 'G', or 'T')
     */
    static char decode_base(int idx);

    /**
     * Decode a k-mer index to a string.
     * @param idx k-mer index (0 to 4^k - 1)
     * @param out Output buffer (must have space for k characters)
     * @param k Number of characters to decode
     */
    static void decode_kmer(int idx, char* out, int k);

    /**
     * Decode a 5-mer index to a string.
     * Backward-compatible wrapper around decode_kmer with k=5.
     * @param idx 5-mer index (0-1023)
     * @param out Output buffer (must have space for 5 characters)
     */
    static inline void decode_5mer(int idx, char* out) {
        decode_kmer(idx, out, 5);
    }

    /**
     * Get the complement of a base index.
     * @param base_idx Base index (0=A, 1=C, 2=G, 3=T)
     * @return Complement base index (A<->T, C<->G)
     */
    static int complement_base(int base_idx);

    /**
     * Compute the reverse complement of a (k+1)-mer given as context + next base.
     * Given forward (k+1)-mer B0 B1 ... B(k-1) B(k) (context=B0..B(k-1), next=B(k)),
     * computes reverse complement:
     *   comp(B(k)) comp(B(k-1)) ... comp(B1) comp(B0)
     * where the new context is comp(B(k))..comp(B1) and new next is comp(B0).
     *
     * @param context_kmer_idx Encoded k-mer context (0 to m_num_kmers-1)
     * @param next_base_idx Next base index (0=A, 1=C, 2=G, 3=T)
     * @param k Markov order (context length)
     * @param revcomp_context_idx [out] Reverse complement k-mer context index
     * @param revcomp_next_idx [out] Reverse complement next base index
     */
    static void revcomp_kmer(int context_kmer_idx, int next_base_idx, int k,
                             int& revcomp_context_idx, int& revcomp_next_idx);

    /**
     * Compute reverse complement of a 6-mer (k=5).
     * Backward-compatible wrapper around revcomp_kmer with k=5.
     */
    static inline void revcomp_6mer(int context_5mer_idx, int next_base_idx,
                                    int& revcomp_context_idx, int& revcomp_next_idx) {
        revcomp_kmer(context_5mer_idx, next_base_idx, 5,
                     revcomp_context_idx, revcomp_next_idx);
    }

    // Accessors
    int get_k() const { return m_k; }
    int get_num_kmers() const { return m_num_kmers; }
    int get_num_bins() const { return m_num_bins; }
    const std::vector<double>& get_breaks() const { return m_breaks; }
    uint64_t get_total_kmers() const { return m_total_kmers; }
    uint64_t get_bin_kmers(int bin_idx) const { return m_per_bin_kmers[bin_idx]; }

    // Data accessors for R interface
    uint64_t get_count(int bin_idx, int context_idx, int base_idx) const {
        return m_counts[bin_idx][context_idx * NUM_BASES + base_idx];
    }
    float get_cdf(int bin_idx, int context_idx, int base_idx) const {
        return m_cdf[bin_idx][context_idx * NUM_BASES + base_idx];
    }

    // Direct access for bulk data transfer
    const std::vector<std::vector<uint64_t>>& get_counts() const {
        return m_counts;
    }
    const std::vector<std::vector<float>>& get_cdf() const {
        return m_cdf;
    }

    /**
     * Serialize the model to a binary file.
     * @param path Output file path
     */
    void save(const std::string& path) const;

    /**
     * Load a model from a binary file.
     * @param path Input file path
     * @return Loaded model
     */
    static StratifiedMarkovModel load(const std::string& path);

private:
    int m_k;          // Markov order (1..8)
    int m_num_kmers;  // = 4^k, computed at init time
    int m_num_bins;
    std::vector<double> m_breaks;

    // Flat vectors: counts[bin][kmer_idx * NUM_BASES + base_idx]
    std::vector<std::vector<uint64_t>> m_counts;

    // Cumulative distribution functions for sampling: cdf[bin][kmer_idx * NUM_BASES + base_idx]
    // For a given context ctx: cdf[bin][ctx*4+0] = P(A), cdf[bin][ctx*4+1] = P(A)+P(C), etc.
    std::vector<std::vector<float>> m_cdf;

    // Statistics
    uint64_t m_total_kmers;
    std::vector<uint64_t> m_per_bin_kmers;

    // Magic number for file format
    static constexpr uint32_t FILE_MAGIC = 0x4D4B5653; // "SVKM" in little-endian
    // Version 2: stores k and num_kmers in header
    static constexpr uint32_t FILE_VERSION = 2;
    // Version 1 for backward-compatible loading of old files
    static constexpr uint32_t FILE_VERSION_1 = 1;
};

#endif // STRATIFIED_MARKOV_MODEL_H_
