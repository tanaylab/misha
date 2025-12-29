/*
 * StratifiedMarkovModel.h
 *
 * Stratified Markov-5 model for synthetic genome generation.
 * Stores transition probabilities for 5-mer -> base transitions,
 * stratified by bins of a track expression (e.g., GC content).
 */

#ifndef STRATIFIED_MARKOV_MODEL_H_
#define STRATIFIED_MARKOV_MODEL_H_

#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Number of possible 5-mers (4^5 = 1024)
constexpr int NUM_5MERS = 1024;

// Number of bases (A, C, G, T)
constexpr int NUM_BASES = 4;

/**
 * StratifiedMarkovModel: A 5th-order Markov model stratified by track values.
 *
 * For each bin (based on a track expression like GC content), stores:
 * - Counts of 6-mers (5-mer context + next base)
 * - Probability distributions P(next_base | 5-mer_context)
 * - Cumulative distribution functions for sampling
 */
class StratifiedMarkovModel {
public:
    StratifiedMarkovModel();

    /**
     * Initialize the model with the specified number of bins.
     * @param num_bins Number of stratification bins
     * @param breaks Bin boundaries (size = num_bins + 1)
     */
    void init(int num_bins, const std::vector<double>& breaks);

    /**
     * Reset all counts to zero.
     */
    void reset_counts();

    /**
     * Increment the count for a 6-mer in the specified bin.
     * @param bin_idx Bin index (0 to num_bins-1)
     * @param context_5mer_idx Encoded 5-mer context (0 to 1023)
     * @param next_base_idx Next base index (0=A, 1=C, 2=G, 3=T)
     */
    void increment_count(int bin_idx, int context_5mer_idx, int next_base_idx);

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
     * Sample the next base given a 5-mer context and bin.
     * @param bin_idx Bin index
     * @param context_5mer_idx Encoded 5-mer context
     * @param random_val Random value in [0, 1)
     * @return Base index (0=A, 1=C, 2=G, 3=T)
     */
    int sample_next_base(int bin_idx, int context_5mer_idx, float random_val) const;

    /**
     * Encode a 5-mer sequence to an index (0-1023).
     * @param seq Pointer to 5 characters (A/C/G/T)
     * @return Index (0-1023), or -1 if any character is not A/C/G/T
     */
    static int encode_5mer(const char* seq);

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
     * Decode a 5-mer index to a string.
     * @param idx 5-mer index (0-1023)
     * @param out Output buffer (must have space for 5 characters)
     */
    static void decode_5mer(int idx, char* out);

    /**
     * Get the complement of a base index.
     * @param base_idx Base index (0=A, 1=C, 2=G, 3=T)
     * @return Complement base index (A↔T, C↔G)
     */
    static int complement_base(int base_idx);

    /**
     * Compute the reverse complement of a 6-mer given as context + next base.
     * Given forward 6-mer B0B1B2B3B4B5 (context=B0-B4, next=B5),
     * computes reverse complement: comp(B5)comp(B4)comp(B3)comp(B2)comp(B1)comp(B0)
     * where the new context is comp(B5)-comp(B1) and new next is comp(B0).
     *
     * @param context_5mer_idx Encoded 5-mer context (0 to 1023)
     * @param next_base_idx Next base index (0=A, 1=C, 2=G, 3=T)
     * @param revcomp_context_idx [out] Reverse complement 5-mer context index
     * @param revcomp_next_idx [out] Reverse complement next base index
     */
    static void revcomp_6mer(int context_5mer_idx, int next_base_idx,
                             int& revcomp_context_idx, int& revcomp_next_idx);

    // Accessors
    int get_num_bins() const { return m_num_bins; }
    const std::vector<double>& get_breaks() const { return m_breaks; }
    uint64_t get_total_kmers() const { return m_total_kmers; }
    uint64_t get_bin_kmers(int bin_idx) const { return m_per_bin_kmers[bin_idx]; }

    // Data accessors for R interface
    uint64_t get_count(int bin_idx, int context_idx, int base_idx) const {
        return m_counts[bin_idx][context_idx][base_idx];
    }
    float get_cdf(int bin_idx, int context_idx, int base_idx) const {
        return m_cdf[bin_idx][context_idx][base_idx];
    }

    // Direct access for bulk data transfer
    const std::vector<std::array<std::array<uint64_t, NUM_BASES>, NUM_5MERS>>& get_counts() const {
        return m_counts;
    }
    const std::vector<std::array<std::array<float, NUM_BASES>, NUM_5MERS>>& get_cdf() const {
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
    int m_num_bins;
    std::vector<double> m_breaks;

    // counts[bin][5mer_idx][base_idx]
    std::vector<std::array<std::array<uint64_t, NUM_BASES>, NUM_5MERS>> m_counts;

    // Cumulative distribution functions for sampling: cdf[bin][5mer_idx][base_idx]
    // cdf[bin][ctx][0] = P(A), cdf[bin][ctx][1] = P(A)+P(C), etc.
    std::vector<std::array<std::array<float, NUM_BASES>, NUM_5MERS>> m_cdf;

    // Statistics
    uint64_t m_total_kmers;
    std::vector<uint64_t> m_per_bin_kmers;

    // Magic number for file format
    static constexpr uint32_t FILE_MAGIC = 0x4D4B5653; // "SVKM" in little-endian
    static constexpr uint32_t FILE_VERSION = 1;
};

#endif // STRATIFIED_MARKOV_MODEL_H_
