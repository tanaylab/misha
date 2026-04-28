/*
 * StratifiedMarkovModel.cpp
 *
 * Implementation of the stratified Markov model with configurable order k.
 */

#include "StratifiedMarkovModel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include "TGLException.h"

// Out-of-line definitions required for ODR-use in C++14 (taking address
// in ofs.write).  Harmless but redundant in C++17.
constexpr uint32_t StratifiedMarkovModel::FILE_MAGIC;
constexpr uint32_t StratifiedMarkovModel::FILE_VERSION;
constexpr uint32_t StratifiedMarkovModel::FILE_VERSION_1;
constexpr int StratifiedMarkovModel::MAX_K;
constexpr int StratifiedMarkovModel::MAX_KMERS;

StratifiedMarkovModel::StratifiedMarkovModel()
    : m_k(0), m_num_kmers(0), m_num_bins(0), m_total_kmers(0) {}

void StratifiedMarkovModel::init(int num_bins, const std::vector<double>& breaks, int k) {
    if (num_bins <= 0) {
        throw std::invalid_argument("num_bins must be positive");
    }
    if (breaks.size() != static_cast<size_t>(num_bins + 1)) {
        throw std::invalid_argument("breaks must have num_bins + 1 elements");
    }
    if (k < 1 || k > MAX_K) {
        throw std::invalid_argument("Markov order k must be in [1, 8]");
    }

    m_k = k;
    m_num_kmers = 1 << (2 * k);  // 4^k
    m_num_bins = num_bins;
    m_breaks = breaks;

    // Allocate count and CDF arrays as flat vectors
    m_counts.resize(num_bins);
    m_cdf.resize(num_bins);
    m_per_bin_kmers.resize(num_bins, 0);

    for (int b = 0; b < num_bins; ++b) {
        m_counts[b].assign(m_num_kmers * NUM_BASES, 0);
        m_cdf[b].assign(m_num_kmers * NUM_BASES, 0.0f);
    }

    // Default per-bin prior: uniform (will be overridden by setter before normalize)
    m_prior.assign(num_bins, {0.25, 0.25, 0.25, 0.25});

    m_total_kmers = 0;
}

void StratifiedMarkovModel::reset_counts() {
    m_total_kmers = 0;
    for (int b = 0; b < m_num_bins; ++b) {
        m_per_bin_kmers[b] = 0;
        std::fill(m_counts[b].begin(), m_counts[b].end(), 0);
        std::fill(m_cdf[b].begin(), m_cdf[b].end(), 0.0f);
    }
}

void StratifiedMarkovModel::increment_count(int bin_idx, int context_kmer_idx, int next_base_idx) {
    if (bin_idx < 0 || bin_idx >= m_num_bins) return;
    if (context_kmer_idx < 0 || context_kmer_idx >= m_num_kmers) return;
    if (next_base_idx < 0 || next_base_idx >= NUM_BASES) return;

    m_counts[bin_idx][context_kmer_idx * NUM_BASES + next_base_idx]++;
    m_per_bin_kmers[bin_idx]++;
    m_total_kmers++;
}

void StratifiedMarkovModel::apply_bin_mapping(const std::vector<int>& bin_map) {
    if (bin_map.empty()) return;
    if (bin_map.size() != static_cast<size_t>(m_num_bins)) {
        throw std::invalid_argument("bin_map size must match num_bins");
    }

    // Create temporary storage for merged counts
    std::vector<std::vector<uint64_t>> new_counts(m_num_bins);
    std::vector<uint64_t> new_per_bin_kmers(m_num_bins, 0);

    int flat_size = m_num_kmers * NUM_BASES;

    // Initialize new counts to zero
    for (int b = 0; b < m_num_bins; ++b) {
        new_counts[b].assign(flat_size, 0);
    }

    // Merge counts according to bin_map
    for (int src_bin = 0; src_bin < m_num_bins; ++src_bin) {
        int tgt_bin = bin_map[src_bin];
        if (tgt_bin < 0 || tgt_bin >= m_num_bins) {
            tgt_bin = src_bin; // Keep original if invalid target
        }

        for (int i = 0; i < flat_size; ++i) {
            new_counts[tgt_bin][i] += m_counts[src_bin][i];
        }
        new_per_bin_kmers[tgt_bin] += m_per_bin_kmers[src_bin];
    }

    // Replace old counts with merged counts
    m_counts = std::move(new_counts);
    m_per_bin_kmers = std::move(new_per_bin_kmers);
}

void StratifiedMarkovModel::set_prior_uniform() {
    for (auto& row : m_prior) {
        row = {0.25, 0.25, 0.25, 0.25};
    }
}

namespace {
void normalize_row_or_uniform(std::array<double, NUM_BASES>& row) {
    double s = 0.0;
    for (int a = 0; a < NUM_BASES; ++a) s += row[a];
    if (s <= 0.0) {
        row = {0.25, 0.25, 0.25, 0.25};
    } else {
        for (int a = 0; a < NUM_BASES; ++a) row[a] /= s;
    }
}
}  // namespace

void StratifiedMarkovModel::set_prior_global(
    const std::array<double, NUM_BASES>& pi) {
    std::array<double, NUM_BASES> row = pi;
    normalize_row_or_uniform(row);
    for (auto& r : m_prior) r = row;
}

void StratifiedMarkovModel::set_prior_explicit(
    const std::vector<std::array<double, NUM_BASES>>& pi_per_bin) {
    if (static_cast<int>(pi_per_bin.size()) != m_num_bins) {
        throw std::invalid_argument("prior matrix row count must match num_bins");
    }
    for (int b = 0; b < m_num_bins; ++b) {
        m_prior[b] = pi_per_bin[b];
        normalize_row_or_uniform(m_prior[b]);
    }
}

int StratifiedMarkovModel::set_prior_from_marginal() {
    int fallback_count = 0;
    for (int b = 0; b < m_num_bins; ++b) {
        std::array<double, NUM_BASES> sums = {0.0, 0.0, 0.0, 0.0};
        for (int ctx = 0; ctx < m_num_kmers; ++ctx) {
            int off = ctx * NUM_BASES;
            for (int a = 0; a < NUM_BASES; ++a) {
                sums[a] += static_cast<double>(m_counts[b][off + a]);
            }
        }
        double total = sums[0] + sums[1] + sums[2] + sums[3];
        if (total <= 0.0) {
            m_prior[b] = {0.25, 0.25, 0.25, 0.25};
            ++fallback_count;
        } else {
            for (int a = 0; a < NUM_BASES; ++a) m_prior[b][a] = sums[a] / total;
        }
    }
    return fallback_count;
}

bool StratifiedMarkovModel::set_prior_from_global_marginal() {
    std::array<double, NUM_BASES> sums = {0.0, 0.0, 0.0, 0.0};
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < m_num_kmers; ++ctx) {
            int off = ctx * NUM_BASES;
            for (int a = 0; a < NUM_BASES; ++a) {
                sums[a] += static_cast<double>(m_counts[b][off + a]);
            }
        }
    }
    double total = sums[0] + sums[1] + sums[2] + sums[3];
    if (total <= 0.0) {
        for (auto& r : m_prior) r = {0.25, 0.25, 0.25, 0.25};
        return false;
    }
    std::array<double, NUM_BASES> pi;
    for (int a = 0; a < NUM_BASES; ++a) pi[a] = sums[a] / total;
    for (auto& r : m_prior) r = pi;
    return true;
}

void StratifiedMarkovModel::normalize_and_build_cdf(double pseudocount) {
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < m_num_kmers; ++ctx) {
            int base_offset = ctx * NUM_BASES;

            // Calculate total count for this context (with pseudocounts)
            double total = 0.0;
            for (int base = 0; base < NUM_BASES; ++base) {
                total += static_cast<double>(m_counts[b][base_offset + base]) + pseudocount;
            }

            // Build CDF
            double cumsum = 0.0;
            for (int base = 0; base < NUM_BASES; ++base) {
                double prob = (static_cast<double>(m_counts[b][base_offset + base]) + pseudocount) / total;
                cumsum += prob;
                m_cdf[b][base_offset + base] = static_cast<float>(cumsum);
            }
            // Ensure last element is exactly 1.0 (avoid floating point issues)
            m_cdf[b][base_offset + NUM_BASES - 1] = 1.0f;
        }
    }
}

int StratifiedMarkovModel::get_bin(double track_value) const {
    if (std::isnan(track_value) || m_num_bins == 0) {
        return -1;
    }

    // Binary search for the bin
    // breaks: [b0, b1, b2, ..., bn] defines bins [b0,b1), [b1,b2), ..., [bn-1,bn]
    // We use left-closed, right-open intervals except for the last bin which is closed

    if (track_value < m_breaks[0] || track_value > m_breaks[m_num_bins]) {
        return -1; // Out of range
    }

    // Binary search
    auto it = std::upper_bound(m_breaks.begin(), m_breaks.end(), track_value);
    int bin = static_cast<int>(it - m_breaks.begin()) - 1;

    // Clamp to valid range
    if (bin < 0) bin = 0;
    if (bin >= m_num_bins) bin = m_num_bins - 1;

    return bin;
}

int StratifiedMarkovModel::sample_next_base(int bin_idx, int context_kmer_idx, float random_val) const {
    if (bin_idx < 0 || bin_idx >= m_num_bins) {
        // Fall back to uniform sampling
        return static_cast<int>(random_val * NUM_BASES);
    }
    if (context_kmer_idx < 0 || context_kmer_idx >= m_num_kmers) {
        return static_cast<int>(random_val * NUM_BASES);
    }

    // Use CDF for sampling
    int base_offset = context_kmer_idx * NUM_BASES;
    for (int base = 0; base < NUM_BASES; ++base) {
        if (random_val < m_cdf[bin_idx][base_offset + base]) {
            return base;
        }
    }
    return NUM_BASES - 1; // Fallback (shouldn't happen with proper CDF)
}

int StratifiedMarkovModel::encode_base(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

char StratifiedMarkovModel::decode_base(int idx) {
    static const char bases[] = "ACGT";
    if (idx < 0 || idx >= NUM_BASES) return 'N';
    return bases[idx];
}

int StratifiedMarkovModel::encode_kmer(const char* seq, int k) {
    int idx = 0;
    for (int i = 0; i < k; ++i) {
        int base = encode_base(seq[i]);
        if (base < 0) return -1;
        idx = (idx << 2) | base;
    }
    return idx;
}

void StratifiedMarkovModel::decode_kmer(int idx, char* out, int k) {
    for (int i = k - 1; i >= 0; --i) {
        out[i] = decode_base(idx & 3);
        idx >>= 2;
    }
}

int StratifiedMarkovModel::complement_base(int base_idx) {
    // A(0) <-> T(3), C(1) <-> G(2)
    return 3 - base_idx;
}

void StratifiedMarkovModel::revcomp_kmer(int context_kmer_idx, int next_base_idx, int k,
                                          int& revcomp_context_idx, int& revcomp_next_idx) {
    // Forward (k+1)-mer: B0 B1 ... B(k-1) B(k)
    //   context = B0 * 4^(k-1) + B1 * 4^(k-2) + ... + B(k-1)
    //   next = B(k)
    //
    // Reverse complement: comp(B(k)) comp(B(k-1)) ... comp(B1) comp(B0)
    //   revcomp_context = comp(B(k)) * 4^(k-1) + comp(B(k-1)) * 4^(k-2) + ... + comp(B1)
    //   revcomp_next = comp(B0)

    // Extract the k bases from context and build reverse complement
    // B0 is most significant (leftmost)
    // We need to extract each base, complement it, and rebuild in reverse order.

    // Extract bases from context_kmer_idx into an array
    // bases[0] = B0 (most significant), bases[k-1] = B(k-1) (least significant)
    int bases[MAX_K + 1];
    {
        int tmp = context_kmer_idx;
        for (int i = k - 1; i >= 0; --i) {
            bases[i] = tmp & 3;
            tmp >>= 2;
        }
    }
    bases[k] = next_base_idx;

    // Build reverse complement:
    //   revcomp_context = comp(B(k)) comp(B(k-1)) ... comp(B1)
    //   revcomp_next = comp(B0)
    revcomp_context_idx = 0;
    for (int i = k; i >= 1; --i) {
        revcomp_context_idx = (revcomp_context_idx << 2) | complement_base(bases[i]);
    }
    revcomp_next_idx = complement_base(bases[0]);
}

void StratifiedMarkovModel::save(const std::string& path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        TGLError("Failed to open file for writing: %s", path.c_str());
    }

    // Write magic and version
    ofs.write(reinterpret_cast<const char*>(&FILE_MAGIC), sizeof(FILE_MAGIC));
    ofs.write(reinterpret_cast<const char*>(&FILE_VERSION), sizeof(FILE_VERSION));

    // Write k and num_kmers (new in version 2)
    ofs.write(reinterpret_cast<const char*>(&m_k), sizeof(m_k));
    ofs.write(reinterpret_cast<const char*>(&m_num_kmers), sizeof(m_num_kmers));

    // Write number of bins
    ofs.write(reinterpret_cast<const char*>(&m_num_bins), sizeof(m_num_bins));

    // Write breaks
    for (const auto& b : m_breaks) {
        ofs.write(reinterpret_cast<const char*>(&b), sizeof(b));
    }

    // Write statistics
    ofs.write(reinterpret_cast<const char*>(&m_total_kmers), sizeof(m_total_kmers));
    for (int i = 0; i < m_num_bins; ++i) {
        ofs.write(reinterpret_cast<const char*>(&m_per_bin_kmers[i]), sizeof(m_per_bin_kmers[i]));
    }

    // Write counts (flat: num_kmers * NUM_BASES uint64_t values per bin)
    for (int b = 0; b < m_num_bins; ++b) {
        ofs.write(reinterpret_cast<const char*>(m_counts[b].data()),
                  m_num_kmers * NUM_BASES * sizeof(uint64_t));
    }

    // Write CDFs (flat: num_kmers * NUM_BASES float values per bin)
    for (int b = 0; b < m_num_bins; ++b) {
        ofs.write(reinterpret_cast<const char*>(m_cdf[b].data()),
                  m_num_kmers * NUM_BASES * sizeof(float));
    }

    if (!ofs) {
        TGLError("Error writing to file: %s", path.c_str());
    }
}

StratifiedMarkovModel StratifiedMarkovModel::load(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        TGLError("Failed to open file for reading: %s", path.c_str());
    }

    StratifiedMarkovModel model;

    // Read and verify magic
    uint32_t magic;
    ifs.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (magic != FILE_MAGIC) {
        TGLError("Invalid file format: %s", path.c_str());
    }

    // Read and verify version
    uint32_t version;
    ifs.read(reinterpret_cast<char*>(&version), sizeof(version));

    if (version == FILE_VERSION_1) {
        // Legacy version 1: hardcoded k=5, num_kmers=1024
        model.m_k = 5;
        model.m_num_kmers = 1024;
    } else if (version == FILE_VERSION) {
        // Version 2: read k and num_kmers from header
        ifs.read(reinterpret_cast<char*>(&model.m_k), sizeof(model.m_k));
        ifs.read(reinterpret_cast<char*>(&model.m_num_kmers), sizeof(model.m_num_kmers));
        if (model.m_k < 1 || model.m_k > MAX_K) {
            TGLError("Invalid Markov order k=%d in file: %s", model.m_k, path.c_str());
        }
        int expected_num_kmers = 1 << (2 * model.m_k);
        if (model.m_num_kmers != expected_num_kmers) {
            TGLError("Inconsistent num_kmers (%d, expected %d for k=%d) in file: %s",
                     model.m_num_kmers, expected_num_kmers, model.m_k, path.c_str());
        }
    } else {
        TGLError("Unsupported file version %u in: %s", version, path.c_str());
    }

    // Read number of bins
    ifs.read(reinterpret_cast<char*>(&model.m_num_bins), sizeof(model.m_num_bins));

    // Read breaks
    model.m_breaks.resize(model.m_num_bins + 1);
    for (int i = 0; i <= model.m_num_bins; ++i) {
        ifs.read(reinterpret_cast<char*>(&model.m_breaks[i]), sizeof(double));
    }

    // Read statistics
    ifs.read(reinterpret_cast<char*>(&model.m_total_kmers), sizeof(model.m_total_kmers));
    model.m_per_bin_kmers.resize(model.m_num_bins);
    for (int i = 0; i < model.m_num_bins; ++i) {
        ifs.read(reinterpret_cast<char*>(&model.m_per_bin_kmers[i]), sizeof(uint64_t));
    }

    int num_kmers = model.m_num_kmers;
    int flat_size = num_kmers * NUM_BASES;

    // Read counts
    // Both version 1 and 2 use the same contiguous layout:
    // [ctx0_base0, ctx0_base1, ..., ctxN_base3]
    model.m_counts.resize(model.m_num_bins);
    for (int b = 0; b < model.m_num_bins; ++b) {
        model.m_counts[b].resize(flat_size);
        ifs.read(reinterpret_cast<char*>(model.m_counts[b].data()),
                 flat_size * sizeof(uint64_t));
    }

    // Read CDFs
    model.m_cdf.resize(model.m_num_bins);
    for (int b = 0; b < model.m_num_bins; ++b) {
        model.m_cdf[b].resize(flat_size);
        ifs.read(reinterpret_cast<char*>(model.m_cdf[b].data()),
                 flat_size * sizeof(float));
    }

    if (!ifs) {
        TGLError("Error reading from file: %s", path.c_str());
    }

    return model;
}
