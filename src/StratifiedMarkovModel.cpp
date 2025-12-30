/*
 * StratifiedMarkovModel.cpp
 *
 * Implementation of the stratified Markov-5 model.
 */

#include "StratifiedMarkovModel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include "TGLException.h"

StratifiedMarkovModel::StratifiedMarkovModel()
    : m_num_bins(0), m_total_kmers(0) {}

void StratifiedMarkovModel::init(int num_bins, const std::vector<double>& breaks) {
    if (num_bins <= 0) {
        throw std::invalid_argument("num_bins must be positive");
    }
    if (breaks.size() != static_cast<size_t>(num_bins + 1)) {
        throw std::invalid_argument("breaks must have num_bins + 1 elements");
    }

    m_num_bins = num_bins;
    m_breaks = breaks;

    // Allocate count and CDF arrays
    m_counts.resize(num_bins);
    m_cdf.resize(num_bins);
    m_per_bin_kmers.resize(num_bins, 0);

    reset_counts();
}

void StratifiedMarkovModel::reset_counts() {
    m_total_kmers = 0;
    for (int b = 0; b < m_num_bins; ++b) {
        m_per_bin_kmers[b] = 0;
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            for (int base = 0; base < NUM_BASES; ++base) {
                m_counts[b][ctx][base] = 0;
                m_cdf[b][ctx][base] = 0.0f;
            }
        }
    }
}

void StratifiedMarkovModel::increment_count(int bin_idx, int context_5mer_idx, int next_base_idx) {
    if (bin_idx < 0 || bin_idx >= m_num_bins) return;
    if (context_5mer_idx < 0 || context_5mer_idx >= NUM_5MERS) return;
    if (next_base_idx < 0 || next_base_idx >= NUM_BASES) return;

    m_counts[bin_idx][context_5mer_idx][next_base_idx]++;
    m_per_bin_kmers[bin_idx]++;
    m_total_kmers++;
}

void StratifiedMarkovModel::apply_bin_mapping(const std::vector<int>& bin_map) {
    if (bin_map.empty()) return;
    if (bin_map.size() != static_cast<size_t>(m_num_bins)) {
        throw std::invalid_argument("bin_map size must match num_bins");
    }

    // Create temporary storage for merged counts
    std::vector<std::array<std::array<uint64_t, NUM_BASES>, NUM_5MERS>> new_counts(m_num_bins);
    std::vector<uint64_t> new_per_bin_kmers(m_num_bins, 0);

    // Initialize new counts to zero
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            for (int base = 0; base < NUM_BASES; ++base) {
                new_counts[b][ctx][base] = 0;
            }
        }
    }

    // Merge counts according to bin_map
    for (int src_bin = 0; src_bin < m_num_bins; ++src_bin) {
        int tgt_bin = bin_map[src_bin];
        if (tgt_bin < 0 || tgt_bin >= m_num_bins) {
            tgt_bin = src_bin; // Keep original if invalid target
        }

        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            for (int base = 0; base < NUM_BASES; ++base) {
                new_counts[tgt_bin][ctx][base] += m_counts[src_bin][ctx][base];
            }
        }
        new_per_bin_kmers[tgt_bin] += m_per_bin_kmers[src_bin];
    }

    // Replace old counts with merged counts
    m_counts = std::move(new_counts);
    m_per_bin_kmers = std::move(new_per_bin_kmers);
}

void StratifiedMarkovModel::normalize_and_build_cdf(double pseudocount) {
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            // Calculate total count for this context (with pseudocounts)
            double total = 0.0;
            for (int base = 0; base < NUM_BASES; ++base) {
                total += static_cast<double>(m_counts[b][ctx][base]) + pseudocount;
            }

            // Build CDF
            double cumsum = 0.0;
            for (int base = 0; base < NUM_BASES; ++base) {
                double prob = (static_cast<double>(m_counts[b][ctx][base]) + pseudocount) / total;
                cumsum += prob;
                m_cdf[b][ctx][base] = static_cast<float>(cumsum);
            }
            // Ensure last element is exactly 1.0 (avoid floating point issues)
            m_cdf[b][ctx][NUM_BASES - 1] = 1.0f;
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

int StratifiedMarkovModel::sample_next_base(int bin_idx, int context_5mer_idx, float random_val) const {
    if (bin_idx < 0 || bin_idx >= m_num_bins) {
        // Fall back to uniform sampling
        return static_cast<int>(random_val * NUM_BASES);
    }
    if (context_5mer_idx < 0 || context_5mer_idx >= NUM_5MERS) {
        return static_cast<int>(random_val * NUM_BASES);
    }

    // Use CDF for sampling
    const auto& cdf = m_cdf[bin_idx][context_5mer_idx];
    for (int base = 0; base < NUM_BASES; ++base) {
        if (random_val < cdf[base]) {
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

int StratifiedMarkovModel::encode_5mer(const char* seq) {
    int idx = 0;
    for (int i = 0; i < 5; ++i) {
        int base = encode_base(seq[i]);
        if (base < 0) return -1;
        idx = (idx << 2) | base;
    }
    return idx;
}

void StratifiedMarkovModel::decode_5mer(int idx, char* out) {
    for (int i = 4; i >= 0; --i) {
        out[i] = decode_base(idx & 3);
        idx >>= 2;
    }
}

int StratifiedMarkovModel::complement_base(int base_idx) {
    // A(0) <-> T(3), C(1) <-> G(2)
    return 3 - base_idx;
}

void StratifiedMarkovModel::revcomp_6mer(int context_5mer_idx, int next_base_idx,
                                          int& revcomp_context_idx, int& revcomp_next_idx) {
    // Forward 6-mer: B0 B1 B2 B3 B4 B5
    //   context = B0*256 + B1*64 + B2*16 + B3*4 + B4
    //   next = B5
    //
    // Reverse complement: comp(B5) comp(B4) comp(B3) comp(B2) comp(B1) comp(B0)
    //   revcomp_context = comp(B5)*256 + comp(B4)*64 + comp(B3)*16 + comp(B2)*4 + comp(B1)
    //   revcomp_next = comp(B0)

    // Extract the 5 bases from context (B0 is most significant)
    int b4 = context_5mer_idx & 3;
    int b3 = (context_5mer_idx >> 2) & 3;
    int b2 = (context_5mer_idx >> 4) & 3;
    int b1 = (context_5mer_idx >> 6) & 3;
    int b0 = (context_5mer_idx >> 8) & 3;
    int b5 = next_base_idx;

    // Compute complements
    int c0 = complement_base(b0);
    int c1 = complement_base(b1);
    int c2 = complement_base(b2);
    int c3 = complement_base(b3);
    int c4 = complement_base(b4);
    int c5 = complement_base(b5);

    // Build reverse complement context: c5 c4 c3 c2 c1
    revcomp_context_idx = (c5 << 8) | (c4 << 6) | (c3 << 4) | (c2 << 2) | c1;
    revcomp_next_idx = c0;
}

void StratifiedMarkovModel::save(const std::string& path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        TGLError("Failed to open file for writing: %s", path.c_str());
    }

    // Write magic and version
    ofs.write(reinterpret_cast<const char*>(&FILE_MAGIC), sizeof(FILE_MAGIC));
    ofs.write(reinterpret_cast<const char*>(&FILE_VERSION), sizeof(FILE_VERSION));

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

    // Write counts
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            ofs.write(reinterpret_cast<const char*>(m_counts[b][ctx].data()),
                      NUM_BASES * sizeof(uint64_t));
        }
    }

    // Write CDFs
    for (int b = 0; b < m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            ofs.write(reinterpret_cast<const char*>(m_cdf[b][ctx].data()),
                      NUM_BASES * sizeof(float));
        }
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
    if (version != FILE_VERSION) {
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

    // Read counts
    model.m_counts.resize(model.m_num_bins);
    for (int b = 0; b < model.m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            ifs.read(reinterpret_cast<char*>(model.m_counts[b][ctx].data()),
                     NUM_BASES * sizeof(uint64_t));
        }
    }

    // Read CDFs
    model.m_cdf.resize(model.m_num_bins);
    for (int b = 0; b < model.m_num_bins; ++b) {
        for (int ctx = 0; ctx < NUM_5MERS; ++ctx) {
            ifs.read(reinterpret_cast<char*>(model.m_cdf[b][ctx].data()),
                     NUM_BASES * sizeof(float));
        }
    }

    if (!ifs) {
        TGLError("Error reading from file: %s", path.c_str());
    }

    return model;
}
