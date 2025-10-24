#include "PWMScorer.h"
#include <algorithm>
#include <limits>
#include <cstring> // For strcmp
#include <cmath>   // For log

// Forward declaration for log_sum_log (defined in util.h)
extern inline void log_sum_log(float& a, float b);

// Compute log-likelihood at position i with spatial weighting
// Combines forward and reverse complement strands using logsumexp if bidirectional
// This is used for TOTAL_LIKELIHOOD mode
static inline float pos_value_with_spat(const DnaPSSM& pssm,
                                        const std::string& target,
                                        size_t i,
                                        char strand_mode,      // 1=plus, -1=minus, 0=both
                                        float spat_log)        // log(spatial_factor)
{
    float has = -std::numeric_limits<float>::infinity();

    // Evaluate forward strand
    float f = 0.f;
    auto it = target.begin() + i;
    pssm.calc_like(it, f);
    has = f;

    // Evaluate reverse strand if bidirectional
    if (pssm.is_bidirect()) {
        float rc = 0.f;
        auto it2 = target.begin() + i;
        pssm.calc_like_rc(it2, rc);
        log_sum_log(has, rc);
    }

    return has + spat_log;
}

// Compute best log-likelihood and its direction at position i with spatial weighting
// For MAX_LIKELIHOOD, MAX_LIKELIHOOD_POS, and MOTIF_COUNT modes
// Returns the best value and sets best_dir to 1 (forward) or -1 (reverse)
// NOTE: For bidirectional PSSMs, ALWAYS checks both strands (ignoring strand_mode)
// to match the behavior of DnaPSSM::max_like_match()
// NOTE: When strand_mode == -1, the target is already reverse-complemented, so we need
// to invert the direction: calc_like() scores the reverse strand, calc_like_rc() scores forward
static inline float pos_value_with_dir(const DnaPSSM& pssm,
                                       const std::string& target,
                                       size_t i,
                                       char strand_mode,
                                       float spat_log,
                                       int& best_dir)
{
    float best_val = -std::numeric_limits<float>::infinity();
    best_dir = 1;

    // For bidirectional PSSMs, always evaluate both strands (matching max_like_match behavior)
    // For non-bidirectional PSSMs, respect strand_mode
    bool check_forward = pssm.is_bidirect() || (strand_mode != -1);
    bool check_reverse = pssm.is_bidirect() || (strand_mode != 1);

    // Evaluate forward strand
    if (check_forward) {
        float f = 0.f;
        auto it = target.begin() + i;
        pssm.calc_like(it, f);
        float val_f = f + spat_log;
        if (val_f > best_val) {
            best_val = val_f;
            best_dir = 1;
        }
    }

    // Evaluate reverse strand
    if (check_reverse) {
        float rc = 0.f;
        auto it2 = target.begin() + i;
        pssm.calc_like_rc(it2, rc);
        float val_rc = rc + spat_log;
        if (val_rc > best_val) {
            best_val = val_rc;
            best_dir = -1;
        }
    }

    return best_val;
}

PWMScorer::PWMScorer(const DnaPSSM& pssm, const std::string& genome_root, bool extend,
                     ScoringMode mode, char strand,
                     const std::vector<float>& spat_factor, int spat_bin_size, float score_thresh)
    : GenomeSeqScorer(genome_root, extend, strand), m_pssm(pssm), m_mode(mode), m_score_thresh(score_thresh)
{
    if (!spat_factor.empty()) {
        m_use_spat = true;
        m_spat_bin_size = std::max(1, spat_bin_size);

        // Precompute log of spatial factors
        m_spat_log_factors.resize(spat_factor.size());
        for (size_t i = 0; i < spat_factor.size(); ++i) {
            m_spat_log_factors[i] = std::log(std::max(1e-30f, spat_factor[i]));
        }
    }
}

PWMScorer::PWMScorer(const DnaPSSM& pssm, GenomeSeqFetch* shared_seqfetch, bool extend,
                     ScoringMode mode, char strand,
                     const std::vector<float>& spat_factor, int spat_bin_size, float score_thresh)
    : GenomeSeqScorer(shared_seqfetch, extend, strand), m_pssm(pssm), m_mode(mode), m_score_thresh(score_thresh)
{
    if (!spat_factor.empty()) {
        m_use_spat = true;
        m_spat_bin_size = std::max(1, spat_bin_size);

        // Precompute log of spatial factors
        m_spat_log_factors.resize(spat_factor.size());
        for (size_t i = 0; i < spat_factor.size(); ++i) {
            m_spat_log_factors[i] = std::log(std::max(1e-30f, spat_factor[i]));
        }
    }
}

void PWMScorer::invalidate_cache()
{
    m_slide.valid = false;
    m_slide.stride = 0;
}

// Get spatial log factor for a given position index
inline float PWMScorer::get_spatial_log_factor(size_t pos_index) const
{
    if (!m_use_spat || m_spat_log_factors.empty()) {
        return 0.0f;
    }
    
    const int bin_sz = std::max(1, m_spat_bin_size);
    size_t bin = pos_index / size_t(bin_sz);
    if (bin >= m_spat_log_factors.size()) {
        bin = m_spat_log_factors.size() - 1;
    }
    return m_spat_log_factors[bin];
}

// Compute position result with strand and direction adjustments
float PWMScorer::compute_position_result(size_t index, size_t target_length, 
                                         size_t motif_length, int direction) const
{
    float pos_result = float(index) + 1.0f; // 1-based
    
    if (m_strand == -1) {
        pos_result = target_length - pos_result - motif_length + 1;
    }
    
    if (m_pssm.is_bidirect()) {
        pos_result = pos_result * direction;
    }
    
    return pos_result;
}

// Count motif hits without spatial weighting
float PWMScorer::count_motif_hits_no_spatial(const std::string& target, size_t motif_length)
{
    if (motif_length <= 0 || target.empty() || target.length() < motif_length) {
        return 0.0f;
    }

    int count = 0;
    size_t max_pos = target.length() - motif_length;

    // Scan forward strand
    if (m_strand != -1) {
        for (size_t i = 0; i <= max_pos; ++i) {
            float logp = 0;
            std::string::const_iterator it = target.begin() + i;
            m_pssm.calc_like(it, logp);
            if (logp >= m_score_thresh) {
                count++;
            }
        }
    }

    // Scan reverse strand if bidirectional
    if (m_pssm.is_bidirect() && m_strand != 1) {
        for (size_t i = 0; i <= max_pos; ++i) {
            float logp_rc = 0;
            std::string::const_iterator it = target.begin() + i;
            m_pssm.calc_like_rc(it, logp_rc);
            if (logp_rc >= m_score_thresh) {
                count++;
            }
        }
    }

    return static_cast<float>(count);
}

// Count motif hits with spatial weighting
float PWMScorer::count_motif_hits_with_spatial(const std::string& target, size_t motif_length)
{
    if (target.length() < motif_length) {
        return 0.0f;
    }

    int count = 0;
    size_t max_pos = target.length() - motif_length;

    // Scan forward strand
    if (m_strand != -1) {
        for (size_t i = 0; i <= max_pos; ++i) {
            int spat_bin = int(i / m_spat_bin_size);
            if (spat_bin >= (int)m_spat_log_factors.size()) {
                spat_bin = m_spat_log_factors.size() - 1;
            }
            float spat_log = m_spat_log_factors[spat_bin];

            float logp = 0;
            std::string::const_iterator it = target.begin() + i;
            m_pssm.calc_like(it, logp);
            float val = logp + spat_log;
            if (val >= m_score_thresh) {
                count++;
            }
        }
    }

    // Scan reverse strand if bidirectional
    if (m_pssm.is_bidirect() && m_strand != 1) {
        for (size_t i = 0; i <= max_pos; ++i) {
            int spat_bin = int(i / m_spat_bin_size);
            if (spat_bin >= (int)m_spat_log_factors.size()) {
                spat_bin = m_spat_log_factors.size() - 1;
            }
            float spat_log = m_spat_log_factors[spat_bin];

            float logp_rc = 0;
            std::string::const_iterator it = target.begin() + i;
            m_pssm.calc_like_rc(it, logp_rc);
            float val_rc = logp_rc + spat_log;
            if (val_rc >= m_score_thresh) {
                count++;
            }
        }
    }

    return static_cast<float>(count);
}

// Get max likelihood position with spatial weighting
float PWMScorer::get_max_likelihood_pos_with_spatial(const std::string& target, size_t motif_length)
{
    if (target.length() < motif_length) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    float best_val = -std::numeric_limits<float>::infinity();
    size_t best_index = 0;
    int best_dir = 1;

    size_t max_i_idx = std::min<size_t>(m_pssm.get_max_range(), target.size() - motif_length);
    size_t min_i_idx = std::min<size_t>(std::max(0, m_pssm.get_min_range()), max_i_idx);

    // For bidirectional PSSMs, always check both strands (matching max_like_match behavior)
    // For non-bidirectional PSSMs, respect m_strand
    bool check_forward = m_pssm.is_bidirect() || (m_strand != -1);
    bool check_reverse = m_pssm.is_bidirect() || (m_strand != 1);

    int pos = 0;
    for (size_t i = min_i_idx; i <= max_i_idx; ++i) {
        int spat_bin = int(pos / m_spat_bin_size);
        if (spat_bin >= (int)m_spat_log_factors.size()) {
            spat_bin = m_spat_log_factors.size() - 1;
        }
        pos++;

        float spat_log = m_spat_log_factors[spat_bin];

        // Forward strand
        if (check_forward) {
            float logp = 0;
            std::string::const_iterator it = target.begin() + i;
            m_pssm.calc_like(it, logp);
            float val = logp + spat_log;
            if (val > best_val) {
                best_val = val;
                best_index = i;
                best_dir = 1;
            }
        }

        // Reverse strand
        if (check_reverse) {
            float logp_rc = 0;
            std::string::const_iterator it2 = target.begin() + i;
            m_pssm.calc_like_rc(it2, logp_rc);
            float val_rc = logp_rc + spat_log;
            if (val_rc > best_val) {
                best_val = val_rc;
                best_index = i;
                best_dir = -1;
            }
        }
    }

    return compute_position_result(best_index, target.length(), motif_length, best_dir);
}

// Score without spatial weighting
float PWMScorer::score_without_spatial(const std::string& target, int64_t motif_length)
{
    if (m_mode == TOTAL_LIKELIHOOD) {
        float energy;
        m_pssm.integrate_like(target, energy);
        return energy;
    }
    
    if (m_mode == MOTIF_COUNT) {
        return count_motif_hits_no_spatial(target, motif_length);
    }
    
    // MAX_LIKELIHOOD or MAX_LIKELIHOOD_POS
    float best_logp;
    int best_dir;
    bool combine_strands = (m_mode == MAX_LIKELIHOOD);
    std::string::const_iterator best_pos = m_pssm.max_like_match(target, best_logp, best_dir, combine_strands);

    if (m_mode == MAX_LIKELIHOOD) {
        return best_logp;
    }
    
    // MAX_LIKELIHOOD_POS
    size_t pos_idx = best_pos - target.begin();
    return compute_position_result(pos_idx, target.length(), motif_length, best_dir);
}

// Score with spatial weighting
float PWMScorer::score_with_spatial(const std::string& target, int64_t motif_length)
{
    if (m_spat_log_factors.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    if (m_mode == TOTAL_LIKELIHOOD) {
        float energy;
        m_pssm.integrate_energy_logspat(target, energy, m_spat_log_factors, m_spat_bin_size);
        return energy;
    }

    if (m_mode == MAX_LIKELIHOOD) {
        float energy;
        m_pssm.integrate_energy_max_logspat(target, energy, m_spat_log_factors, m_spat_bin_size);
        return energy;
    }

    if (m_mode == MOTIF_COUNT) {
        return count_motif_hits_with_spatial(target, motif_length);
    }

    // MAX_LIKELIHOOD_POS
    return get_max_likelihood_pos_with_spatial(target, motif_length);
}

// Try to advance the sliding window by stride positions
float PWMScorer::try_slide_window(const std::string& target,
                                  const GInterval& original_interval,
                                  const GInterval& expanded_interval,
                                  size_t i_min, size_t i_max, size_t motif_len, size_t stride)
{
    const size_t W = m_slide.window_size;
    const size_t first_incoming_pos = m_slide.pos0 + W;
    const size_t tlen = target.size();
    char strand_mode = m_strand;

    if (m_mode == TOTAL_LIKELIHOOD) {
        // Pop outgoing values
        for (size_t k = 0; k < stride; ++k) {
            m_slide.rlse.pop_front();
        }

        // Push incoming values
        for (size_t k = 0; k < stride; ++k) {
            size_t incoming_i = i_max - (stride - 1 - k);
            if (incoming_i + motif_len > tlen) {
                return std::numeric_limits<float>::quiet_NaN(); // Signal failure
            }
            float val = pos_value_with_spat(m_pssm, target, incoming_i, strand_mode,
                                            get_spatial_log_factor(first_incoming_pos + k));
            m_slide.rlse.push(val);
        }

        // Update cache state
        m_slide.pos0 += stride;
        m_slide.last_interval_start = original_interval.start;
        m_slide.last_interval_end = original_interval.end;
        m_slide.last_i_min = i_min;
        m_slide.last_i_max = i_max;
        m_slide.stride = stride;

        return (float)m_slide.rlse.value();
    }

    if (m_mode == MAX_LIKELIHOOD) {
        // Pop outgoing values - advance base by stride
        m_slide.rmax.pop_front(m_slide.rmax.base_genomic_pos + stride);

        // Push incoming values (log-summed strands, like TOTAL_LIKELIHOOD)
        for (size_t k = 0; k < stride; ++k) {
            size_t incoming_i = i_max - (stride - 1 - k);
            if (incoming_i + motif_len > tlen) {
                return std::numeric_limits<float>::quiet_NaN();
            }
            float val = pos_value_with_spat(m_pssm, target, incoming_i, strand_mode,
                                            get_spatial_log_factor(first_incoming_pos + k));
            int64_t genomic_pos = expanded_interval.start + int64_t(incoming_i);
            m_slide.rmax.push(val, 1, genomic_pos);  // Direction doesn't matter for MAX_LIKELIHOOD
        }

        // Update cache state
        m_slide.pos0 += stride;
        m_slide.last_interval_start = original_interval.start;
        m_slide.last_interval_end = original_interval.end;
        m_slide.last_i_min = i_min;
        m_slide.last_i_max = i_max;
        m_slide.stride = stride;

        return m_slide.rmax.value();
    }

    if (m_mode == MAX_LIKELIHOOD_POS) {
        // Pop outgoing values - advance base by stride
        m_slide.rmax.pop_front(m_slide.rmax.base_genomic_pos + stride);

        // Push incoming values with strand direction and genomic position
        for (size_t k = 0; k < stride; ++k) {
            size_t incoming_i = i_max - (stride - 1 - k);
            if (incoming_i + motif_len > tlen) {
                return std::numeric_limits<float>::quiet_NaN();
            }
            int best_dir = 1;
            float val = pos_value_with_dir(m_pssm, target, incoming_i, strand_mode,
                                           get_spatial_log_factor(first_incoming_pos + k),
                                           best_dir);
            int64_t genomic_pos = expanded_interval.start + int64_t(incoming_i);
            m_slide.rmax.push(val, best_dir, genomic_pos);
        }

        // Update cache state
        m_slide.pos0 += stride;
        m_slide.last_interval_start = original_interval.start;
        m_slide.last_interval_end = original_interval.end;
        m_slide.last_i_min = i_min;
        m_slide.last_i_max = i_max;
        m_slide.stride = stride;

        // MAX_LIKELIHOOD_POS: Get genomic position directly from deque
        int64_t best_genomic_pos = m_slide.rmax.argmax_genomic_position();
        int best_dir = m_slide.rmax.argmax_direction();

        // Convert genomic position to index in target string
        size_t pos_in_target = size_t(best_genomic_pos - expanded_interval.start);

        // Use compute_position_result like the non-sliding implementation
        return compute_position_result(pos_in_target, target.length(), motif_len, best_dir);
    }

    if (m_mode == MOTIF_COUNT) {
        // Pop outgoing values
        for (size_t k = 0; k < stride && !m_slide.hits.empty(); ++k) {
            m_slide.hit_count -= m_slide.hits.front();
            m_slide.hits.pop_front();
        }

        // Push incoming values
        for (size_t k = 0; k < stride; ++k) {
            size_t incoming_i = i_max - (stride - 1 - k);
            if (incoming_i + motif_len > tlen) {
                return std::numeric_limits<float>::quiet_NaN();
            }
            int best_dir = 1;
            float v = pos_value_with_dir(m_pssm, target, incoming_i, strand_mode,
                                         get_spatial_log_factor(first_incoming_pos + k),
                                         best_dir);
            uint8_t hit = (v >= m_score_thresh) ? 1 : 0;
            m_slide.hits.push_back(hit);
            m_slide.hit_count += hit;
        }

        // Update cache state
        m_slide.pos0 += stride;
        m_slide.last_interval_start = original_interval.start;
        m_slide.last_interval_end = original_interval.end;
        m_slide.last_i_min = i_min;
        m_slide.last_i_max = i_max;
        m_slide.stride = stride;

        return (float)m_slide.hit_count;
    }

    return std::numeric_limits<float>::quiet_NaN();
}

// Initialize/seed the sliding window
float PWMScorer::seed_sliding_window(const std::string& target,
                                     const GInterval& original_interval,
                                     const GInterval& expanded_interval,
                                     size_t i_min, size_t i_max, size_t motif_len)
{
    char strand_mode = m_strand;

    // Compute all values for the window
    std::vector<float> vals;
    std::vector<int> dirs;  // strand directions
    std::vector<int64_t> genomic_positions;  // genomic coordinates for each position
    vals.reserve(i_max - i_min + 1);
    dirs.reserve(i_max - i_min + 1);
    genomic_positions.reserve(i_max - i_min + 1);

    size_t pos0 = 0;

    // Use appropriate helper based on mode
    if (m_mode == TOTAL_LIKELIHOOD || m_mode == MAX_LIKELIHOOD) {
        // Use log-sum for TOTAL_LIKELIHOOD and MAX_LIKELIHOOD
        // MAX_LIKELIHOOD finds the max of the combined (log-summed) strand likelihoods
        for (size_t i = i_min; i <= i_max; ++i) {
            float v = pos_value_with_spat(m_pssm, target, i, strand_mode,
                                          get_spatial_log_factor(pos0 + (i - i_min)));
            vals.push_back(v);
            dirs.push_back(1);  // Direction doesn't matter for these modes
            genomic_positions.push_back(expanded_interval.start + int64_t(i));
        }
    } else {
        // Use strand-aware helper for MAX_LIKELIHOOD_POS and MOTIF_COUNT
        // These modes need to track which strand gave the best score
        for (size_t i = i_min; i <= i_max; ++i) {
            int best_dir = 1;
            float v = pos_value_with_dir(m_pssm, target, i, strand_mode,
                                         get_spatial_log_factor(pos0 + (i - i_min)),
                                         best_dir);
            vals.push_back(v);
            dirs.push_back(best_dir);
            genomic_positions.push_back(expanded_interval.start + int64_t(i));
        }
    }

    size_t W = vals.size();

    // Store geometry (using original interval for cache coherency)
    m_slide.valid = true;
    m_slide.chromid = original_interval.chromid;
    m_slide.strand_mode = strand_mode;
    m_slide.last_interval_start = original_interval.start;
    m_slide.last_interval_end = original_interval.end;
    m_slide.last_i_min = i_min;
    m_slide.last_i_max = i_max;
    m_slide.window_size = W;
    m_slide.pos0 = 0;
    m_slide.stride = 0;

    // Seed aggregators based on mode
    if (m_mode == TOTAL_LIKELIHOOD) {
        m_slide.rlse.init(vals);
        return (float)m_slide.rlse.value();
    }

    if (m_mode == MAX_LIKELIHOOD) {
        m_slide.rmax.init(vals, dirs, genomic_positions);
        return m_slide.rmax.value();
    }

    if (m_mode == MAX_LIKELIHOOD_POS) {
        m_slide.rmax.init(vals, dirs, genomic_positions);

        // MAX_LIKELIHOOD_POS: Get genomic position directly from deque
        int64_t best_genomic_pos = m_slide.rmax.argmax_genomic_position();
        int best_dir = m_slide.rmax.argmax_direction();

        // Convert genomic position to index in target string
        size_t pos_in_target = size_t(best_genomic_pos - expanded_interval.start);

        // Use compute_position_result like the non-sliding implementation
        return compute_position_result(pos_in_target, target.length(), motif_len, best_dir);
    }

    if (m_mode == MOTIF_COUNT) {
        m_slide.hits.clear();
        m_slide.hit_count = 0;
        for (float v : vals) {
            uint8_t hit = (v >= m_score_thresh) ? 1 : 0;
            m_slide.hits.push_back(hit);
            m_slide.hit_count += hit;
        }
        return (float)m_slide.hit_count;
    }

    return std::numeric_limits<float>::quiet_NaN();
}

// Try to use sliding window optimization
float PWMScorer::score_with_sliding_window(const std::string& target,
                                           const GInterval& original_interval,
                                           const GInterval& expanded_interval,
                                           size_t i_min, size_t i_max, size_t motif_len)
{
    char strand_mode = m_strand;

    // Calculate stride based on original interval movement
    size_t stride = 0;
    if (m_slide.valid) {
        int64_t step_start = original_interval.start - m_slide.last_interval_start;
        int64_t step_end = original_interval.end - m_slide.last_interval_end;
        if (step_start > 0 && step_start == step_end) {
            stride = static_cast<size_t>(step_start);
        }
    }

    // Check if we can slide
    bool stride_ok = m_slide.valid &&
                     stride > 0 &&
                     stride <= m_slide.window_size &&
                     (m_slide.stride == 0 || stride == m_slide.stride);

    bool can_slide =
        m_slide.valid &&
        m_slide.chromid == original_interval.chromid &&
        m_slide.strand_mode == strand_mode &&
        stride_ok &&
        (i_min == m_slide.last_i_min) &&
        (i_max == m_slide.last_i_max);

    // Try to slide if possible
    if (can_slide) {
        float result = try_slide_window(target, original_interval, expanded_interval, i_min, i_max, motif_len, stride);
        if (!std::isnan(result)) {
            return result;
        }
        // If sliding failed, fall through to re-seeding
    }

    // Seed or re-seed the window
    return seed_sliding_window(target, original_interval, expanded_interval, i_min, i_max, motif_len);
}

float PWMScorer::score_interval(const GInterval& interval, const GenomeChromKey& chromkey)
{
    // Calculate expanded interval to include full motif coverage
    int64_t motif_length = m_pssm.size();
    GInterval expanded_interval = calculate_expanded_interval(interval, chromkey, motif_length);
    expanded_interval.strand = m_strand;

    // Check if interval is too small for motif
    if (!m_extend && (expanded_interval.end - expanded_interval.start) < motif_length) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    try {
        // Read sequence
        std::vector<char> seq;
        m_seqfetch_ptr->read_interval(expanded_interval, chromkey, seq);
        std::string target(seq.begin(), seq.end());

        // Sliding window optimization for contiguous intervals
        // Enabled for all modes when spatial weighting is disabled, except for MAX_LIKELIHOOD_POS
        if (!m_use_spat && m_mode != MAX_LIKELIHOOD_POS) {
            const size_t tlen = target.size();
            const size_t motif_len = m_pssm.size();

            if (tlen >= motif_len) {
                // Calculate allowed start range
                size_t i_min = std::max(0, m_pssm.get_min_range());
                size_t i_max = std::min<size_t>(m_pssm.get_max_range(), tlen - motif_len);

                // Ensure valid range
                if (i_min > i_max) {
                    i_min = 0;
                }

                return score_with_sliding_window(target, interval, expanded_interval, i_min, i_max, motif_len);
            }
        }

        // Standard scoring without sliding window optimization
        if (m_use_spat) {
            return score_with_spatial(target, motif_length);
        } else {
            return score_without_spatial(target, motif_length);
        }

    } catch (TGLException &e) {
        return std::numeric_limits<float>::quiet_NaN();
    }
}

DnaPSSM PWMScorer::create_pssm_from_matrix(SEXP matrix)
{
    DnaPSSM pssm;
    int row_num = Rf_nrows(matrix);
    pssm.resize(row_num);

    // Extract column names from matrix
    SEXP dimnames = Rf_getAttrib(matrix, R_DimNamesSymbol);
    if (Rf_isNull(dimnames) || VECTOR_ELT(dimnames, 1) == R_NilValue) {
        rdb::verror("PWM matrix must have column names 'A', 'C', 'G', and 'T'");
    }

    SEXP colnames = VECTOR_ELT(dimnames, 1);
    if (!Rf_isString(colnames) || Rf_length(colnames) != 4) {
        rdb::verror("PWM matrix must have exactly 4 columns labeled 'A', 'C', 'G', and 'T'");
    }

    // Map column names to indices
    int a_idx = -1, c_idx = -1, g_idx = -1, t_idx = -1;
    for (int i = 0; i < 4; i++) {
        const char *name = CHAR(STRING_ELT(colnames, i));
        if (strcmp(name, "A") == 0) a_idx = i;
        else if (strcmp(name, "C") == 0) c_idx = i;
        else if (strcmp(name, "G") == 0) g_idx = i;
        else if (strcmp(name, "T") == 0) t_idx = i;
    }

    // Verify we found all required indices
    if (a_idx == -1 || c_idx == -1 || g_idx == -1 || t_idx == -1) {
        rdb::verror("PWM matrix must have columns labeled 'A', 'C', 'G', and 'T'");
    }

    double *matrix_ptr = REAL(matrix);

    // Extract probabilities for each position
    for (int i = 0; i < row_num; i++) {
        float pa = matrix_ptr[i + a_idx * row_num];  // Get A probability for position i
        float pc = matrix_ptr[i + c_idx * row_num];  // Get C probability for position i
        float pg = matrix_ptr[i + g_idx * row_num];  // Get G probability for position i
        float pt = matrix_ptr[i + t_idx * row_num];  // Get T probability for position i
        
        pssm[i] = DnaProbVec(pa, pc, pg, pt);
    }

    return pssm;
}
