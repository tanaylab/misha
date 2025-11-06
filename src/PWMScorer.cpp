#include "PWMScorer.h"
#include <algorithm>
#include <limits>
#include <cstring> // For strcmp
#include <cmath>   // For log

// Forward declaration for log_sum_log (defined in util.h)
extern inline void log_sum_log(float& a, float b);

// Helpers that return scores in ORIGINAL-genome strand meaning
// When strand_mode == -1, target is reverse-complemented, so we need to invert the functions
static inline float score_forward_original(const DnaPSSM& pssm,
                                           const std::string& target,
                                           size_t i,
                                           char strand_mode) {
    float s = 0.f;
    auto it = target.begin() + i;
    if (strand_mode == -1) {
        // target is reverse-complemented -> forward(original) == rc(target)
        pssm.calc_like_rc(it, s);
    } else {
        pssm.calc_like(it, s);
    }
    return s;
}

static inline float score_reverse_original(const DnaPSSM& pssm,
                                           const std::string& target,
                                           size_t i,
                                           char strand_mode) {
    float s = 0.f;
    auto it = target.begin() + i;
    if (strand_mode == -1) {
        // target is reverse-complemented -> reverse(original) == forward(target)
        pssm.calc_like(it, s);
    } else {
        pssm.calc_like_rc(it, s);
    }
    return s;
}

// Compute log-likelihood at position i with spatial weighting
// Combines forward and reverse complement strands using logsumexp if bidirectional
// This is used for TOTAL_LIKELIHOOD mode
static inline float pos_value_with_spat(const DnaPSSM& pssm,
                                        const std::string& target,
                                        size_t i,
                                        char strand_mode,      // 1=plus, -1=minus
                                        float spat_log)        // log(spatial_factor)
{
    // Score in ORIGINAL-genome meaning to match non-sliding implementations
    float val = -std::numeric_limits<float>::infinity();
    if (pssm.is_bidirect()) {
        // LSE union of forward & reverse
        float fwd = score_forward_original(pssm, target, i, strand_mode);
        float rev = score_reverse_original(pssm, target, i, strand_mode);
        val = fwd;
        log_sum_log(val, rev);
    } else {
        // Single selected strand
        val = (strand_mode == -1)
              ? score_reverse_original(pssm, target, i, strand_mode)
              : score_forward_original(pssm, target, i, strand_mode);
    }
    return val + spat_log;
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
        float val_f = score_forward_original(pssm, target, i, strand_mode) + spat_log;
        if (val_f > best_val) {
            best_val = val_f;
            best_dir = 1;
        }
    }

    // Evaluate reverse strand
    if (check_reverse) {
        float val_rc = score_reverse_original(pssm, target, i, strand_mode) + spat_log;
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
    m_spat_slide.valid = false;
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

    // Decide which strands to scan:
    // - If bidirect, scan both strands regardless of m_strand
    // - If not bidirect, scan the strand indicated by m_strand (1 or -1)
    const bool scan_both = m_pssm.is_bidirect();
    const bool scan_fwd  = scan_both || (m_strand == 1);
    const bool scan_rev  = scan_both || (m_strand == -1);

    int count = 0;
    size_t max_pos = target.length() - motif_length;

    // Per-position union (LSE) across permitted strands
    for (size_t i = 0; i <= max_pos; ++i) {
        float comb = -std::numeric_limits<float>::infinity();
        bool have = false;
        if (scan_fwd) {
            comb = score_forward_original(m_pssm, target, i, m_strand);
            have = true;
        }
        if (scan_rev) {
            float rc = score_reverse_original(m_pssm, target, i, m_strand);
            if (have) log_sum_log(comb, rc); else { comb = rc; have = true; }
        }
        if (comb >= m_score_thresh) count++;
    }

    return static_cast<float>(count);
}

// Count motif hits with spatial weighting
float PWMScorer::count_motif_hits_with_spatial(const std::string& target, size_t motif_length)
{
    if (target.length() < motif_length) {
        return 0.0f;
    }

    // Decide which strands to scan (same rules as non-spatial)
    const bool scan_both = m_pssm.is_bidirect();
    const bool scan_fwd  = scan_both || (m_strand == 1);
    const bool scan_rev  = scan_both || (m_strand == -1);

    int count = 0;
    size_t max_pos = target.length() - motif_length;

    // Per-position union (LSE) across permitted strands, with spatial offset
    for (size_t i = 0; i <= max_pos; ++i) {
        float spat_log = get_spatial_log_factor(i);
        float comb = -std::numeric_limits<float>::infinity();
        bool have = false;
        if (scan_fwd) {
            comb = score_forward_original(m_pssm, target, i, m_strand);
            have = true;
        }
        if (scan_rev) {
            float rc = score_reverse_original(m_pssm, target, i, m_strand);
            if (have) log_sum_log(comb, rc); else { comb = rc; have = true; }
        }
        if (comb + spat_log >= m_score_thresh) count++;
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
        // With END-only extension policy (correct behavior):
        // Plus strand: new anchors appear at high target indices (close to i_max).
        // Minus strand: target is reverse-complemented, so new anchors appear at low indices (close to i_min).
        //
        // Difference is in deque management due to RC reversal:
        // Plus strand: target forward, old at front, new at back → pop front, push back
        // Minus strand: target RC'd, old at back, new at front → pop back, push front

        if (m_strand == -1) {
            // Minus strand: pop from back (lowest genomic starts), push new anchors to front.
            for (size_t k = 0; k < stride; ++k) {
                m_slide.rlse.pop_back();
            }

            // New anchors on minus strand reside at the lowest target indices (i_min, i_min+1, ...),
            // so push them in reverse order to keep the deque in ascending target-index order.
            for (size_t k = stride; k > 0; --k) {
                size_t incoming_i = i_min + (k - 1);
                if (incoming_i + motif_len > tlen) {
                    return std::numeric_limits<float>::quiet_NaN();
                }
                float val = pos_value_with_spat(m_pssm, target, incoming_i, strand_mode,
                                                get_spatial_log_factor(first_incoming_pos + (stride - k)));
                m_slide.rlse.push_front(val);
            }
        } else {
            // Plus strand: pop from front, push to back
            for (size_t k = 0; k < stride; ++k) {
                m_slide.rlse.pop_front();
            }

            // Push incoming values from the right side
            for (size_t k = 0; k < stride; ++k) {
                size_t incoming_i = i_max - (stride - 1 - k);
                if (incoming_i + motif_len > tlen) {
                    return std::numeric_limits<float>::quiet_NaN();
                }
                float val = pos_value_with_spat(m_pssm, target, incoming_i, strand_mode,
                                                get_spatial_log_factor(first_incoming_pos + k));
                m_slide.rlse.push(val);
            }
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
        // Uses genomic positions, so same for both strands
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
        bool is_minus = (strand_mode == -1);

        const bool scan_both = m_pssm.is_bidirect();
        const bool scan_fwd  = scan_both || (strand_mode == 1);
        const bool scan_rev  = scan_both || (strand_mode == -1);

        if (is_minus) {
            // Minus strand: pop from back, push to front
            for (size_t k = 0; k < stride && !m_slide.hits.empty(); ++k) {
                m_slide.hit_count -= m_slide.hits.back();
                m_slide.hits.pop_back();
            }

            // Push incoming union hits at front (new minus-strand anchors are near i_min)
            for (size_t k = stride; k > 0; --k) {
                size_t incoming_i = i_min + (k - 1);
                if (incoming_i + motif_len > tlen) {
                    return std::numeric_limits<float>::quiet_NaN();
                }
                float spat_log = get_spatial_log_factor(first_incoming_pos + (stride - k));

                float comb = -std::numeric_limits<float>::infinity();
                bool have = false;
                if (scan_fwd) {
                    float f = score_forward_original(m_pssm, target, incoming_i, strand_mode);
                    comb = f; have = true;
                }
                if (scan_rev) {
                    float rc = score_reverse_original(m_pssm, target, incoming_i, strand_mode);
                    if (have) log_sum_log(comb, rc); else { comb = rc; have = true; }
                }

                uint8_t hit = ((comb + spat_log) >= m_score_thresh) ? 1 : 0;
                m_slide.hits.push_front(hit);
                m_slide.hit_count += hit;
            }
        } else {
            // Plus strand: pop from front, push to back
            for (size_t k = 0; k < stride && !m_slide.hits.empty(); ++k) {
                m_slide.hit_count -= m_slide.hits.front();
                m_slide.hits.pop_front();
            }

            // Push incoming union hits (0/1 per position)
            for (size_t k = 0; k < stride; ++k) {
                size_t incoming_i = i_max - (stride - 1 - k);
                if (incoming_i + motif_len > tlen) {
                    return std::numeric_limits<float>::quiet_NaN();
                }
                float spat_log = get_spatial_log_factor(first_incoming_pos + k);

                float comb = -std::numeric_limits<float>::infinity();
                bool have = false;
                if (scan_fwd) {
                    float f = score_forward_original(m_pssm, target, incoming_i, strand_mode);
                    comb = f; have = true;
                }
                if (scan_rev) {
                    float rc = score_reverse_original(m_pssm, target, incoming_i, strand_mode);
                    if (have) log_sum_log(comb, rc); else { comb = rc; have = true; }
                }

                uint8_t hit = ((comb + spat_log) >= m_score_thresh) ? 1 : 0;
                m_slide.hits.push_back(hit);
                m_slide.hit_count += hit;
            }
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
    } else if (m_mode == MAX_LIKELIHOOD_POS) {
        // Need best score and direction per position
        for (size_t i = i_min; i <= i_max; ++i) {
            int best_dir = 1;
            float v = pos_value_with_dir(m_pssm, target, i, strand_mode,
                                         get_spatial_log_factor(pos0 + (i - i_min)),
                                         best_dir);
            vals.push_back(v);
            dirs.push_back(best_dir);
            genomic_positions.push_back(expanded_interval.start + int64_t(i));
        }
    } else { // MOTIF_COUNT: handled separately below (per-position union)
        // no-op: we don't need vals/dirs here
    }

    size_t W = (m_mode == MOTIF_COUNT) ? (i_max - i_min + 1) : vals.size();

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

        const bool scan_both = m_pssm.is_bidirect();
        const bool scan_fwd  = scan_both || (strand_mode == 1);
        const bool scan_rev  = scan_both || (strand_mode == -1);

        size_t pos0_local = 0;
        for (size_t i = i_min; i <= i_max; ++i, ++pos0_local) {
            float spat_log = get_spatial_log_factor(pos0_local);

            float comb = -std::numeric_limits<float>::infinity();
            bool have = false;
            if (scan_fwd) {
                float f = score_forward_original(m_pssm, target, i, strand_mode);
                comb = f; have = true;
            }
            if (scan_rev) {
                float rc = score_reverse_original(m_pssm, target, i, strand_mode);
                if (have) log_sum_log(comb, rc); else { comb = rc; have = true; }
            }

            uint8_t hit = ((comb + spat_log) >= m_score_thresh) ? 1 : 0;
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

        const size_t tlen = target.size();
        const size_t motif_len = m_pssm.size();

        if (tlen >= motif_len) {
            // Calculate allowed start range
            size_t i_min = std::max(0, m_pssm.get_min_range());
            size_t i_max = std::min<size_t>(m_pssm.get_max_range(), tlen - motif_len);

            // Clamp the scanning window to anchors whose starts fall inside the iterator
            const int64_t interval_len = interval.end - interval.start;
            if (interval_len > 0) {
                const int64_t extra_left = std::max<int64_t>(0, interval.start - expanded_interval.start);
                const int64_t max_valid = static_cast<int64_t>(tlen - motif_len);

                auto clamp_index = [&](int64_t idx) -> size_t {
                    if (idx < 0) {
                        return 0;
                    }
                    if (idx > max_valid) {
                        return static_cast<size_t>(std::max<int64_t>(0, max_valid));
                    }
                    return static_cast<size_t>(idx);
                };

                // For all strand modes, motif anchors start at extra_left in the fetched sequence
                // (since we extend END only, extra_left is always 0 when extend=true)
                int64_t desired_min = extra_left;
                int64_t desired_max = desired_min + interval_len - 1;

                size_t clamped_min = clamp_index(desired_min);
                size_t clamped_max = clamp_index(desired_max);

                if (clamped_min <= clamped_max) {
                    i_min = std::max(i_min, clamped_min);
                    i_max = std::min(i_max, clamped_max);
                } else {
                    i_min = clamped_min;
                    i_max = clamped_min;
                }
            }

            // Ensure valid range
            if (i_min > i_max) {
                i_min = 0;
            }

            // Try spatial sliding window optimization
            if (m_use_spat) {
                size_t stride = 0;
                if (can_use_spatial_sliding(interval, expanded_interval, i_min, i_max, motif_len, stride)) {
                    if (stride == 0) {
                        // Need to seed
                        spat_seed(target, expanded_interval, i_min, i_max, motif_len);
                    } else {
                        // Can slide - loop stride times to handle stride>1
                        for (size_t s = 0; s < stride; ++s) {
                            spat_slide_once(target, expanded_interval, i_min, i_max, motif_len);
                            if (!m_spat_slide.valid) {
                                // Sliding failed, reseed
                                spat_seed(target, expanded_interval, i_min, i_max, motif_len);
                                break;
                            }
                        }
                    }

                    // Update cache state
                    m_spat_slide.last_interval_start = interval.start;
                    m_spat_slide.last_interval_end = interval.end;
                    m_spat_slide.last_i_min = i_min;
                    m_spat_slide.last_i_max = i_max;

                    // Return answer based on mode
                    switch (m_mode) {
                        case TOTAL_LIKELIHOOD:   return spat_answer_TOTAL();
                        case MAX_LIKELIHOOD:     return spat_answer_MAX();
                        case MAX_LIKELIHOOD_POS: return spat_answer_MAXPOS(target, expanded_interval, motif_len, i_min);
                        case MOTIF_COUNT:        return spat_answer_COUNT();
                        default: break;
                    }
                }

                // Fallback to standard spatial scoring
                return score_with_spatial(target, motif_length);
            }

            // Non-spatial sliding window optimization (original code)
            if (!m_use_spat && m_mode != MAX_LIKELIHOOD_POS) {
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

// ============================================================================
// Spatial Sliding Window Implementation
// ============================================================================

// Ring buffer and bin helpers
inline size_t PWMScorer::ring_idx_from_j(size_t j) const {
    return (m_spat_slide.head + j) % m_spat_slide.W;
}

inline size_t PWMScorer::j_from_ring_idx(size_t ridx) const {
    const size_t W = m_spat_slide.W;
    return (ridx + W - m_spat_slide.head) % W;
}

// Map j to spatial bin index, clamped to last bin (bins-1)
inline size_t PWMScorer::bin_of_j(size_t j) const {
    size_t b = j / (size_t)m_spat_slide.B;
    if (b >= m_spat_slide.bins) b = m_spat_slide.bins - 1;
    return b;
}

// Compute motif(s) at target offset `i_in_target` (0-based in expanded interval)
void PWMScorer::compute_motif_at(const std::string& target, size_t i_in_target,
                                 float& fwd, float& rc, uint8_t& has_f, uint8_t& has_r) {
    has_f = 0; has_r = 0;
    fwd = -std::numeric_limits<float>::infinity();
    rc = -std::numeric_limits<float>::infinity();

    // Respect strand / bidirect rules
    const bool check_fwd = (m_pssm.is_bidirect() || m_strand != -1);
    const bool check_rc  = (m_pssm.is_bidirect() || m_strand != +1);

    if (check_fwd) {
        fwd = score_forward_original(m_pssm, target, i_in_target, m_strand);
        has_f = 1;
    }
    if (check_rc) {
        rc = score_reverse_original(m_pssm, target, i_in_target, m_strand);
        has_r = 1;
    }
}

// COUNT helper: is this motif score a hit in bin b?
inline bool PWMScorer::is_hit(float m, size_t b) const {
    if (!std::isfinite(m)) return false;
    const float Lb = m_spat_log_factors[b];
    return (m + Lb) >= m_score_thresh;
}

// TOTAL per-bin LSE maintenance
void PWMScorer::bin_total_recompute(size_t b) {
    SpatSlideCache& S = m_spat_slide;
    const size_t j0 = b * (size_t)S.B;
    const size_t j1 = (b == S.bins - 1) ? (S.W - 1) : (std::min(S.W, (b+1)*(size_t)S.B) - 1);

    double a = -std::numeric_limits<double>::infinity();
    // find joint anchor (max over all strands present in bin)
    for (size_t j = j0; j <= j1; ++j) {
        const size_t ridx = ring_idx_from_j(j);
        if (S.has_fwd[ridx]) a = std::max(a, (double)S.motif_fwd[ridx]);
        if (S.has_rc[ridx])  a = std::max(a, (double)S.motif_rc[ridx]);
    }
    if (!std::isfinite(a)) a = -std::numeric_limits<double>::infinity(); // empty bin

    double sf = 0.0, sr = 0.0;
    if (std::isfinite(a)) {
        for (size_t j = j0; j <= j1; ++j) {
            const size_t ridx = ring_idx_from_j(j);
            if (S.has_fwd[ridx]) sf += std::exp((double)S.motif_fwd[ridx] - a);
            if (S.has_rc[ridx])  sr += std::exp((double)S.motif_rc[ridx]  - a);
        }
    }
    S.bin_anchor[b] = a;
    S.bin_sum_fwd[b] = sf;
    S.bin_sum_rc[b]  = sr;
    S.bin_dirty[b] = 0;
}

// Add a single motif value to bin b; rebase if it becomes new anchor
void PWMScorer::bin_total_add(size_t b, float m, bool is_rc) {
    SpatSlideCache& S = m_spat_slide;
    double& a  = S.bin_anchor[b];
    double& sf = S.bin_sum_fwd[b];
    double& sr = S.bin_sum_rc[b];
    if (!std::isfinite(m)) return;

    if (!std::isfinite(a)) {
        a = m;
        if (is_rc) sr = 1.0; else sf = 1.0;
        return;
    }

    if ((double)m <= a) {
        double e = std::exp((double)m - a);
        if (is_rc) sr += e; else sf += e;
    } else {
        // new anchor -> rescale
        double scale = std::exp(a - (double)m);
        sf *= scale; sr *= scale; a = m;
        if (is_rc) sr += 1.0; else sf += 1.0;
    }
}

// Remove a single motif value from bin b
void PWMScorer::bin_total_remove(size_t b, float m, bool is_rc) {
    SpatSlideCache& S = m_spat_slide;
    if (!std::isfinite(m)) return;
    double& a  = S.bin_anchor[b];
    double& sf = S.bin_sum_fwd[b];
    double& sr = S.bin_sum_rc[b];

    if (!std::isfinite(a)) return;

    const double e = std::exp((double)m - a);
    if (is_rc) { sr -= e; if (sr < 0) sr = 0; }
    else       { sf -= e; if (sf < 0) sf = 0; }

    // If we removed (one of) the anchor elements, sums may get tiny/inaccurate — mark dirty.
    if (e > 0.5 || sf + sr < 1e-12) S.bin_dirty[b] = 1;
}

// MAX / MAX_POS per-bin maintenance
void PWMScorer::bin_max_maybe_recompute(size_t b) {
    SpatSlideCache& S = m_spat_slide;
    const size_t j0 = b * (size_t)S.B;
    const size_t j1 = (b == S.bins - 1) ? (S.W - 1) : (std::min(S.W, (b+1)*(size_t)S.B) - 1);

    float best = -std::numeric_limits<float>::infinity();
    int best_idx = -1;
    int best_dir = +1;

    // For MAX_LIKELIHOOD with bidirectional, use LSE of strands (matching legacy behavior)
    // For MAX_LIKELIHOOD_POS, pick the best single strand to track direction
    const bool use_lse = (m_mode == MAX_LIKELIHOOD && m_pssm.is_bidirect());

    for (size_t j = j0; j <= j1; ++j) {
        const size_t ridx = ring_idx_from_j(j);

        if (use_lse && S.has_fwd[ridx] && S.has_rc[ridx]) {
            // Combine forward and reverse using log-sum-exp
            float combined = S.motif_fwd[ridx];
            log_sum_log(combined, S.motif_rc[ridx]);
            if (combined > best) {
                best = combined;
                best_idx = (int)ridx;
                best_dir = +1;  // Direction doesn't matter for MAX_LIKELIHOOD
            }
        } else {
            // MAX_POS or single strand: pick the best strand
            if (S.has_fwd[ridx] && S.motif_fwd[ridx] > best) {
                best = S.motif_fwd[ridx];
                best_idx = (int)ridx;
                best_dir = +1;
            }
            if (S.has_rc[ridx] && S.motif_rc[ridx] > best) {
                best = S.motif_rc[ridx];
                best_idx = (int)ridx;
                best_dir = -1;
            }
        }
    }
    S.bin_max[b] = { best, best_idx, best_dir };
}

void PWMScorer::bin_max_consider(size_t b, float m, int ridx, int dir) {
    SpatSlideCache& S = m_spat_slide;
    if (m > S.bin_max[b].val) {
        S.bin_max[b] = { m, ridx, dir };
    }
}

// Initialize the spatial sliding cache
void PWMScorer::spat_seed(const std::string& target, const GInterval& expd,
                          size_t i_min, size_t i_max, size_t motif_len) {
    SpatSlideCache& S = m_spat_slide;
    S.valid = false;
    S.head = 0;

    S.W = i_max - i_min + 1;
    S.B = m_spat_bin_size;
    const size_t spat_len = m_spat_log_factors.size();
    S.bins = std::min((size_t)std::ceil((double)S.W / (double)S.B), spat_len);

    // Allocate contiguous buffers
    S.motif_fwd.assign(S.W, -std::numeric_limits<float>::infinity());
    S.motif_rc.assign(S.W, -std::numeric_limits<float>::infinity());
    S.has_fwd.assign(S.W, 0);
    S.has_rc.assign(S.W, 0);

    S.bin_anchor.assign(S.bins, -std::numeric_limits<double>::infinity());
    S.bin_sum_fwd.assign(S.bins, 0.0);
    S.bin_sum_rc.assign(S.bins, 0.0);
    S.bin_dirty.assign(S.bins, 0);

    S.bin_max.assign(S.bins, SpatSlideCache::BinMax{});
    S.bin_hits_fwd.assign(S.bins, 0);
    S.bin_hits_rc.assign(S.bins, 0);

    // Fill motifs at j=0..W-1
    for (size_t j = 0; j < S.W; ++j) {
        const size_t i = i_min + j; // offset inside target/expanded
        if (i + motif_len > target.size()) {
            break; // Safety check
        }
        float fwd, rc;
        uint8_t hf, hr;
        compute_motif_at(target, i, fwd, rc, hf, hr);
        S.motif_fwd[j] = fwd;
        S.motif_rc[j] = rc;
        S.has_fwd[j] = hf;
        S.has_rc[j] = hr;
    }

    // Build per-bin aggregates
    for (size_t b = 0; b < S.bins; ++b) {
        bin_total_recompute(b);      // sets anchor/sums from the current j-range in bin b
        bin_max_maybe_recompute(b);  // sets max/idx/dir for bin b

        // COUNT: compute initial counts
        int cf = 0, cr = 0;
        const size_t j_start = b * (size_t)S.B;
        const size_t j_end   = std::min(S.W, (b+1)*(size_t)S.B);
        const float Lb = m_spat_log_factors[std::min(b, S.bins-1)];
        for (size_t j = j_start; j < j_end; ++j) {
            if (S.has_fwd[j] && (S.motif_fwd[j] + Lb) >= m_score_thresh) ++cf;
            if (S.has_rc[j]  && (S.motif_rc[j]  + Lb) >= m_score_thresh) ++cr;
        }
        S.bin_hits_fwd[b] = cf;
        S.bin_hits_rc[b] = cr;
    }

    S.strand_mode = m_strand;
    S.chromid = expd.chromid;
    S.valid = true;
}

// Slide the spatial window by one position
void PWMScorer::spat_slide_once(const std::string& target, const GInterval& expd,
                                size_t i_min, size_t i_max, size_t motif_len) {
    SpatSlideCache& S = m_spat_slide;
    const size_t W = S.W;
    const int B = S.B;

    // Strand-specific sliding direction:
    // Plus strand: target is forward, sliding forward means j=0 goes out, j=W-1 comes in
    // Minus strand: target is RC'd, sliding forward means j=W-1 goes out, j=0 comes in
    const bool is_minus = (m_strand == -1);

    // 1) OUTGOING position
    const size_t j_out = is_minus ? (W - 1) : 0;
    {
        const size_t ridx_out = ring_idx_from_j(j_out);
        const size_t b_out = bin_of_j(j_out);

        // TOTAL remove
        if (S.has_fwd[ridx_out]) bin_total_remove(b_out, S.motif_fwd[ridx_out], false);
        if (S.has_rc[ridx_out])  bin_total_remove(b_out, S.motif_rc[ridx_out],  true);

        // MAX remove if outgoing was argmax (needs re-scan)
        if (S.bin_max[b_out].idx == (int)ridx_out) {
            bin_max_maybe_recompute(b_out);
        }

        // COUNT decrement hits if any (use L_b)
        const float Lb = m_spat_log_factors[b_out];
        if (S.has_fwd[ridx_out] && (S.motif_fwd[ridx_out] + Lb) >= m_score_thresh) {
            S.bin_hits_fwd[b_out]--;
        }
        if (S.has_rc[ridx_out] && (S.motif_rc[ridx_out] + Lb) >= m_score_thresh) {
            S.bin_hits_rc[b_out]--;
        }
    }

    // 2) BOUNDARY MOVERS
    // Plus strand: When j→j-1 (sliding forward), boundaries at j=B, 2B, 3B, ... cross to adjacent bin
    // Minus strand: When j→j+1 (sliding forward with RC), boundaries at j=B-1, 2B-1, ... cross to adjacent bin
    // Restrict to respect clamping
    const size_t max_mover_j = std::min(W, S.bins * (size_t)B);

    // Collect boundary positions that will cross bins after head advance
    std::vector<size_t> boundary_j_vals;
    if (is_minus) {
        // For minus strand, check j=B-1, 2B-1, 3B-1, ... (they become j=B, 2B, 3B after head decrement)
        for (size_t j_old = (size_t)B - 1; j_old < max_mover_j && j_old < W; j_old += (size_t)B) {
            // Skip the outgoing position (j_out = W-1) to avoid double-processing
            if (j_old + 1 == W) {
                break;
            }
            boundary_j_vals.push_back(j_old);
        }
    } else {
        // For plus strand, check j=B, 2B, 3B, ... (they become j=B-1, 2B-1, 3B-1 after head increment)
        for (size_t j_old = (size_t)B; j_old < max_mover_j; j_old += (size_t)B) {
            boundary_j_vals.push_back(j_old);
        }
    }

    for (size_t j_old : boundary_j_vals) {
        const size_t ridx = ring_idx_from_j(j_old);
        const size_t b_old = bin_of_j(j_old);
        const size_t j_new = is_minus ? ((j_old + 1) % W) : ((j_old + W - 1) % W);
        const size_t b_new = bin_of_j(j_new);

        if (b_new == b_old) continue; // no bin change

        // move TOTAL contributions: old bin -> new bin
        if (S.has_fwd[ridx]) {
            bin_total_remove(b_old, S.motif_fwd[ridx], false);
            bin_total_add(b_new,  S.motif_fwd[ridx], false);
        }
        if (S.has_rc[ridx]) {
            bin_total_remove(b_old, S.motif_rc[ridx],  true);
            bin_total_add(b_new,  S.motif_rc[ridx],  true);
        }

        // move MAX participation
        if (S.bin_max[b_old].idx == (int)ridx) {
            bin_max_maybe_recompute(b_old);
        }

        // consider as candidate in new bin (motif-only)
        // For MAX_LIKELIHOOD with bidirectional, use LSE; for MAX_POS, consider each strand
        const bool use_lse = (m_mode == MAX_LIKELIHOOD && m_pssm.is_bidirect());
        if (use_lse && S.has_fwd[ridx] && S.has_rc[ridx]) {
            float combined = S.motif_fwd[ridx];
            log_sum_log(combined, S.motif_rc[ridx]);
            bin_max_consider(b_new, combined, (int)ridx, +1);
        } else {
            if (S.has_fwd[ridx]) bin_max_consider(b_new, S.motif_fwd[ridx], (int)ridx, +1);
            if (S.has_rc[ridx])  bin_max_consider(b_new, S.motif_rc[ridx],  (int)ridx, -1);
        }

        // COUNT: re-evaluate hit against new bin threshold
        const float L_old = m_spat_log_factors[b_old];
        const float L_new = m_spat_log_factors[b_new];
        if (S.has_fwd[ridx]) {
            const bool was = (S.motif_fwd[ridx] + L_old) >= m_score_thresh;
            const bool now = (S.motif_fwd[ridx] + L_new) >= m_score_thresh;
            if (was && !now) S.bin_hits_fwd[b_old]--;
            else if (!was && now) S.bin_hits_fwd[b_new]++;
            else if (was && now) { S.bin_hits_fwd[b_old]--; S.bin_hits_fwd[b_new]++; }
        }
        if (S.has_rc[ridx]) {
            const bool was = (S.motif_rc[ridx] + L_old) >= m_score_thresh;
            const bool now = (S.motif_rc[ridx] + L_new) >= m_score_thresh;
            if (was && !now) S.bin_hits_rc[b_old]--;
            else if (!was && now) S.bin_hits_rc[b_new]++;
            else if (was && now) { S.bin_hits_rc[b_old]--; S.bin_hits_rc[b_new]++; }
        }
    }

    // 3) ADVANCE HEAD (j relabel)
    // Plus strand: j→j-1 means head advances forward (increment)
    // Minus strand: j→j+1 means head advances backward (decrement) due to RC reversal
    if (is_minus) {
        S.head = (S.head + W - 1) % W;  // decrement
    } else {
        S.head = (S.head + 1) % W;       // increment
    }

    // 4) INCOMING position
    const size_t j_in = is_minus ? 0 : (W - 1);
    {
        const size_t ridx_in = ring_idx_from_j(j_in);
        const size_t i_in_target = i_min + j_in; // new sequence offset
        if (i_in_target + motif_len > target.size()) {
            // Safety: mark invalid if we're out of bounds
            S.valid = false;
            return;
        }

        float fwd, rc;
        uint8_t hf, hr;
        compute_motif_at(target, i_in_target, fwd, rc, hf, hr);

        S.motif_fwd[ridx_in] = fwd;
        S.has_fwd[ridx_in] = hf;
        S.motif_rc[ridx_in]  = rc;
        S.has_rc[ridx_in]  = hr;

        const size_t b_in = bin_of_j(j_in);

        // TOTAL add
        if (hf) bin_total_add(b_in, fwd, false);
        if (hr) bin_total_add(b_in, rc,  true);

        // MAX consider
        // For MAX_LIKELIHOOD with bidirectional, use LSE; for MAX_POS, consider each strand
        const bool use_lse = (m_mode == MAX_LIKELIHOOD && m_pssm.is_bidirect());
        if (use_lse && hf && hr) {
            float combined = fwd;
            log_sum_log(combined, rc);
            bin_max_consider(b_in, combined, (int)ridx_in, +1);
        } else {
            if (hf) bin_max_consider(b_in, fwd, (int)ridx_in, +1);
            if (hr) bin_max_consider(b_in, rc,  (int)ridx_in, -1);
        }

        // COUNT add
        const float Lb = m_spat_log_factors[b_in];
        if (hf && (fwd + Lb) >= m_score_thresh) S.bin_hits_fwd[b_in]++;
        if (hr && (rc  + Lb) >= m_score_thresh) S.bin_hits_rc[b_in]++;
    }
}

// Answer functions
float PWMScorer::spat_answer_TOTAL() {
    SpatSlideCache& S = m_spat_slide;

    // Recompute dirty bins (rare)
    for (size_t b = 0; b < S.bins; ++b) {
        if (S.bin_dirty[b]) bin_total_recompute(b);
    }

    // A = max_b (a_b + L_b)
    double A = -std::numeric_limits<double>::infinity();
    for (size_t b = 0; b < S.bins; ++b) {
        const double Lb = (double)m_spat_log_factors[b];
        A = std::max(A, S.bin_anchor[b] + Lb);
    }
    if (!std::isfinite(A)) return -std::numeric_limits<float>::infinity();

    // sum_b exp(a_b + L_b - A) * (s_f + s_r)
    double Ssum = 0.0;
    for (size_t b = 0; b < S.bins; ++b) {
        const double Lb = (double)m_spat_log_factors[b];
        const double w  = std::exp( (S.bin_anchor[b] + Lb) - A );
        Ssum += w * (S.bin_sum_fwd[b] + S.bin_sum_rc[b]);
    }
    if (Ssum <= 0.0) return (float)A; // degenerate
    return (float)(A + std::log(Ssum));
}

float PWMScorer::spat_answer_MAX() {
    SpatSlideCache& S = m_spat_slide;
    float best = -std::numeric_limits<float>::infinity();
    for (size_t b = 0; b < S.bins; ++b) {
        const float cand = S.bin_max[b].val + m_spat_log_factors[b];
        if (cand > best) best = cand;
    }
    return best;
}

float PWMScorer::spat_answer_MAXPOS(const std::string& target,
                                    const GInterval& expd, size_t motif_len, size_t i_min) {
    SpatSlideCache& S = m_spat_slide;

    // Find best bin then best overall (include spatial log)
    float best_val = -std::numeric_limits<float>::infinity();
    int best_idx = -1;
    int best_dir = +1;
    for (size_t b = 0; b < S.bins; ++b) {
        const float cand = S.bin_max[b].val + m_spat_log_factors[b];
        if (cand > best_val) {
            best_val = cand;
            best_idx = S.bin_max[b].idx;
            best_dir = S.bin_max[b].dir;
        }
    }
    if (best_idx < 0) return std::numeric_limits<float>::quiet_NaN();

    // Convert ring index to relative j, then to absolute target index
    const size_t j = j_from_ring_idx((size_t)best_idx);
    const size_t target_idx = i_min + j;

    // Use existing position computation
    return compute_position_result(target_idx, target.length(), motif_len, best_dir);
}

float PWMScorer::spat_answer_COUNT() const {
    const SpatSlideCache& S = m_spat_slide;
    int total = 0;
    for (size_t b = 0; b < S.bins; ++b) {
        total += (S.bin_hits_fwd[b] + S.bin_hits_rc[b]);
    }
    return (float)total;
}

// Check if we can use spatial sliding window optimization
bool PWMScorer::can_use_spatial_sliding(const GInterval& orig, const GInterval& expd,
                                        size_t i_min, size_t i_max, size_t motif_len,
                                        size_t& stride) const {
    // Must have spatial enabled
    if (!m_use_spat) return false;

    // Mode must be supported
    if (m_mode != TOTAL_LIKELIHOOD && m_mode != MAX_LIKELIHOOD &&
        m_mode != MAX_LIKELIHOOD_POS && m_mode != MOTIF_COUNT) {
        return false;
    }

    // Check for environment variable to disable
    const char* disable_env = std::getenv("MISHA_DISABLE_SPATIAL_SLIDING");
    if (disable_env && std::string(disable_env) == "1") {
        return false;
    }

    // Minus strand now supported with reverse sliding direction
    // (no restriction needed)

    // MAX modes now correctly use LSE for bidirectional motifs
    // (no restriction needed)

    const size_t W = i_max - i_min + 1;
    if (W == 0) return false;

    // If not valid, we can't slide (but we can seed)
    if (!m_spat_slide.valid) {
        stride = 0;
        return true; // we'll seed
    }

    // Check geometry consistency
    if (m_spat_slide.chromid != orig.chromid) return false;
    if (m_spat_slide.strand_mode != m_strand) return false;
    if (m_spat_slide.W != W) return false;

    // Calculate stride from original interval movement
    int64_t step_start = orig.start - m_spat_slide.last_interval_start;
    int64_t step_end = orig.end - m_spat_slide.last_interval_end;

    // Must be consistent forward movement
    if (step_start != step_end || step_start <= 0) return false;

    stride = static_cast<size_t>(step_start);

    // Stride must be reasonable (not larger than window)
    if (stride > W) return false;

    // Limit maximum stride to avoid excessive looping.
    // This is a heuristic: for very large strides (e.g., >100bp), the cost of
    // looping `spat_slide_once` `stride` times can exceed the cost of a full re-seed.
    // Falling back to reseeding is safer and often faster in such cases.
    if (stride > 100) return false;

    // Check i_min, i_max consistency
    if (i_min != m_spat_slide.last_i_min || i_max != m_spat_slide.last_i_max) {
        return false;
    }

    return true;
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
