#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include "GlmVarProcessor.h"
#include "TrackExpressionVars.h"
#include "GenomeTrack1D.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "rdbutils.h"

using namespace rdb;
using namespace std;

void GlmVarProcessor::prepare_batch(TrackExpressionVars::Track_vars &track_vars)
{
    m_glm_vars.clear();
    for (auto &var : track_vars) {
        if (var.val_func == TrackExpressionVars::Track_var::GLM_PREDICT) {
            m_glm_vars.push_back(&var);
        }
    }
}

void GlmVarProcessor::process_glm_vars(
    TrackExpressionVars::Track_vars &track_vars,
    const GInterval &interval,
    unsigned idx)
{
    for (auto *var : m_glm_vars) {
        process_single_glm_var(*var, interval, idx);
    }
}

void GlmVarProcessor::process_single_glm_var(
    TrackExpressionVars::Track_var &var,
    const GInterval &interval,
    unsigned idx)
{
    using GlmEntry = TrackExpressionVars::Track_var::GlmEntry;
    using GlmScalingGroup = TrackExpressionVars::Track_var::GlmScalingGroup;

    int64_t center = (interval.start + interval.end) / 2;

    // ---- Multi-selector: compound bin = sum_m b_m * stride_m, NaN propagation ----
    int b = 0;
    const int nsel = (int)var.glm_selector_fixedbins.size();
    bool selector_failed = false;
    for (int m_sel = 0; m_sel < nsel; ++m_sel) {
        GenomeTrackFixedBin *sel = var.glm_selector_fixedbins[m_sel];
        unsigned bs = var.glm_selector_bin_sizes[m_sel];
        if (!sel || bs == 0) { selector_failed = true; break; }
        int64_t sel_bin = center / (int64_t)bs;
        int64_t raw_count;
        const float *ptr = sel->get_mmap_bins_ptr(sel_bin, 1, raw_count);
        if (!ptr || raw_count <= 0 || !std::isfinite(ptr[0])) { selector_failed = true; break; }
        int bm = var.glm_selector_binfinders[m_sel].val2bin((double)ptr[0]);
        if (bm < 0) { selector_failed = true; break; }
        b += bm * (int)var.glm_selector_strides[m_sel];
    }
    if (selector_failed) {
        var.var[idx] = NAN;
        return;
    }
    // For nsel == 0 (no selector), b stays 0 — valid since glm_num_bins == 1.

    double result = var.glm_bias[b];
    bool has_interactions = !var.glm_interactions.empty();
    bool has_kernel_bins = !var.glm_kernel_bins.empty();
    int N = (int)var.glm_entries.size();
    int M = (int)var.glm_interactions.size();

    // Bin-major layout: stride-1 access across entries within the same bin
    const double *bin_weights = var.glm_all_weights.data() + (int64_t)b * N;
    const double *bin_inter_weights = M > 0 ? var.glm_all_inter_weights.data() + (int64_t)b * M : nullptr;

    // ---- Pre-compute interaction-referenced entries (once per var, reuse across positions) ----
    const void *var_interactions_id = &var.glm_interactions;
    if (has_interactions && (m_inter_referenced.size() != (size_t)N || m_inter_referenced_owner != var_interactions_id)) {
        m_inter_referenced.assign(N, false);
        for (const auto &inter : var.glm_interactions) {
            m_inter_referenced[inter.entry_i] = true;
            m_inter_referenced[inter.entry_j] = true;
        }
        m_inter_referenced_owner = var_interactions_id;
    }

    // ---- Main effects ----
    // Track groups provide cache-friendly ordering: read one track's super-window,
    // then process all its scaling groups. Each scaling group aggregates+scales once.
    if (!var.glm_track_groups.empty()) {
        using GlmTrackGroup = TrackExpressionVars::Track_var::GlmTrackGroup;

        for (const GlmTrackGroup &tg : var.glm_track_groups) {

            // Per-bin skip: if all entries in this track group have zero weight
            // and none are interaction-referenced, skip the track read entirely.
            if (var.glm_num_bins > 1) {
                bool tg_needed = false;
                for (int sg_idx : tg.scaling_group_indices) {
                    for (int eidx : var.glm_scaling_groups[sg_idx].entry_indices) {
                        if (bin_weights[var.glm_entries[eidx].weight_offset] != 0.0 ||
                            (has_interactions && m_inter_referenced[eidx])) {
                            tg_needed = true;
                            goto tg_check_done;
                        }
                    }
                }
                tg_check_done:
                if (!tg_needed) continue;
            }

            const GlmEntry &track_rep = var.glm_entries[tg.track_entry_idx];

            // Pre-read the super-window for this track (touches the cache once)
            const float *super_ptr = nullptr;
            int64_t super_sbin = 0;
            int64_t super_count = 0;
            unsigned bin_size = track_rep.cached_bin_size;

            if (track_rep.cached_fixedbin && bin_size > 0) {
                int64_t sw_start = center + tg.min_sshift;
                int64_t sw_end = center + tg.max_eshift;
                super_sbin = sw_start / (int64_t)bin_size;
                int64_t super_ebin = (int64_t)std::ceil(sw_end / (double)bin_size);
                super_ptr = track_rep.cached_fixedbin->get_mmap_bins_ptr(
                    super_sbin, super_ebin - super_sbin, super_count);
            }

            for (int sg_idx : tg.scaling_group_indices) {
                const GlmScalingGroup &group = var.glm_scaling_groups[sg_idx];
                const GlmEntry &rep = var.glm_entries[group.representative_idx];

                // Per-bin skip at scaling group level
                if (var.glm_num_bins > 1) {
                    bool sg_needed = false;
                    for (int eidx : group.entry_indices) {
                        if (bin_weights[var.glm_entries[eidx].weight_offset] != 0.0 ||
                            (has_interactions && m_inter_referenced[eidx])) {
                            sg_needed = true;
                            break;
                        }
                    }
                    if (!sg_needed) continue;
                }

                double raw;

                // Try to aggregate from the pre-read super-window buffer
                if (super_ptr && rep.cached_fixedbin && bin_size > 0) {
                    int64_t win_start = center + rep.sshift;
                    int64_t win_end = center + rep.eshift;
                    int64_t sbin = win_start / (int64_t)bin_size;
                    int64_t ebin = (int64_t)std::ceil(win_end / (double)bin_size);

                    // Compute offsets into the super-window buffer
                    int64_t buf_off = sbin - super_sbin;
                    int64_t buf_end = ebin - super_sbin;
                    if (buf_off < 0) buf_off = 0;
                    if (buf_end > super_count) buf_end = super_count;
                    int64_t num_bins = buf_end - buf_off;

                    if (num_bins > 0) {
                        const float *bins = super_ptr + buf_off;
                        if (rep.inner_is_lse) {
                            double max_val = -numeric_limits<double>::infinity();
                            for (int64_t i = 0; i < num_bins; i++) {
                                float v = bins[i];
                                if (!std::isnan(v) && !std::isinf(v)) {
                                    if ((double)v > max_val) max_val = (double)v;
                                }
                            }
                            if (!std::isfinite(max_val)) {
                                raw = numeric_limits<double>::quiet_NaN();
                            } else {
                                double exp_sum = 0.0;
                                for (int64_t i = 0; i < num_bins; i++) {
                                    float v = bins[i];
                                    if (!std::isnan(v) && !std::isinf(v)) {
                                        exp_sum += std::exp((double)v - max_val);
                                    }
                                }
                                raw = max_val + std::log(exp_sum);
                            }
                        } else {
                            double acc = 0.0;
                            bool has_val = false;
                            for (int64_t i = 0; i < num_bins; i++) {
                                float v = bins[i];
                                if (!std::isnan(v) && !std::isinf(v)) {
                                    acc += (double)v;
                                    has_val = true;
                                }
                            }
                            raw = has_val ? acc : numeric_limits<double>::quiet_NaN();
                        }
                    } else {
                        raw = numeric_limits<double>::quiet_NaN();
                    }
                } else {
                    // Fallback: sparse track or no mmap
                    int64_t win_start = center + rep.sshift;
                    int64_t win_end = center + rep.eshift;
                    raw = aggregate_window(rep, win_start, win_end, rep.inner_is_lse);
                }

                // All-NaN window → entire position is NaN
                if (std::isnan(raw)) {
                    var.var[idx] = NAN;
                    return;
                }
                if (!std::isfinite(raw)) raw = 0.0;
                double scaled = apply_scaling(raw, rep.scaling, var.glm_scale_factor);

                if (has_interactions) {
                    for (int idx : group.entry_indices) {
                        var.glm_scaled_cache[idx] = scaled;
                    }
                }

                for (int idx : group.entry_indices) {
                    const GlmEntry &entry = var.glm_entries[idx];
                    double transformed = scaled;
                    if (entry.transform.enabled) {
                        transformed = apply_transform(scaled, entry.transform);
                    }
                    if (!std::isfinite(transformed)) transformed = 0.0;
                    result += bin_weights[entry.weight_offset] * transformed;
                }
            }
        }
    } else if (!var.glm_scaling_groups.empty()) {
        // Scaling groups without track groups
        for (const GlmScalingGroup &group : var.glm_scaling_groups) {
            // Per-bin skip
            if (var.glm_num_bins > 1) {
                bool sg_needed = false;
                for (int eidx : group.entry_indices) {
                    if (bin_weights[var.glm_entries[eidx].weight_offset] != 0.0 ||
                        (has_interactions && m_inter_referenced[eidx])) {
                        sg_needed = true;
                        break;
                    }
                }
                if (!sg_needed) continue;
            }
            const GlmEntry &rep = var.glm_entries[group.representative_idx];
            int64_t win_start = center + rep.sshift;
            int64_t win_end = center + rep.eshift;
            double raw = aggregate_window(rep, win_start, win_end, rep.inner_is_lse);
            if (std::isnan(raw)) {
                var.var[idx] = NAN;
                return;
            }
            if (!std::isfinite(raw)) raw = 0.0;
            double scaled = apply_scaling(raw, rep.scaling, var.glm_scale_factor);
            if (has_interactions) {
                for (int idx : group.entry_indices) {
                    var.glm_scaled_cache[idx] = scaled;
                }
            }
            for (int idx : group.entry_indices) {
                const GlmEntry &entry = var.glm_entries[idx];
                double transformed = scaled;
                if (entry.transform.enabled) {
                    transformed = apply_transform(scaled, entry.transform);
                }
                if (!std::isfinite(transformed)) transformed = 0.0;
                result += bin_weights[entry.weight_offset] * transformed;
            }
        }
    } else {
        // Fallback: per-entry loop (kernel_bins case or no groups built)
        for (int i = 0; i < N; i++) {
            // Per-bin skip
            if (var.glm_num_bins > 1 && bin_weights[i] == 0.0 &&
                !(has_interactions && m_inter_referenced[i])) continue;

            const GlmEntry &entry = var.glm_entries[i];

            double raw;
            if (!has_kernel_bins) {
                int64_t win_start = center + entry.sshift;
                int64_t win_end = center + entry.eshift;
                raw = aggregate_window(entry, win_start, win_end, entry.inner_is_lse);
            } else {
                // Kernel sub-bin smoothing
                int64_t win_center = center + (entry.sshift + entry.eshift) / 2;
                int B = (int)var.glm_kernel_bins.size();

                if (entry.inner_is_lse) {
                    double acc = -numeric_limits<double>::infinity();
                    for (int kb = 0; kb < B; kb++) {
                        double bin_val = read_single_bin(entry, win_center + (int64_t)var.glm_kernel_bins[kb]);
                        if (std::isfinite(bin_val) && entry.kernel_weights[kb] > 0) {
                            lse_accumulate(acc, std::log(entry.kernel_weights[kb]) + bin_val);
                        }
                    }
                    raw = std::isfinite(acc) ? acc : numeric_limits<double>::quiet_NaN();
                } else {
                    raw = 0.0;
                    bool has_val = false;
                    for (int kb = 0; kb < B; kb++) {
                        double bin_val = read_single_bin(entry, win_center + (int64_t)var.glm_kernel_bins[kb]);
                        if (std::isfinite(bin_val)) {
                            raw += entry.kernel_weights[kb] * bin_val;
                            has_val = true;
                        }
                    }
                    if (!has_val) raw = numeric_limits<double>::quiet_NaN();
                }
            }

            if (std::isnan(raw)) {
                var.var[idx] = NAN;
                return;
            }
            if (!std::isfinite(raw)) raw = 0.0;
            double scaled = apply_scaling(raw, entry.scaling, var.glm_scale_factor);

            if (has_interactions) {
                var.glm_scaled_cache[i] = scaled;
            }

            double transformed = scaled;
            if (entry.transform.enabled) {
                transformed = apply_transform(scaled, entry.transform);
            }
            if (!std::isfinite(transformed)) transformed = 0.0;
            result += bin_weights[entry.weight_offset] * transformed;
        }
    }

    // ---- Interactions ----
    for (const auto &inter : var.glm_interactions) {
        double product;

        if (!has_kernel_bins) {
            // Product of cached scaled values, re-normalized by scale_factor
            product = var.glm_scaled_cache[inter.entry_i]
                    * var.glm_scaled_cache[inter.entry_j]
                    / var.glm_scale_factor;
        } else {
            // Per-bin product, kernel-weighted (design doc §3.3):
            // product = Σ_b kernel[b] × scaled_i_at_b × scaled_j_at_b / scale_factor
            // Uses entry_i's kernel weights as the shared kernel[b].
            const GlmEntry &ei = var.glm_entries[inter.entry_i];
            const GlmEntry &ej = var.glm_entries[inter.entry_j];
            int64_t win_center_i = center + (ei.sshift + ei.eshift) / 2;
            int64_t win_center_j = center + (ej.sshift + ej.eshift) / 2;
            int B = (int)var.glm_kernel_bins.size();

            product = 0.0;
            bool inter_has_val = false;
            for (int kb = 0; kb < B; kb++) {
                double vi = read_single_bin(ei, win_center_i + (int64_t)var.glm_kernel_bins[kb]);
                double vj = read_single_bin(ej, win_center_j + (int64_t)var.glm_kernel_bins[kb]);
                if (!std::isfinite(vi) || !std::isfinite(vj)) continue;
                double si = apply_scaling(vi, ei.scaling, var.glm_scale_factor);
                double sj = apply_scaling(vj, ej.scaling, var.glm_scale_factor);

                product += ei.kernel_weights[kb] * si * sj / var.glm_scale_factor;
                inter_has_val = true;
            }
            if (!inter_has_val) {
                var.var[idx] = NAN;
                return;
            }
        }

        // Apply interaction-specific transform
        if (inter.transform.enabled) {
            product = apply_transform(product, inter.transform);
        }

        // NaN/Inf → 0
        if (!std::isfinite(product)) product = 0.0;

        result += bin_inter_weights[inter.weight_offset] * product;
    }

    var.var[idx] = result;
}

double GlmVarProcessor::aggregate_window(
    const TrackExpressionVars::Track_var::GlmEntry &entry,
    int64_t start, int64_t end,
    bool is_lse)
{
    // Clamp window to valid coordinate range (mirrors iterator modifier clipping)
    if (start < 0) start = 0;
    if (start >= end) return 0.0;

    if (entry.cached_fixedbin) {
        unsigned bin_size = entry.cached_bin_size;
        if (bin_size == 0) return 0.0;

        int64_t sbin = start / bin_size;
        int64_t ebin = (int64_t)std::ceil(end / (double)bin_size);
        int64_t num_bins = ebin - sbin;
        if (num_bins <= 0) return 0.0;

        int64_t raw_count;
        const float *raw = entry.cached_fixedbin->get_mmap_bins_ptr(sbin, num_bins, raw_count);
        if (!raw) {
            raw_count = entry.cached_fixedbin->read_bins_bulk(sbin, num_bins, const_cast<std::vector<float>&>(m_raw_bins));
            raw = m_raw_bins.data();
        }

        if (raw_count <= 0) return 0.0;

        if (is_lse) {
            // Two-pass numerically stable LSE
            double max_val = -numeric_limits<double>::infinity();
            for (int64_t i = 0; i < raw_count; i++) {
                float v = raw[i];
                if (!std::isnan(v) && !std::isinf(v)) {
                    if ((double)v > max_val) max_val = (double)v;
                }
            }
            if (!std::isfinite(max_val)) return numeric_limits<double>::quiet_NaN();

            double exp_sum = 0.0;
            for (int64_t i = 0; i < raw_count; i++) {
                float v = raw[i];
                if (!std::isnan(v) && !std::isinf(v)) {
                    exp_sum += std::exp((double)v - max_val);
                }
            }
            return max_val + std::log(exp_sum);
        } else {
            // Sum
            double acc = 0.0;
            bool has_val = false;
            for (int64_t i = 0; i < raw_count; i++) {
                float v = raw[i];
                if (!std::isnan(v) && !std::isinf(v)) {
                    acc += (double)v;
                    has_val = true;
                }
            }
            return has_val ? acc : numeric_limits<double>::quiet_NaN();
        }
    }

    if (entry.cached_sparse) {
        const GIntervals &intervals = entry.cached_sparse->get_intervals();
        const vector<float> &vals = entry.cached_sparse->get_vals();

        if (is_lse) {
            double acc = -numeric_limits<double>::infinity();
            for (size_t i = 0; i < intervals.size(); i++) {
                if (intervals[i].start >= end) break;
                if (intervals[i].end <= start) continue;
                float v = vals[i];
                if (!std::isnan(v) && !std::isinf(v)) {
                    lse_accumulate(acc, (double)v);
                }
            }
            return std::isfinite(acc) ? acc : numeric_limits<double>::quiet_NaN();
        } else {
            double acc = 0.0;
            bool has_val = false;
            for (size_t i = 0; i < intervals.size(); i++) {
                if (intervals[i].start >= end) break;
                if (intervals[i].end <= start) continue;
                float v = vals[i];
                if (!std::isnan(v) && !std::isinf(v)) {
                    acc += (double)v;
                    has_val = true;
                }
            }
            return has_val ? acc : numeric_limits<double>::quiet_NaN();
        }
    }

    // Track not available for this chromosome
    return numeric_limits<double>::quiet_NaN();
}

double GlmVarProcessor::read_single_bin(
    const TrackExpressionVars::Track_var::GlmEntry &entry,
    int64_t pos)
{
    if (entry.cached_fixedbin) {
        unsigned bin_size = entry.cached_bin_size;
        if (bin_size == 0) return numeric_limits<double>::quiet_NaN();

        int64_t bin_idx = pos / bin_size;
        int64_t raw_count;
        const float *raw = entry.cached_fixedbin->get_mmap_bins_ptr(bin_idx, 1, raw_count);
        if (!raw || raw_count <= 0)
            return numeric_limits<double>::quiet_NaN();
        float v = raw[0];
        if (std::isnan(v) || std::isinf(v))
            return numeric_limits<double>::quiet_NaN();
        return (double)v;
    }

    if (entry.cached_sparse) {
        // For sparse tracks, read a single-bin-width window
        const GIntervals &intervals = entry.cached_sparse->get_intervals();
        const vector<float> &vals = entry.cached_sparse->get_vals();
        for (size_t i = 0; i < intervals.size(); i++) {
            if (intervals[i].start > pos) break;
            if (intervals[i].end <= pos) continue;
            float v = vals[i];
            if (!std::isnan(v) && !std::isinf(v))
                return (double)v;
        }
        return numeric_limits<double>::quiet_NaN();
    }

    return numeric_limits<double>::quiet_NaN();
}
