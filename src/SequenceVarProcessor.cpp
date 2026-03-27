#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <string>

#include "SequenceVarProcessor.h"
#include "KmerCounter.h"
#include "PWMScorer.h"
#include "PWMEditDistanceScorer.h"
#include "PWMLseEditDistanceScorer.h"
#include "MaskedBpCounter.h"
#include "rdbutils.h"

using namespace rdb;
using namespace std;

namespace {

// Aggregation modes for edit distance multi-part scoring
enum class EditDistAggMode {
	MIN_EDITS,     // PWM_EDIT_DISTANCE, PWM_EDIT_DISTANCE_LSE
	MIN_EDITS_POS, // PWM_EDIT_DISTANCE_POS, PWM_EDIT_DISTANCE_LSE_POS
	MAX_EDITS      // PWM_MAX_EDIT_DISTANCE
};

// Per-part scoring result for edit distance aggregation
struct EditDistPartResult {
	double score;        // return value of score_interval (edits or position depending on mode)
	double min_edits;    // from get_last_min_edits() (used by POS modes)
	double pwm_max_logp; // from get_last_pwm_max_logp() (used by MAX_EDITS mode)
	int64_t part_start;  // start coordinate of the scored part
};

// Map val_func to aggregation mode
EditDistAggMode get_agg_mode(TrackExpressionVars::Track_var::Val_func func) {
	switch (func) {
		case TrackExpressionVars::Track_var::PWM_EDIT_DISTANCE_POS:
		case TrackExpressionVars::Track_var::PWM_EDIT_DISTANCE_LSE_POS:
			return EditDistAggMode::MIN_EDITS_POS;
		case TrackExpressionVars::Track_var::PWM_MAX_EDIT_DISTANCE:
			return EditDistAggMode::MAX_EDITS;
		default:
			return EditDistAggMode::MIN_EDITS;
	}
}

// Aggregate edit distance results from multiple unmasked parts.
// seq_interval_start: start of the original (pre-filter) interval, for position offset correction.
double aggregate_edit_distance_parts(
	const vector<EditDistPartResult> &parts,
	EditDistAggMode mode,
	int64_t seq_interval_start)
{
	switch (mode) {
		case EditDistAggMode::MIN_EDITS: {
			double min_edits = numeric_limits<double>::quiet_NaN();
			for (const auto &r : parts) {
				if (std::isnan(min_edits) || (!std::isnan(r.score) && r.score < min_edits)) {
					min_edits = r.score;
				}
			}
			return min_edits;
		}
		case EditDistAggMode::MIN_EDITS_POS: {
			double best_edits = numeric_limits<double>::quiet_NaN();
			double best_pos = numeric_limits<double>::quiet_NaN();
			for (const auto &r : parts) {
				if (std::isnan(r.min_edits)) continue;
				if (std::isnan(best_edits) || r.min_edits < best_edits) {
					best_edits = r.min_edits;
					// Offset position: score is relative to the sub-interval,
					// adjust to be relative to the original seq_interval.
					// score is signed (positive=fwd, negative=rev for bidirect).
					double offset = static_cast<double>(r.part_start - seq_interval_start);
					best_pos = (r.score > 0) ? r.score + offset :
					           (r.score < 0) ? r.score - offset : r.score;
				}
			}
			return best_pos;
		}
		case EditDistAggMode::MAX_EDITS: {
			bool found = false;
			double best_logp = -numeric_limits<double>::infinity();
			double best_edits = numeric_limits<double>::quiet_NaN();
			for (const auto &r : parts) {
				if (!found || r.pwm_max_logp > best_logp) {
					best_logp = r.pwm_max_logp;
					best_edits = r.score;
					found = true;
				}
			}
			return found ? best_edits : numeric_limits<double>::quiet_NaN();
		}
	}
	return numeric_limits<double>::quiet_NaN(); // unreachable
}

} // anonymous namespace

void SequenceVarProcessor::classify_track_vars(
	TrackExpressionVars::Track_vars &track_vars)
{
	m_kmer_vtracks.clear();
	m_pwm_vtracks.clear();
	m_masked_vtracks.clear();
	m_pwm_edit_distance_vtracks.clear();
	m_pwm_lse_edit_distance_vtracks.clear();
	m_any_filtered = false;

	for (TrackExpressionVars::Track_vars::iterator ivar = track_vars.begin(); ivar != track_vars.end(); ++ivar) {
		if (TrackExpressionVars::is_pwm_lse_edit_distance_function(ivar->val_func)) {
			m_pwm_lse_edit_distance_vtracks.push_back(&*ivar);
		} else if (TrackExpressionVars::is_pwm_edit_distance_function(ivar->val_func)) {
			m_pwm_edit_distance_vtracks.push_back(&*ivar);
		} else if (TrackExpressionVars::is_pwm_function(ivar->val_func)) {
			m_pwm_vtracks.push_back(&*ivar);
		} else if (TrackExpressionVars::is_kmer_function(ivar->val_func)) {
			m_kmer_vtracks.push_back(&*ivar);
		} else if (TrackExpressionVars::is_masked_function(ivar->val_func)) {
			m_masked_vtracks.push_back(&*ivar);
		}
	}

	// Check if any sequence vtrack has a filter
	for (TrackExpressionVars::Track_var* ivar : m_pwm_vtracks) {
		if (ivar->filter) {
			m_any_filtered = true;
			break;
		}
	}
	if (!m_any_filtered) {
		for (TrackExpressionVars::Track_var* ivar : m_kmer_vtracks) {
			if (ivar->filter) {
				m_any_filtered = true;
				break;
			}
		}
	}
	if (!m_any_filtered) {
		for (TrackExpressionVars::Track_var* ivar : m_masked_vtracks) {
			if (ivar->filter) {
				m_any_filtered = true;
				break;
			}
		}
	}

	m_classification_done = true;
}

void SequenceVarProcessor::process_sequence_vars(
	TrackExpressionVars::Track_vars &track_vars,
	const GInterval &interval,
	unsigned idx)
{
	// Classify track variables once; reuse on subsequent calls
	if (!m_classification_done) {
		classify_track_vars(track_vars);
	}

	// Batch process if we have multiple sequence vtracks (threshold: 4+) AND no filters
	// (batch processing with filters is complex due to variable-length segments)
	if (!m_any_filtered && m_kmer_vtracks.size() + m_pwm_vtracks.size() + m_masked_vtracks.size() >= 4) {
		batch_process_sequence_vtracks(m_kmer_vtracks, m_pwm_vtracks, interval, idx);
		// Process masked vtracks individually even in batch mode (simpler than batching)
		for (TrackExpressionVars::Track_var* ivar : m_masked_vtracks) {
			const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;
			ivar->var[idx] = ivar->masked_counter->score_interval(seq_interval, m_iu.get_chromkey());
		}
		// Process edit distance vtracks individually (no batching support)
		for (TrackExpressionVars::Track_var* ivar : m_pwm_edit_distance_vtracks) {
			const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

			if (!ivar->pwm_edit_distance_scorer) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				continue;
			}

			if (ivar->filter) {
				vector<GInterval> unmasked_parts;
				ivar->filter->subtract(seq_interval, unmasked_parts);

				if (unmasked_parts.empty()) {
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
					continue;
				}

				if (unmasked_parts.size() == 1) {
					ivar->var[idx] = ivar->pwm_edit_distance_scorer->score_interval(unmasked_parts[0], m_iu.get_chromkey());
					continue;
				}

				// Multiple unmasked parts - score each and aggregate
				EditDistAggMode mode = get_agg_mode(ivar->val_func);
				vector<EditDistPartResult> part_results;
				part_results.reserve(unmasked_parts.size());
				for (const auto& part : unmasked_parts) {
					double score = ivar->pwm_edit_distance_scorer->score_interval(part, m_iu.get_chromkey());
					part_results.push_back({
						score,
						ivar->pwm_edit_distance_scorer->get_last_min_edits(),
						ivar->pwm_edit_distance_scorer->get_last_pwm_max_logp(),
						part.start
					});
				}
				ivar->var[idx] = aggregate_edit_distance_parts(part_results, mode, seq_interval.start);
			} else {
				ivar->var[idx] = ivar->pwm_edit_distance_scorer->score_interval(seq_interval, m_iu.get_chromkey());
			}
		}
		// Process LSE edit distance vtracks individually (no batching support)
		for (TrackExpressionVars::Track_var* ivar : m_pwm_lse_edit_distance_vtracks) {
			const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

			if (!ivar->pwm_lse_edit_distance_scorer) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				continue;
			}

			if (ivar->filter) {
				vector<GInterval> unmasked_parts;
				ivar->filter->subtract(seq_interval, unmasked_parts);

				if (unmasked_parts.empty()) {
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
					continue;
				}

				if (unmasked_parts.size() == 1) {
					ivar->var[idx] = ivar->pwm_lse_edit_distance_scorer->score_interval(unmasked_parts[0], m_iu.get_chromkey());
					continue;
				}

				// Multiple unmasked parts - score each and aggregate
				EditDistAggMode mode = get_agg_mode(ivar->val_func);
				vector<EditDistPartResult> part_results;
				part_results.reserve(unmasked_parts.size());
				for (const auto& part : unmasked_parts) {
					double score = ivar->pwm_lse_edit_distance_scorer->score_interval(part, m_iu.get_chromkey());
					part_results.push_back({
						score,
						ivar->pwm_lse_edit_distance_scorer->get_last_min_edits(),
						0.0, // LSE does not use pwm_max_logp
						part.start
					});
				}
				ivar->var[idx] = aggregate_edit_distance_parts(part_results, mode, seq_interval.start);
			} else {
				ivar->var[idx] = ivar->pwm_lse_edit_distance_scorer->score_interval(seq_interval, m_iu.get_chromkey());
			}
		}
	} else {
		// Process individually if too few or has filters
		process_individual_sequence_vars(m_kmer_vtracks, m_pwm_vtracks, m_masked_vtracks, m_pwm_edit_distance_vtracks, m_pwm_lse_edit_distance_vtracks, interval, idx);
	}
}

void SequenceVarProcessor::batch_process_sequence_vtracks(
	vector<TrackExpressionVars::Track_var*> &kmer_vtracks,
	vector<TrackExpressionVars::Track_var*> &pwm_vtracks,
	const GInterval &interval,
	unsigned idx)
{
	// Group kmers by: (interval, extend_params, strand)
	// Key: "chromid:start-end:max_extension:strand"
	unordered_map<string, vector<TrackExpressionVars::Track_var*>> kmer_groups;

	for (TrackExpressionVars::Track_var* var : kmer_vtracks) {
		if (!var->kmer_counter) continue;

		const GInterval &base_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
		char strand = var->kmer_counter->get_strand();
		int64_t extension = var->kmer_counter->get_extend() ? (int64_t)(var->kmer_counter->get_kmer().length() - 1) : 0;

		string key = to_string(base_interval.chromid) + ":" +
		             to_string(base_interval.start) + "-" +
		             to_string(base_interval.end) + ":" +
		             to_string(extension) + ":" +
		             to_string((int)strand);
		kmer_groups[key].push_back(var);
	}

	// Process each kmer group with HashMap optimization
	for (const auto& [key, group] : kmer_groups) {
		if (group.empty()) continue;

		TrackExpressionVars::Track_var* first = group[0];
		const GInterval &base_interval = first->seq_imdf1d ? first->seq_imdf1d->interval : interval;
		char strand = first->kmer_counter->get_strand();
		int64_t extension = first->kmer_counter->get_extend() ? (int64_t)(first->kmer_counter->get_kmer().length() - 1) : 0;

		// Calculate fetch interval with extension
		GInterval fetch_interval = base_interval;
		if (extension > 0) {
			fetch_interval.start = max((int64_t)0, base_interval.start - extension);
			fetch_interval.end = min((int64_t)m_iu.get_chromkey().get_chrom_size(base_interval.chromid),
			                         base_interval.end + extension);
		}

		// Build HashMap for forward and/or reverse strands
		unordered_map<string, vector<size_t>> fwd_positions;
		unordered_map<string, vector<size_t>> rev_positions;

		// Get max kmer length for this group
		size_t max_kmer_len = 0;
		for (TrackExpressionVars::Track_var* var : group) {
			max_kmer_len = max(max_kmer_len, var->kmer_counter->get_kmer().length());
		}

		// Fetch and process forward strand if needed
		if (strand == 0 || strand == 1) {
			GInterval fwd_interval = fetch_interval;
			fwd_interval.strand = 1;
			vector<char> seq;
			m_shared_seqfetch.read_interval(fwd_interval, m_iu.get_chromkey(), seq);
			string sequence(seq.begin(), seq.end());
			transform(sequence.begin(), sequence.end(), sequence.begin(),
			          [](unsigned char c) { return toupper(c); });

			// Build HashMap: scan once, record all kmer positions
			for (size_t pos = 0; pos + max_kmer_len <= sequence.length(); pos++) {
				for (size_t k = 1; k <= max_kmer_len && pos + k <= sequence.length(); k++) {
					string kmer = sequence.substr(pos, k);
					fwd_positions[kmer].push_back(pos);
				}
			}
		}

		// Fetch and process reverse strand if needed
		if (strand == 0 || strand == -1) {
			GInterval rev_interval = fetch_interval;
			rev_interval.strand = -1;
			vector<char> seq;
			m_shared_seqfetch.read_interval(rev_interval, m_iu.get_chromkey(), seq);
			string sequence(seq.begin(), seq.end());
			transform(sequence.begin(), sequence.end(), sequence.begin(),
			          [](unsigned char c) { return toupper(c); });

			// Build HashMap for reverse strand
			for (size_t pos = 0; pos + max_kmer_len <= sequence.length(); pos++) {
				for (size_t k = 1; k <= max_kmer_len && pos + k <= sequence.length(); k++) {
					string kmer = sequence.substr(pos, k);
					rev_positions[kmer].push_back(pos);
				}
			}
		}

		// Distribute to all kmer counters in this group
		for (TrackExpressionVars::Track_var* var : group) {
			const GInterval &var_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
			const string &target_kmer = var->kmer_counter->get_kmer();

			var->var[idx] = var->kmer_counter->score_from_positions(
				fwd_positions[target_kmer],
				rev_positions[target_kmer],
				var_interval,
				fetch_interval);
		}
	}

	// Process PWM vtracks using score_interval
	// The shared cache prevents redundant disk I/O
	for (TrackExpressionVars::Track_var* var : pwm_vtracks) {
		const GInterval &seq_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
		var->var[idx] = var->pwm_scorer->score_interval(seq_interval, m_iu.get_chromkey());
	}
}

void SequenceVarProcessor::process_individual_sequence_vars(
	vector<TrackExpressionVars::Track_var*> &kmer_vtracks,
	vector<TrackExpressionVars::Track_var*> &pwm_vtracks,
	vector<TrackExpressionVars::Track_var*> &masked_vtracks,
	vector<TrackExpressionVars::Track_var*> &pwm_edit_distance_vtracks,
	vector<TrackExpressionVars::Track_var*> &pwm_lse_edit_distance_vtracks,
	const GInterval &interval,
	unsigned idx)
{
	// Process PWM vtracks
	for (TrackExpressionVars::Track_var* ivar : pwm_vtracks) {
		const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

		// Check if filter applies
		if (ivar->filter) {
			vector<GInterval> unmasked_parts;
			ivar->filter->subtract(seq_interval, unmasked_parts);

			if (unmasked_parts.empty()) {
				// Completely masked
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			} else if (unmasked_parts.size() == 1) {
				// Single unmasked part
				ivar->var[idx] = ivar->pwm_scorer->score_interval(unmasked_parts[0], m_iu.get_chromkey());
			} else {
				// Multiple unmasked parts - score each and aggregate
				// For PWM: sum scores (additive in log space), max, or count
				double result = 0.0;
				bool first = true;

				for (const auto& part : unmasked_parts) {
					double part_score = ivar->pwm_scorer->score_interval(part, m_iu.get_chromkey());

					if (ivar->val_func == TrackExpressionVars::Track_var::PWM_MAX) {
						if (first || part_score > result) {
							result = part_score;
						}
					} else if (ivar->val_func == TrackExpressionVars::Track_var::PWM_MAX_POS) {
						// For position, take the position with max score
						if (first || part_score > result) {
							result = part_score;
						}
					} else if (ivar->val_func == TrackExpressionVars::Track_var::PWM_COUNT) {
						// Sum counts
						result += part_score;
					} else {
						// PWM (total likelihood) - sum in log space
						result += part_score;
					}
					first = false;
				}

				ivar->var[idx] = result;
			}
		} else {
			// No filter
			ivar->var[idx] = ivar->pwm_scorer->score_interval(seq_interval, m_iu.get_chromkey());
		}
	}

	// Process kmer vtracks
	for (TrackExpressionVars::Track_var* ivar : kmer_vtracks) {
		const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

		if (ivar->kmer_counter) {
			// Check if filter applies
			if (ivar->filter) {
				vector<GInterval> unmasked_parts;
				ivar->filter->subtract(seq_interval, unmasked_parts);

				if (unmasked_parts.empty()) {
					// Completely masked
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				} else if (unmasked_parts.size() == 1) {
					// Single unmasked part
					ivar->var[idx] = ivar->kmer_counter->score_interval(unmasked_parts[0], m_iu.get_chromkey());
				} else {
					// Multiple unmasked parts - score each and aggregate
					double total_count = 0.0;
					int64_t total_bases = 0;

					for (const auto& part : unmasked_parts) {
						double part_score = ivar->kmer_counter->score_interval(part, m_iu.get_chromkey());

						if (ivar->val_func == TrackExpressionVars::Track_var::KMER_COUNT) {
							// Sum counts
							total_count += part_score;
						} else {
							// KMER_FRAC - need to weight by length
							total_count += part_score * part.range();
							total_bases += part.range();
						}
					}

					if (ivar->val_func == TrackExpressionVars::Track_var::KMER_FRAC && total_bases > 0) {
						ivar->var[idx] = total_count / total_bases;
					} else {
						ivar->var[idx] = total_count;
					}
				}
			} else {
				// No filter
				ivar->var[idx] = ivar->kmer_counter->score_interval(seq_interval, m_iu.get_chromkey());
			}
		} else {
			ivar->var[idx] = numeric_limits<double>::quiet_NaN();
		}
	}

	// Process masked vtracks
	for (TrackExpressionVars::Track_var* ivar : masked_vtracks) {
		const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

		if (ivar->masked_counter) {
			if (ivar->filter) {
				// Filter application: score unmasked parts and aggregate
				vector<GInterval> unmasked_parts;
				ivar->filter->subtract(seq_interval, unmasked_parts);

				if (unmasked_parts.empty()) {
					// Completely masked
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				} else if (unmasked_parts.size() == 1) {
					// Single unmasked part
					ivar->var[idx] = ivar->masked_counter->score_interval(unmasked_parts[0], m_iu.get_chromkey());
				} else {
					// Multiple unmasked parts - score each and aggregate
					double total_masked_count = 0.0;
					int64_t total_bases = 0;

					for (const auto& part : unmasked_parts) {
						double part_score = ivar->masked_counter->score_interval(part, m_iu.get_chromkey());

						if (ivar->val_func == TrackExpressionVars::Track_var::MASKED_COUNT) {
							// Sum counts
							total_masked_count += part_score;
						} else {
							// MASKED_FRAC - need to weight by length
							total_masked_count += part_score * part.range();
							total_bases += part.range();
						}
					}

					if (ivar->val_func == TrackExpressionVars::Track_var::MASKED_FRAC && total_bases > 0) {
						ivar->var[idx] = total_masked_count / total_bases;
					} else {
						ivar->var[idx] = total_masked_count;
					}
				}
			} else {
				// No filter
				ivar->var[idx] = ivar->masked_counter->score_interval(seq_interval, m_iu.get_chromkey());
			}
		} else {
			ivar->var[idx] = numeric_limits<double>::quiet_NaN();
		}
	}

	// Process PWM edit distance vtracks
	for (TrackExpressionVars::Track_var* ivar : pwm_edit_distance_vtracks) {
		const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

		if (!ivar->pwm_edit_distance_scorer) {
			ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			continue;
		}

		auto evaluate_interval = [&](const GInterval& interval_to_score) -> double {
			return ivar->pwm_edit_distance_scorer->score_interval(interval_to_score, m_iu.get_chromkey());
		};

		if (ivar->filter) {
			vector<GInterval> unmasked_parts;
			ivar->filter->subtract(seq_interval, unmasked_parts);

			if (unmasked_parts.empty()) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				continue;
			}

			if (unmasked_parts.size() == 1) {
				ivar->var[idx] = evaluate_interval(unmasked_parts[0]);
				continue;
			}

			// Multiple unmasked parts - score each and aggregate
			EditDistAggMode mode = get_agg_mode(ivar->val_func);
			vector<EditDistPartResult> part_results;
			part_results.reserve(unmasked_parts.size());
			for (const auto& part : unmasked_parts) {
				double score = evaluate_interval(part);
				part_results.push_back({
					score,
					ivar->pwm_edit_distance_scorer->get_last_min_edits(),
					ivar->pwm_edit_distance_scorer->get_last_pwm_max_logp(),
					part.start
				});
			}
			ivar->var[idx] = aggregate_edit_distance_parts(part_results, mode, seq_interval.start);
		} else {
			ivar->var[idx] = evaluate_interval(seq_interval);
		}
	}

	// Process LSE edit distance vtracks
	for (TrackExpressionVars::Track_var* ivar : pwm_lse_edit_distance_vtracks) {
		const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : interval;

		if (!ivar->pwm_lse_edit_distance_scorer) {
			ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			continue;
		}

		auto evaluate_interval = [&](const GInterval& interval_to_score) -> double {
			return ivar->pwm_lse_edit_distance_scorer->score_interval(interval_to_score, m_iu.get_chromkey());
		};

		if (ivar->filter) {
			vector<GInterval> unmasked_parts;
			ivar->filter->subtract(seq_interval, unmasked_parts);

			if (unmasked_parts.empty()) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				continue;
			}

			if (unmasked_parts.size() == 1) {
				ivar->var[idx] = evaluate_interval(unmasked_parts[0]);
				continue;
			}

			// Multiple unmasked parts - score each and aggregate
			EditDistAggMode mode = get_agg_mode(ivar->val_func);
			vector<EditDistPartResult> part_results;
			part_results.reserve(unmasked_parts.size());
			for (const auto& part : unmasked_parts) {
				double score = evaluate_interval(part);
				part_results.push_back({
					score,
					ivar->pwm_lse_edit_distance_scorer->get_last_min_edits(),
					0.0, // LSE does not use pwm_max_logp
					part.start
				});
			}
			ivar->var[idx] = aggregate_edit_distance_parts(part_results, mode, seq_interval.start);
		} else {
			ivar->var[idx] = evaluate_interval(seq_interval);
		}
	}
}

