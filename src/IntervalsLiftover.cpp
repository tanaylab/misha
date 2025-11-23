#include "port.h"

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "AggregationHelpers.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gintervs_liftover(SEXP _src_intervs, SEXP _chain, SEXP _src_overlap_policy, SEXP _tgt_overlap_policy, SEXP _canonic, SEXP _include_metadata, SEXP _value_col, SEXP _multi_target_agg, SEXP _nth_param, SEXP _na_rm, SEXP _min_n, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_src_overlap_policy) || Rf_length(_src_overlap_policy) != 1)
			verror("Source overlap policy argument is not a string");

		if (!Rf_isString(_tgt_overlap_policy) || Rf_length(_tgt_overlap_policy) != 1)
			verror("Target overlap policy argument is not a string");

		bool canonic = Rf_asLogical(_canonic);
		bool include_metadata = Rf_asLogical(_include_metadata);

		// Parse aggregation parameters
		bool use_aggregation = false;
		const char *value_col_name = NULL;
		AggregationConfig agg_cfg;
		agg_cfg.na_rm = true;
		agg_cfg.min_n = -1;
		agg_cfg.nth_index = -1;
		SEXP value_col_data = R_NilValue;
		int value_col_idx = -1;

		if (!Rf_isNull(_value_col)) {
			use_aggregation = true;
			if (!Rf_isString(_value_col) || Rf_length(_value_col) != 1)
				verror("value_col must be a single string");

			value_col_name = CHAR(STRING_ELT(_value_col, 0));

			// Parse multi_target_agg
			if (!Rf_isString(_multi_target_agg) || Rf_length(_multi_target_agg) != 1)
				verror("multi_target_agg must be a single string");
			const char *agg_type_str = CHAR(STRING_ELT(_multi_target_agg, 0));
			agg_cfg.type = parse_aggregation_type(agg_type_str);

			// Parse na.rm
			if (!Rf_isLogical(_na_rm) || Rf_length(_na_rm) != 1)
				verror("na.rm must be a single logical value");
			agg_cfg.na_rm = Rf_asLogical(_na_rm);

			// Parse min_n
			if (!Rf_isInteger(_min_n) || Rf_length(_min_n) != 1)
				verror("min_n must be a single integer");
			int min_n_val = INTEGER(_min_n)[0];
			if (min_n_val != NA_INTEGER)
				agg_cfg.min_n = min_n_val;

			// Parse nth_index for nth aggregation
			if (!Rf_isInteger(_nth_param) || Rf_length(_nth_param) != 1)
				verror("nth_param must be a single integer");
			int nth_val = INTEGER(_nth_param)[0];
			if (nth_val != NA_INTEGER)
				agg_cfg.nth_index = nth_val;
		}

		IntervUtils iu(_envir);
		const char *src_overlap_policy = CHAR(STRING_ELT(_src_overlap_policy, 0));
		const char *tgt_overlap_policy = CHAR(STRING_ELT(_tgt_overlap_policy, 0));
		std::string effective_tgt_policy = tgt_overlap_policy;
		if (!strcmp(tgt_overlap_policy, "auto"))
			effective_tgt_policy = "auto_score";

		ChainIntervals chain_intervs;
		vector<string> src_id2chrom;

		iu.convert_rchain_intervs(_chain, chain_intervs, src_id2chrom);

		// Handle target overlaps first
		chain_intervs.sort_by_tgt();
		chain_intervs.handle_tgt_overlaps(effective_tgt_policy, iu.get_chromkey(), src_id2chrom);
		chain_intervs.set_tgt_overlap_policy(effective_tgt_policy);

		// Handle source overlaps
		chain_intervs.sort_by_src();
		chain_intervs.handle_src_overlaps(src_overlap_policy, iu.get_chromkey(), src_id2chrom);

		// Build auxiliary structures for efficient source interval mapping
		chain_intervs.buildSrcAux();

		GenomeChromKey src_chromkey;
		for (vector<string>::const_iterator ichrom = src_id2chrom.begin(); ichrom != src_id2chrom.end(); ++ichrom)
			src_chromkey.add_chrom(*ichrom, numeric_limits<int64_t>::max());

		GIntervals src_intervs1d;
		GIntervals2D src_intervs2d;
		iu.convert_rintervs(_src_intervs, &src_intervs1d, &src_intervs2d, false, &src_chromkey);
		src_intervs1d.sort();
		src_intervs2d.sort();

		// Extract value column if aggregation is enabled
		vector<double> src_values;
		if (use_aggregation) {
			// Find the value column in the source intervals dataframe
			SEXP col_names = Rf_getAttrib(_src_intervs, R_NamesSymbol);
			if (Rf_isNull(col_names))
				verror("Source intervals must have column names when value_col is specified");

			int ncols = Rf_length(_src_intervs);
			value_col_idx = -1;
			for (int i = 0; i < ncols; ++i) {
				if (!strcmp(CHAR(STRING_ELT(col_names, i)), value_col_name)) {
					value_col_idx = i;
					value_col_data = VECTOR_ELT(_src_intervs, i);
					break;
				}
			}

			if (value_col_idx < 0)
				verror("value_col '%s' not found in source intervals", value_col_name);

			// Extract values as doubles
			int n_rows = Rf_length(value_col_data);
			src_values.resize(n_rows);

			if (Rf_isReal(value_col_data)) {
				double *vals = REAL(value_col_data);
				for (int i = 0; i < n_rows; ++i)
					src_values[i] = vals[i];
			} else if (Rf_isInteger(value_col_data)) {
				int *vals = INTEGER(value_col_data);
				for (int i = 0; i < n_rows; ++i)
					src_values[i] = (vals[i] == NA_INTEGER) ? numeric_limits<double>::quiet_NaN() : static_cast<double>(vals[i]);
			} else {
				verror("value_col must be numeric or integer");
			}
		}

		vector<int> src_indices;
		vector<int64_t> chain_ids;
		vector<double> scores;
		vector<double> values;  // Track values for aggregation
		GIntervals tgt_intervs1d;
		GIntervals2D tgt_intervs2d;

		// 1D intervals
		if (!src_intervs1d.empty()) {
			GIntervals tmp_tgt_intervs;
			vector<ChainMappingMetadata> tmp_metadata;
			ChainIntervals::const_iterator hint = chain_intervs.begin();

			for (GIntervals::const_iterator isrc_interv = src_intervs1d.begin(); isrc_interv != src_intervs1d.end(); ++isrc_interv) {
				tmp_metadata.clear();
				hint = chain_intervs.map_interval(*isrc_interv, tmp_tgt_intervs, hint, &tmp_metadata);
				if (!tmp_tgt_intervs.empty()) {
					int src_idx = iu.get_orig_interv_idx(*isrc_interv);
					tgt_intervs1d.insert(tgt_intervs1d.end(), tmp_tgt_intervs.begin(), tmp_tgt_intervs.end());
					src_indices.insert(src_indices.end(), tmp_tgt_intervs.size(), src_idx + 1);
					for (size_t i = 0; i < tmp_metadata.size(); ++i) {
						chain_ids.push_back(tmp_metadata[i].chain_id);
						scores.push_back(tmp_metadata[i].score);
					}
					if (use_aggregation) {
						// Track value for this source interval
						double val = src_values[src_idx];
						values.insert(values.end(), tmp_tgt_intervs.size(), val);
					}
					iu.verify_max_data_size(tgt_intervs1d.size(), "Result");
				}
			}
		}
		// 2D intervals
		else {
			GInterval src_intervs[2];
			GIntervals tmp_tgt_intervs[2];
			ChainIntervals::const_iterator hints[2] = { chain_intervs.begin(), chain_intervs.begin() };

			for (GIntervals2D::const_iterator isrc_interv = src_intervs2d.begin(); isrc_interv != src_intervs2d.end(); ++isrc_interv) {
				src_intervs[0].chromid = isrc_interv->chromid1();
				src_intervs[0].start = isrc_interv->start1();
				src_intervs[0].end = isrc_interv->end1();
				src_intervs[1].chromid = isrc_interv->chromid2();
				src_intervs[1].start = isrc_interv->start2();
				src_intervs[1].end = isrc_interv->end2();

				hints[0] = chain_intervs.map_interval(src_intervs[0], tmp_tgt_intervs[0], hints[0]);
				hints[1] = chain_intervs.map_interval(src_intervs[1], tmp_tgt_intervs[1], hints[1]);

				for (GIntervals::const_iterator iinterv1 = tmp_tgt_intervs[0].begin(); iinterv1 != tmp_tgt_intervs[0].end(); ++iinterv1) {
					for (GIntervals::const_iterator iinterv2 = tmp_tgt_intervs[1].begin(); iinterv2 != tmp_tgt_intervs[1].end(); ++iinterv2) {
						tgt_intervs2d.push_back(GInterval2D(iinterv1->chromid, iinterv1->start, iinterv1->end, iinterv2->chromid, iinterv2->start, iinterv2->end));
						src_indices.push_back(iu.get_orig_interv_idx(*isrc_interv) + 1);
						iu.verify_max_data_size(tgt_intervs2d.size(), "Result");
					}
				}
			}
		}

		// Canonicalization: merge adjacent intervals within the same intervalID and chain_id
		if (canonic && !tgt_intervs1d.empty()) {
			// Create index pairs for sorting: (index, intervalID, chain_id, chromid, start)
			vector<size_t> order(tgt_intervs1d.size());
			for (size_t i = 0; i < order.size(); ++i)
				order[i] = i;

			// Sort by intervalID, then chain_id, then chromid, then start
			sort(order.begin(), order.end(), [&](size_t a, size_t b) {
				if (src_indices[a] != src_indices[b])
					return src_indices[a] < src_indices[b];
				if (chain_ids[a] != chain_ids[b])
					return chain_ids[a] < chain_ids[b];
				if (tgt_intervs1d[a].chromid != tgt_intervs1d[b].chromid)
					return tgt_intervs1d[a].chromid < tgt_intervs1d[b].chromid;
				return tgt_intervs1d[a].start < tgt_intervs1d[b].start;
			});

			// Merge adjacent intervals within same intervalID and chain_id
			GIntervals merged_intervs;
			vector<int> merged_indices;
			vector<int64_t> merged_chain_ids;
			vector<double> merged_scores;
			vector<double> merged_values;

			for (size_t i = 0; i < order.size(); ++i) {
				size_t idx = order[i];
				const GInterval &interval = tgt_intervs1d[idx];
				int intervalID = src_indices[idx];
				int64_t chainID = chain_ids[idx];
				double score = scores[idx];
				double val = use_aggregation ? values[idx] : 0.0;

				if (merged_intervs.empty() ||
				    merged_indices.back() != intervalID ||
				    merged_chain_ids.back() != chainID ||
				    merged_intervs.back().chromid != interval.chromid ||
				    merged_intervs.back().end != interval.start) {
					// Start new interval
					merged_intervs.push_back(interval);
					merged_indices.push_back(intervalID);
					merged_chain_ids.push_back(chainID);
					merged_scores.push_back(score);
					if (use_aggregation)
						merged_values.push_back(val);
				} else {
					// Merge with previous (extend end)
					// Values stay the same (same source interval)
					merged_intervs.back().end = interval.end;
				}
			}

			tgt_intervs1d = merged_intervs;
			src_indices = merged_indices;
			chain_ids = merged_chain_ids;
			scores = merged_scores;
			if (use_aggregation)
				values = merged_values;
		}

		// assemble the answer
		SEXP answer;
		unsigned num_interv_cols;

		if (!tgt_intervs1d.empty()) {
			// Add extra columns for intervalID, chain_id, score (if include_metadata), and value (if use_aggregation)
			int extra_cols = 2;  // intervalID, chain_id
			if (include_metadata) extra_cols++;  // score
			if (use_aggregation) extra_cols++;  // value
			answer = iu.convert_intervs(&tgt_intervs1d, GInterval::NUM_COLS + extra_cols);
			num_interv_cols = GInterval::NUM_COLS;
		} else if (!tgt_intervs2d.empty()) {
			answer = iu.convert_intervs(&tgt_intervs2d, GInterval2D::NUM_COLS + 1);
			num_interv_cols = GInterval2D::NUM_COLS;
		} else
			return R_NilValue;

		// Add intervalID, chain_id, and optionally score and value columns for 1D intervals
		if (!tgt_intervs1d.empty()) {
			SEXP rsrc_indices;
			SEXP rchain_ids;
			SEXP rscores = R_NilValue;
			SEXP rvalues = R_NilValue;
			SEXP col_names = Rf_getAttrib(answer, R_NamesSymbol);
			rprotect(col_names);

			rprotect(rsrc_indices = RSaneAllocVector(INTSXP, src_indices.size()));
			rprotect(rchain_ids = RSaneAllocVector(REALSXP, chain_ids.size()));

			for (size_t i = 0; i < src_indices.size(); ++i)
				INTEGER(rsrc_indices)[i] = src_indices[i];

			for (size_t i = 0; i < chain_ids.size(); ++i)
				REAL(rchain_ids)[i] = (double)chain_ids[i];

			SET_STRING_ELT(col_names, num_interv_cols, Rf_mkChar("intervalID"));
			SET_VECTOR_ELT(answer, num_interv_cols, rsrc_indices);
			SET_STRING_ELT(col_names, num_interv_cols + 1, Rf_mkChar("chain_id"));
			SET_VECTOR_ELT(answer, num_interv_cols + 1, rchain_ids);

			int col_idx = num_interv_cols + 2;
			int num_protected = 2;  // rsrc_indices, rchain_ids

			if (include_metadata) {
				rprotect(rscores = RSaneAllocVector(REALSXP, scores.size()));
				for (size_t i = 0; i < scores.size(); ++i)
					REAL(rscores)[i] = scores[i];
				SET_STRING_ELT(col_names, col_idx, Rf_mkChar("score"));
				SET_VECTOR_ELT(answer, col_idx, rscores);
				col_idx++;
				num_protected++;
			}

			if (use_aggregation) {
				rprotect(rvalues = RSaneAllocVector(REALSXP, values.size()));
				for (size_t i = 0; i < values.size(); ++i)
					REAL(rvalues)[i] = values[i];
				SET_STRING_ELT(col_names, col_idx, Rf_mkChar(value_col_name));
				SET_VECTOR_ELT(answer, col_idx, rvalues);
				num_protected++;
			}

			runprotect(num_protected);
			return answer;
		}

		// Add intervalID column for 2D intervals
        SEXP rsrc_indices;
        SEXP col_names = Rf_getAttrib(answer, R_NamesSymbol);
        rprotect(col_names);

        rprotect(rsrc_indices = RSaneAllocVector(INTSXP, src_indices.size()));

		for (vector<int>::const_iterator iindex = src_indices.begin(); iindex != src_indices.end(); ++iindex)
			INTEGER(rsrc_indices)[iindex - src_indices.begin()] = *iindex;

        SET_STRING_ELT(col_names, num_interv_cols, Rf_mkChar("intervalID"));
        SET_VECTOR_ELT(answer, num_interv_cols, rsrc_indices);

        runprotect(1); // col_names
        return answer;

	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
