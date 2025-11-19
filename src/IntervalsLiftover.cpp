#include "port.h"

#include <algorithm>
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gintervs_liftover(SEXP _src_intervs, SEXP _chain, SEXP _src_overlap_policy, SEXP _tgt_overlap_policy, SEXP _canonic, SEXP _include_metadata, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_src_overlap_policy) || Rf_length(_src_overlap_policy) != 1)
			verror("Source overlap policy argument is not a string");

		if (!Rf_isString(_tgt_overlap_policy) || Rf_length(_tgt_overlap_policy) != 1)
			verror("Target overlap policy argument is not a string");

		bool canonic = Rf_asLogical(_canonic);
		bool include_metadata = Rf_asLogical(_include_metadata);

		IntervUtils iu(_envir);
		const char *src_overlap_policy = CHAR(STRING_ELT(_src_overlap_policy, 0));
		const char *tgt_overlap_policy = CHAR(STRING_ELT(_tgt_overlap_policy, 0));
		std::string effective_tgt_policy = tgt_overlap_policy;

		// Strategy enum for clustering policies
		enum ClusterStrategy { STRAT_UNION, STRAT_SUM, STRAT_MAX };
		bool use_clustering = false;
		ClusterStrategy strat = STRAT_UNION;

		if (!strcmp(tgt_overlap_policy, "best_cluster_union") || !strcmp(tgt_overlap_policy, "best_source_cluster")) {
			effective_tgt_policy = "keep"; // Load ALL chains to resolve later
			use_clustering = true;
			strat = STRAT_UNION;
		} else if (!strcmp(tgt_overlap_policy, "best_cluster_sum")) {
			effective_tgt_policy = "keep";
			use_clustering = true;
			strat = STRAT_SUM;
		} else if (!strcmp(tgt_overlap_policy, "best_cluster_max")) {
			effective_tgt_policy = "keep";
			use_clustering = true;
			strat = STRAT_MAX;
		} else if (!strcmp(tgt_overlap_policy, "auto")) {
			effective_tgt_policy = "auto_score";
		}

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

		vector<int> src_indices;
		vector<int64_t> chain_ids;
		vector<double> scores;
		GIntervals tgt_intervs1d;
		GIntervals2D tgt_intervs2d;

		// 1D intervals
		if (!src_intervs1d.empty()) {
			GIntervals tmp_tgt_intervs;
			vector<ChainMappingMetadata> tmp_metadata;
			ChainIntervals::const_iterator hint = chain_intervs.begin();

			for (GIntervals::const_iterator isrc_interv = src_intervs1d.begin(); isrc_interv != src_intervs1d.end(); ++isrc_interv) {
				tmp_tgt_intervs.clear();
				tmp_metadata.clear();
				hint = chain_intervs.map_interval(*isrc_interv, tmp_tgt_intervs, hint, &tmp_metadata);

				// Clustering Strategy: cluster by source overlap (and same chain_id), keep heaviest cluster
				if (use_clustering && tmp_tgt_intervs.size() > 1) {
					// A. Zip intervals and metadata into candidates
					struct Candidate {
						GInterval interv;
						ChainMappingMetadata meta;
					};
					size_t n = tmp_tgt_intervs.size();
					vector<Candidate> candidates(n);
					for (size_t i = 0; i < n; ++i) {
						candidates[i].interv = tmp_tgt_intervs[i];
						if (i < tmp_metadata.size()) candidates[i].meta = tmp_metadata[i];
					}

					// B. Union-Find data structure for clustering
					vector<size_t> parent(n), rank_uf(n, 0);
					for (size_t i = 0; i < n; ++i) parent[i] = i;

					// Find with path compression
					std::function<size_t(size_t)> find_root = [&](size_t x) -> size_t {
						if (parent[x] != x) parent[x] = find_root(parent[x]);
						return parent[x];
					};

					// Union by rank
					auto unite = [&](size_t x, size_t y) {
						size_t rx = find_root(x), ry = find_root(y);
						if (rx == ry) return;
						if (rank_uf[rx] < rank_uf[ry]) std::swap(rx, ry);
						parent[ry] = rx;
						if (rank_uf[rx] == rank_uf[ry]) rank_uf[rx]++;
					};

					// C. Connect candidates with same chain_id
					unordered_map<int64_t, vector<size_t>> chain_id_groups;
					for (size_t i = 0; i < n; ++i) {
						chain_id_groups[candidates[i].meta.chain_id].push_back(i);
					}
					for (auto& kv : chain_id_groups) {
						vector<size_t>& group = kv.second;
						for (size_t i = 1; i < group.size(); ++i) {
							unite(group[0], group[i]);
						}
					}

					// D. Connect candidates that overlap in source space (sweep line)
					// Sort indices by source start
					vector<size_t> sorted_idx(n);
					for (size_t i = 0; i < n; ++i) sorted_idx[i] = i;
					sort(sorted_idx.begin(), sorted_idx.end(), [&](size_t a, size_t b) {
						return candidates[a].meta.start_src < candidates[b].meta.start_src;
					});

					// Sweep line to find overlaps
					int64_t max_end = -1;
					size_t max_end_idx = 0;
					for (size_t i = 0; i < n; ++i) {
						size_t idx = sorted_idx[i];
						if (candidates[idx].meta.start_src < max_end) {
							// Overlap with previous - connect to the one that set max_end
							unite(idx, max_end_idx);
						}
						if (candidates[idx].meta.end_src > max_end) {
							max_end = candidates[idx].meta.end_src;
							max_end_idx = idx;
						}
					}

					// E. Group candidates by their root (connected component)
					unordered_map<size_t, vector<size_t>> components;
					for (size_t i = 0; i < n; ++i) {
						components[find_root(i)].push_back(i);
					}

					// F. Calculate mass for each component and find the best
					vector<size_t> best_component_indices;
					int64_t best_mass = -1;
					int64_t best_min_start = INT64_MAX;  // For tie-breaking

					for (auto& kv : components) {
						vector<size_t>& comp = kv.second;

						// Sort component members by source start for mass calculation
						sort(comp.begin(), comp.end(), [&](size_t a, size_t b) {
							return candidates[a].meta.start_src < candidates[b].meta.start_src;
						});

						int64_t mass = 0;
						int64_t min_start = candidates[comp[0]].meta.start_src;  // First element after sort

						if (strat == STRAT_UNION) {
							// Union: count unique source bp
							int64_t union_end = -1;
							for (size_t idx : comp) {
								int64_t s = candidates[idx].meta.start_src;
								int64_t e = candidates[idx].meta.end_src;
								if (e > union_end) {
									mass += e - std::max(s, union_end);
									union_end = e;
								}
							}
						} else if (strat == STRAT_SUM) {
							// Sum: total of all segment lengths
							for (size_t idx : comp) {
								mass += candidates[idx].meta.end_src - candidates[idx].meta.start_src;
							}
						} else if (strat == STRAT_MAX) {
							// Max: largest single segment
							for (size_t idx : comp) {
								int64_t len = candidates[idx].meta.end_src - candidates[idx].meta.start_src;
								mass = std::max(mass, len);
							}
						}

						// Prefer higher mass, or earlier start position for ties
						if (mass > best_mass || (mass == best_mass && min_start < best_min_start)) {
							best_mass = mass;
							best_min_start = min_start;
							best_component_indices = comp;
						}
					}

					// G. Apply result: keep only the best component
					tmp_tgt_intervs.clear();
					tmp_metadata.clear();

					for (size_t idx : best_component_indices) {
						tmp_tgt_intervs.push_back(candidates[idx].interv);
						tmp_metadata.push_back(candidates[idx].meta);
					}
				}

				if (!tmp_tgt_intervs.empty()) {
					tgt_intervs1d.insert(tgt_intervs1d.end(), tmp_tgt_intervs.begin(), tmp_tgt_intervs.end());
					src_indices.insert(src_indices.end(), tmp_tgt_intervs.size(), iu.get_orig_interv_idx(*isrc_interv) + 1);
					for (size_t i = 0; i < tmp_metadata.size(); ++i) {
						chain_ids.push_back(tmp_metadata[i].chain_id);
						scores.push_back(tmp_metadata[i].score);
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

			for (size_t i = 0; i < order.size(); ++i) {
				size_t idx = order[i];
				const GInterval &interval = tgt_intervs1d[idx];
				int intervalID = src_indices[idx];
				int64_t chainID = chain_ids[idx];
				double score = scores[idx];

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
				} else {
					// Merge with previous (extend end)
					merged_intervs.back().end = interval.end;
				}
			}

			tgt_intervs1d = merged_intervs;
			src_indices = merged_indices;
			chain_ids = merged_chain_ids;
			scores = merged_scores;
		}

		// assemble the answer
		SEXP answer;
		unsigned num_interv_cols;

		if (!tgt_intervs1d.empty()) {
			// Add extra columns for intervalID, chain_id, and optionally score
			int extra_cols = include_metadata ? 3 : 2;
			answer = iu.convert_intervs(&tgt_intervs1d, GInterval::NUM_COLS + extra_cols);
			num_interv_cols = GInterval::NUM_COLS;
		} else if (!tgt_intervs2d.empty()) {
			answer = iu.convert_intervs(&tgt_intervs2d, GInterval2D::NUM_COLS + 1);
			num_interv_cols = GInterval2D::NUM_COLS;
		} else
			return R_NilValue;

		// Add intervalID, chain_id, and optionally score columns for 1D intervals
		if (!tgt_intervs1d.empty()) {
			SEXP rsrc_indices;
			SEXP rchain_ids;
			SEXP rscores = R_NilValue;
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

			if (include_metadata) {
				rprotect(rscores = RSaneAllocVector(REALSXP, scores.size()));
				for (size_t i = 0; i < scores.size(); ++i)
					REAL(rscores)[i] = scores[i];
				SET_STRING_ELT(col_names, num_interv_cols + 2, Rf_mkChar("score"));
				SET_VECTOR_ELT(answer, num_interv_cols + 2, rscores);
				runprotect(3); // rsrc_indices, rchain_ids, rscores
			} else {
				runprotect(2); // rsrc_indices, rchain_ids
			}
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
