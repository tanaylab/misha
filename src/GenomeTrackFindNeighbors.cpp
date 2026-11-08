/*
 * GenomeTrackFindNeighbors.cpp
 *
 *  Created on: Nov 9, 2010
 *      Author: hoichman
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "SegmentFinder.h"
#include "StatQuadTree.h"

using namespace std;
using namespace rdb;

// Neighbors are ordered by query, then by distance, then by the neighbor's coordinates, so which
// neighbors are chosen among those at the same distance, and their order, does not depend on the
// input order. Only neighbors with the same coordinates keep their input order (by row).
struct IntervNeighbor {
	int64_t        id1;
	int64_t        id2;
	int64_t        dist;
	int64_t        start2;
	int64_t        end2;

	IntervNeighbor(int64_t _id1, int64_t _id2, int64_t _dist, int64_t _start2 = 0, int64_t _end2 = 0) :
		id1(_id1), id2(_id2), dist(_dist), start2(_start2), end2(_end2) {}

	bool operator<(const IntervNeighbor &o) const {
		if (id1 != o.id1)
			return id1 < o.id1;
		if (llabs(dist) != llabs(o.dist))
			return llabs(dist) < llabs(o.dist);
		if (start2 != o.start2)
			return start2 < o.start2;
		if (end2 != o.end2)
			return end2 < o.end2;
		return id2 < o.id2;
	}
};

struct IntervNeighbor2D {
	int64_t        id1;
	int64_t        id2;
	int64_t        dist1;
	int64_t        dist2;
	int64_t        start1, end1, start2, end2;  // of the neighbor

	IntervNeighbor2D(int64_t _id1, int64_t _id2, int64_t _dist1, int64_t _dist2, const Rectangle &rect = Rectangle(0, 0, 0, 0)) :
		id1(_id1), id2(_id2), dist1(_dist1), dist2(_dist2), start1(rect.x1), end1(rect.x2), start2(rect.y1), end2(rect.y2) {}

	bool operator<(const IntervNeighbor2D &o) const {
		if (id1 != o.id1)
			return id1 < o.id1;
		if (llabs(dist1 + dist2) != llabs(o.dist1 + o.dist2))
			return llabs(dist1 + dist2) < llabs(o.dist1 + o.dist2);
		if (start1 != o.start1)
			return start1 < o.start1;
		if (end1 != o.end1)
			return end1 < o.end1;
		if (start2 != o.start2)
			return start2 < o.start2;
		if (end2 != o.end2)
			return end2 < o.end2;
		return id2 < o.id2;
	}
};

// The intervals2 of one chromosome at distance 0 from a query - those that overlap or touch it,
// start <= query end and end >= query start - in the order neighbors are kept: (start, end, row).
// A query that overlaps K of them would otherwise walk all K in the nearest-neighbor iterator to
// choose the maxneighbors with the smallest coordinates; here each one kept costs O(log n), through
// a tree of the largest end over the sorted intervals.
class ZeroDistNeighbors {
public:
	struct Item {
		int64_t start;
		int64_t end;
		int64_t row;
		bool operator<(const Item &o) const {
			return start < o.start || (start == o.start && (end < o.end || (end == o.end && row < o.row)));
		}
	};

	void clear() { m_items.clear(); m_maxend.clear(); m_cap = 0; }
	bool empty() const { return m_items.empty() && !m_cap; }

	void build(vector<Item> &&items) {
		m_items = std::move(items);
		sort(m_items.begin(), m_items.end());
		m_cap = 1;
		while (m_cap < m_items.size())
			m_cap *= 2;
		m_maxend.assign(2 * m_cap, numeric_limits<int64_t>::min());
		for (size_t i = 0; i < m_items.size(); ++i)
			m_maxend[m_cap + i] = m_items[i].end;
		for (size_t i = m_cap - 1; i > 0; --i)
			m_maxend[i] = max(m_maxend[2 * i], m_maxend[2 * i + 1]);
	}

	// Calls take(item) for the first n intervals at distance 0 from [qstart, qend), in order
	template <class F> void first(int64_t qstart, int64_t qend, int n, F take) const {
		for (size_t i = find(1, 0, m_cap, 0, qstart); n > 0 && i < m_items.size() && m_items[i].start <= qend; i = find(1, 0, m_cap, i + 1, qstart), --n)
			take(m_items[i]);
	}

private:
	vector<Item>    m_items;
	vector<int64_t> m_maxend;
	size_t          m_cap = 0;

	// the first index >= from whose end is >= qstart, or the number of items
	size_t find(size_t node, size_t lo, size_t hi, size_t from, int64_t qstart) const {
		if (hi <= from || m_maxend[node] < qstart)
			return m_items.size();
		if (hi - lo == 1)
			return lo;
		size_t mid = lo + (hi - lo) / 2;
		size_t res = find(2 * node, lo, mid, from, qstart);
		return res < m_items.size() ? res : find(2 * node + 1, mid, hi, from, qstart);
	}
};

// The same for the intervals2 of one chromosome pair: those at dist1 = dist2 = 0 from a query
// (overlapping or touching it on both axes), in the order (start1, end1, start2, end2, row). The
// tree holds, for each run of the sorted rectangles, the largest end1, the smallest start2 and the
// largest end2, and a run that cannot hold a match is skipped whole. The rectangles do not overlap
// (verify_no_overlaps), so the runs a query does not prune are mostly ones that it overlaps.
class ZeroDistNeighbors2D {
public:
	struct Item {
		int64_t start1, end1, start2, end2;
		int64_t row;
		bool operator<(const Item &o) const {
			if (start1 != o.start1) return start1 < o.start1;
			if (end1 != o.end1) return end1 < o.end1;
			if (start2 != o.start2) return start2 < o.start2;
			if (end2 != o.end2) return end2 < o.end2;
			return row < o.row;
		}
	};

	void clear() { m_items.clear(); m_nodes.clear(); m_cap = 0; }
	bool empty() const { return m_items.empty() && !m_cap; }

	void build(vector<Item> &&items) {
		m_items = std::move(items);
		sort(m_items.begin(), m_items.end());
		m_cap = 1;
		while (m_cap < m_items.size())
			m_cap *= 2;
		m_nodes.assign(2 * m_cap, Node());
		for (size_t i = 0; i < m_items.size(); ++i)
			m_nodes[m_cap + i] = Node{m_items[i].end1, m_items[i].start2, m_items[i].end2};
		for (size_t i = m_cap - 1; i > 0; --i) {
			const Node &l = m_nodes[2 * i], &r = m_nodes[2 * i + 1];
			m_nodes[i] = Node{max(l.max_end1, r.max_end1), min(l.min_start2, r.min_start2), max(l.max_end2, r.max_end2)};
		}
	}

	// Calls take(item) for the first n rectangles at distance 0 from the query, in order
	template <class F> void first(const Rectangle &q, int n, F take) const {
		for (size_t i = find(1, 0, m_cap, 0, q); n > 0 && i < m_items.size(); i = find(1, 0, m_cap, i + 1, q), --n)
			take(m_items[i]);
	}

private:
	struct Node {
		int64_t max_end1 = numeric_limits<int64_t>::min();
		int64_t min_start2 = numeric_limits<int64_t>::max();
		int64_t max_end2 = numeric_limits<int64_t>::min();
	};

	vector<Item> m_items;
	vector<Node> m_nodes;
	size_t       m_cap = 0;

	// the first index >= from of a rectangle at distance 0 from q, or the number of items
	size_t find(size_t node, size_t lo, size_t hi, size_t from, const Rectangle &q) const {
		const Node &nd = m_nodes[node];
		if (hi <= from || lo >= m_items.size() || m_items[lo].start1 > q.x2 || nd.max_end1 < q.x1 || nd.min_start2 > q.y2 || nd.max_end2 < q.y1)
			return m_items.size();
		if (hi - lo == 1)
			return lo;
		size_t mid = lo + (hi - lo) / 2;
		size_t res = find(2 * node, lo, mid, from, q);
		return res < m_items.size() ? res : find(2 * node + 1, mid, hi, from, q);
	}
};

// Helper function to extract strand values from R data frame
int get_strand_value(const SEXP intervals_df, int64_t row_index, bool has_strand_col) {
    if (!has_strand_col) {
        return 1;  // Default to + strand
    }
    
    // Find strand column index
    SEXP colnames = Rf_getAttrib(intervals_df, R_NamesSymbol);
    int strand_col_idx = -1;
    for (int i = 0; i < Rf_length(colnames); i++) {
        if (strcmp(CHAR(STRING_ELT(colnames, i)), "strand") == 0) {
            strand_col_idx = i;
            break;
        }
    }
    
    if (strand_col_idx == -1) {
        return 1;
    }
    
    SEXP strand_col = VECTOR_ELT(intervals_df, strand_col_idx);
    
    // Bounds check
    if (row_index < 0 || row_index >= Rf_length(strand_col)) {
        return 1;  // Default if out of bounds
    }
    
    const char *strand_str = NULL;
    if (Rf_isFactor(strand_col)) {
        // A factor is an integer vector of level codes: read the level, not the code
        int code = INTEGER(strand_col)[row_index];
        SEXP levels = Rf_getAttrib(strand_col, R_LevelsSymbol);
        if (code != NA_INTEGER && code >= 1 && code <= Rf_length(levels))
            strand_str = CHAR(STRING_ELT(levels, code - 1));
    } else if (Rf_isString(strand_col)) {
        strand_str = CHAR(STRING_ELT(strand_col, row_index));
    }

    if (Rf_isFactor(strand_col) || Rf_isString(strand_col)) {
        // Handle character or factor strand column (e.g., "+", "-")
        if (strand_str != NULL) {
            if (strcmp(strand_str, "+") == 0 || strcmp(strand_str, "1") == 0) {
                return 1;
            } else if (strcmp(strand_str, "-") == 0 || strcmp(strand_str, "-1") == 0) {
                return -1;
            } else {
                return 1;  // Default to + strand
            }
        } else {
            return 1;  // Default if NA
        }
    } else if (Rf_isInteger(strand_col)) {
        int val = INTEGER(strand_col)[row_index];
        if (val != -1 && val != 1) {
            rdb::verror("Invalid strand value %d at row %lld. Strand must be -1 or 1.", val, (long long)row_index + 1);
        }
        return val;
    } else if (Rf_isReal(strand_col)) {
        int val = (int)REAL(strand_col)[row_index];
        if (val != -1 && val != 1) {
            rdb::verror("Invalid strand value %d at row %lld. Strand must be -1 or 1.", val, (long long)row_index + 1);
        }
        return val;
    } else if (Rf_isLogical(strand_col)) {
        // Handle logical strand column
        int logical_val = LOGICAL(strand_col)[row_index];
        return logical_val ? 1 : -1;
    } else {
        return 1;  // Default if can't parse
    }
}

// Helper function to check if intervals data frame has strand column
bool has_strand_column(const SEXP intervals_df) {
    if (intervals_df == R_NilValue) {
        return false;
    }
    
    SEXP colnames = Rf_getAttrib(intervals_df, R_NamesSymbol);
    if (colnames == R_NilValue) {
        return false;
    }
    
    for (int i = 0; i < Rf_length(colnames); i++) {
        if (strcmp(CHAR(STRING_ELT(colnames, i)), "strand") == 0) {
            return true;
        }
    }
    return false;
}

// A distance limit argument (mindist, maxdist, ...): a single number, not NA. -Inf and Inf
// (and anything beyond +-1e18, far past any genome) mean no limit; casting them to int64_t
// is undefined.
static int64_t get_dist_arg(SEXP x, const char *name)
{
	if (!is_num_vector(x) || Rf_length(x) != 1)
		verror("%s argument is not numeric", name);
	double v = num_elt(x, 0);
	if (std::isnan(v))
		verror("%s argument is NA", name);
	return (int64_t)std::max(-1e18, std::min(1e18, v));
}

extern "C" {

SEXP gfind_neighbors(SEXP _intervs1, SEXP _intervs2, SEXP _maxneighbors, SEXP _distrange_start, SEXP _distrange_end,
					 SEXP _distrange_start1, SEXP _distrange_end1, SEXP _distrange_start2, SEXP _distrange_end2,
					 SEXP _na_if_notfound, SEXP _report_progress, SEXP _use_intervals1_strand, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;
		GIntervals intervals1d[2];
		GIntervals2D intervals2d[2];

		IntervUtils iu(_envir);
		_intervs1 = iu.convert_rintervs(_intervs1, &intervals1d[0], &intervals2d[0]);
		_intervs2 = iu.convert_rintervs(_intervs2, &intervals1d[1], &intervals2d[1]);

		unsigned type_mask1 = iu.get_rintervs_type_mask(_intervs1);
		unsigned type_mask2 = iu.get_rintervs_type_mask(_intervs2);

		if (type_mask1 == (IntervUtils::INTERVS1D | IntervUtils::INTERVS2D) ||
			type_mask2 == (IntervUtils::INTERVS1D | IntervUtils::INTERVS2D))
			verror("Intervals cannot be of dual intervals type");

		if (type_mask1 != type_mask2) 
			verror("Cannot intermix 1D and 2D intervals");

		if (!is_num_vector(_maxneighbors) || Rf_length(_maxneighbors) != 1)
			verror("maxneighbors argument is not numeric");
		int maxneighbors = Rf_isInteger(_maxneighbors) ? INTEGER(_maxneighbors)[0] : (int)REAL(_maxneighbors)[0];
		if (maxneighbors < 1) 
			verror("maxneighbors must be greater or equal to 1");

		if (!is_true_or_false(_na_if_notfound))
			verror("na.if.notfound argument is not logical");
		bool na_if_notfound = LOGICAL(_na_if_notfound)[0];

		if (!is_true_or_false(_use_intervals1_strand))
			verror("use_intervals1_strand argument is not logical");
		bool use_intervals1_strand = LOGICAL(_use_intervals1_strand)[0];


		Progress_reporter progress;

		if (type_mask1 & IntervUtils::INTERVS1D) {
			int64_t distance_start = get_dist_arg(_distrange_start, "mindist");
			int64_t distance_end = get_dist_arg(_distrange_end, "maxdist");

			if (distance_start > distance_end)
				verror("mindist exceeds maxdist");

			int64_t maxdist = max(llabs(distance_start), llabs(distance_end));
			int64_t mindist = min(llabs(distance_start), llabs(distance_end));

			if (intervals1d[0].empty() || (intervals1d[1].empty() && !na_if_notfound))
				return R_NilValue;

			vector<IntervNeighbor> result;

			intervals1d[0].sort();
			// Ties do not depend on this order (see IntervNeighbor); stable to keep it the same on
			// every platform all the same.
			stable_sort(intervals1d[1].begin(), intervals1d[1].end(), GIntervalsFetcher1D::compare_by_start_coord);

			int chromid = -1;
			SegmentFinder<GInterval> segment_finder;
			SegmentFinder<GInterval>::NNIterator inn(&segment_finder);
			ZeroDistNeighbors zero_dist;  // built for chromid when a query first needs it

			if (LOGICAL(_report_progress)[0])
				progress.init(intervals1d[0].size(), 100);

			for (GIntervals::const_iterator iinterv = intervals1d[0].begin(); iinterv < intervals1d[0].end(); ++iinterv) {
				// Get strand for current query interval
				int query_strand = 1;  // default
				if (use_intervals1_strand) {
					bool intervals1_has_strand = has_strand_column(_intervs1);
					query_strand = get_strand_value(_intervs1, iu.get_orig_interv_idx(*iinterv), intervals1_has_strand);
				}

				if (iinterv->chromid != chromid) {
					chromid = iinterv->chromid;
					segment_finder.reset(0, iu.get_chromkey().get_chrom_size(chromid));

					for (intervals1d[1].begin_chrom_iter(chromid); !intervals1d[1].isend_chrom(); intervals1d[1].next_in_chrom()) {
						const GInterval &interv = intervals1d[1].cur_interval();
						segment_finder.insert(interv);
					}
					zero_dist.clear();
				}

				if ((distance_start > 0 && distance_end > 0) || (distance_start < 0 && distance_end < 0)) {
					Segment excluded_area(iinterv->start - mindist + 1, iinterv->end + mindist - 1);
					inn.begin(*iinterv, excluded_area);
				} else
					inn.begin(*iinterv);

				int num_neighbors = 0;
				int64_t cutoff = 0;  // |dist| of the maxneighbors-th neighbor found
				size_t first = result.size();
				bool zero_dist_ties = false;
				while (!inn.is_end()) {
					int64_t dist;
					
					if (use_intervals1_strand) {
						// Use query (intervals1) strand for distance calculation
						dist = iinterv->dist2interv_with_query_strand(*inn, query_strand);
					} else {
						// Original behavior - use target (intervals2) strand if available
						dist = iinterv->dist2interv(*inn);
					}
					
					if (llabs(dist) > maxdist)
						break;
					// The iterator yields the neighbors by |dist|. Once maxneighbors are found, those
					// at the distance of the last one are taken too, and the ties are cut below.
					if (num_neighbors >= maxneighbors && llabs(dist) > cutoff)
						break;
					if (dist >= distance_start && dist <= distance_end) {
						result.push_back(IntervNeighbor(iu.get_orig_interv_idx(*iinterv), iu.get_orig_interv_idx(*inn), dist, inn->start, inn->end));
						if (++num_neighbors == maxneighbors)
							cutoff = llabs(dist);
						// the rows kept: tied candidates beyond maxneighbors are cut below
						iu.verify_max_data_size(first + min(num_neighbors, maxneighbors));
						// maxneighbors at distance 0: the ties can be all the intervals2 that overlap
						// the query, taken below in their order instead of walked here
						if (num_neighbors == maxneighbors && cutoff == 0) {
							zero_dist_ties = true;
							break;
						}
					}
					inn.next();
				}

				// Distance 0 is in [mindist, maxdist] here (it was reached), so every interval at
				// distance 0 is a candidate, and the ones kept are the first maxneighbors in
				// (start, end, row) order. Their distance is 0 with any strand.
				if (zero_dist_ties) {
					if (zero_dist.empty()) {
						vector<ZeroDistNeighbors::Item> items;
						for (intervals1d[1].begin_chrom_iter(chromid); !intervals1d[1].isend_chrom(); intervals1d[1].next_in_chrom()) {
							const GInterval &interv = intervals1d[1].cur_interval();
							items.push_back({interv.start, interv.end, (int64_t)iu.get_orig_interv_idx(interv)});
						}
						zero_dist.build(std::move(items));
					}
					result.erase(result.begin() + first, result.end());
					int64_t id1 = iu.get_orig_interv_idx(*iinterv);
					zero_dist.first(iinterv->start, iinterv->end, maxneighbors, [&](const ZeroDistNeighbors::Item &item) {
						result.push_back(IntervNeighbor(id1, item.row, 0, item.start, item.end));
						iu.verify_max_data_size(result.size());
					});
				}

				// neighbors tied at the cutoff distance: keep those with the smallest coordinates
				if (num_neighbors > maxneighbors) {
					sort(result.begin() + first, result.end());
					result.erase(result.begin() + first + maxneighbors, result.end());
					num_neighbors = maxneighbors;
				}

				if (na_if_notfound && !num_neighbors) {
					result.push_back(IntervNeighbor(iu.get_orig_interv_idx(*iinterv), -1, 0));
					iu.verify_max_data_size(result.size());
				}

				if (LOGICAL(_report_progress)[0])
					progress.report(1);
				check_interrupt();
			}
			if (LOGICAL(_report_progress)[0])
				progress.report_last();

			if (result.empty())
				return R_NilValue;

			sort(result.begin(), result.end());

			SEXP answer = iu.create_data_frame(result.size(), Rf_length(_intervs1) + Rf_length(_intervs2) + 1);
			SEXP rdists;
			vector<SEXP> src_cols1;
			vector<SEXP> src_cols2;
			vector<SEXP> tgt_cols;

			iu.define_data_frame_cols(_intervs1, src_cols1, answer, tgt_cols, 0);
			iu.define_data_frame_cols(_intervs2, src_cols2, answer, tgt_cols, Rf_length(_intervs1));
			rprotect(rdists = RSaneAllocVector(REALSXP, result.size()));

			for (uint64_t i = 0; i < result.size(); ++i) {
				const IntervNeighbor &r = result[i];
				iu.copy_data_frame_row(src_cols1, r.id1, tgt_cols, i, 0);

				if (r.id2 >= 0) {
					iu.copy_data_frame_row(src_cols2, r.id2, tgt_cols, i, Rf_length(_intervs1));
					REAL(rdists)[i] = (double)r.dist;
				} else {
					for (int j = 0; j < Rf_length(_intervs2); ++j)
						iu.set_data_frame_val_nan(tgt_cols, i, Rf_length(_intervs1) + j);
					REAL(rdists)[i] = NA_REAL;
				}
			}
			SET_VECTOR_ELT(answer, Rf_length(_intervs1) + Rf_length(_intervs2), rdists);
			SEXP colnames = rprotect_ptr(Rf_getAttrib(answer, R_NamesSymbol));
			SET_STRING_ELT(colnames, Rf_length(_intervs1) + Rf_length(_intervs2), Rf_mkChar("dist"));
			// colnames only - the last thing protected. create_data_frame() and
			// define_data_frame_cols() left their own objects below it, and "answer"
			// must stay protected among them.
			runprotect(1);

			return answer;
		} else {   // 2D
			// For now, error if use_intervals1_strand with 2D
			if (use_intervals1_strand) {
				verror("use_intervals1_strand not yet supported for 2D intervals");
			}

			int64_t distance_start1 = get_dist_arg(_distrange_start1, "mindist1");
			int64_t distance_end1 = get_dist_arg(_distrange_end1, "maxdist1");

			if (distance_start1 > distance_end1)
				verror("mindist1 exceeds maxdist1");

			int64_t distance_start2 = get_dist_arg(_distrange_start2, "mindist2");
			int64_t distance_end2 = get_dist_arg(_distrange_end2, "maxdist2");

			if (distance_start2 > distance_end2)
				verror("mindist2 exceeds maxdist2");

			if (intervals2d[0].empty() || (intervals2d[1].empty() && !na_if_notfound) || distance_end1 < 0 || distance_end2 < 0)
				return R_NilValue;

			vector<IntervNeighbor2D> result;

			intervals2d[0].sort();
			intervals2d[1].sort();
			intervals2d[1].verify_no_overlaps(iu.get_chromkey());

			int chromid1 = -1;
			int chromid2 = -1;
			StatQuadTree<Rectangle_val<uint64_t>, uint64_t> qtree;
			StatQuadTree<Rectangle_val<uint64_t>, uint64_t>::NNIterator inn(&qtree);
			ZeroDistNeighbors2D zero_dist;  // built for the chromosome pair when a query first needs it

			if (LOGICAL(_report_progress)[0])
				progress.init(intervals2d[0].size(), 100);

			for (GIntervals2D::const_iterator iinterv = intervals2d[0].begin(); iinterv < intervals2d[0].end(); ++iinterv) {
				if (iinterv->chromid1() != chromid1 || iinterv->chromid2() != chromid2) {
					chromid1 = iinterv->chromid1();
					chromid2 = iinterv->chromid2();
					qtree.reset(0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2));

					for (intervals2d[1].begin_chrom_iter(chromid1, chromid2); !intervals2d[1].isend_chrom(); intervals2d[1].next_in_chrom()) {
						const GInterval2D &interv = intervals2d[1].cur_interval();
						StatQuadTree<Rectangle_val<uint64_t>, uint64_t>::ValueType rect(interv, iu.get_orig_interv_idx(interv));
						qtree.insert(rect);
					}
					zero_dist.clear();
				}

				if (distance_start1 > 0 && distance_start2 > 0) {
					Rectangle excluded_area(iinterv->x1 - distance_start1, iinterv->y1 - distance_start2,
											iinterv->x2 + distance_start1, iinterv->y2 + distance_start2);
					inn.begin(*iinterv, excluded_area);
				} else
					inn.begin(*iinterv);

				int num_neighbors = 0;
				int64_t cutoff = 0;  // the largest dist1 + dist2 of the first maxneighbors neighbors found
				size_t first = result.size();
				bool zero_dist_ties = false;
				while (!inn.is_end()) {
					int64_t dist1 = iinterv->xdist(*inn);
					int64_t dist2 = iinterv->ydist(*inn);

					if (dist1 > distance_end1 && dist2 > distance_end2) 
						break;

					// The iterator orders the neighbors by dist1 + dist2 plus one for each axis on which
					// the rectangles do not overlap (manhattan_dist with touch_is_at_dist_one), so a
					// neighbor found later can be nearer by dist1 + dist2 than one found earlier, by up
					// to 2. Once maxneighbors are found, every neighbor at dist1 + dist2 <= cutoff comes
					// before that order passes cutoff + 2; the nearest are chosen from them below.
					if (num_neighbors >= maxneighbors && iinterv->manhattan_dist(*inn, true) > cutoff + 2)
						break;

					if (dist1 >= distance_start1 && dist1 <= distance_end1 && dist2 >= distance_start2 && dist2 <= distance_end2) {
						result.push_back(IntervNeighbor2D(iu.get_orig_interv_idx(*iinterv), inn->v, dist1, dist2, *inn));
						if (num_neighbors < maxneighbors)
							cutoff = max(cutoff, dist1 + dist2);
						num_neighbors++;
						// the rows kept: candidates beyond maxneighbors are cut below
						iu.verify_max_data_size(first + min(num_neighbors, maxneighbors));
						// maxneighbors at distance 0: taken below in their order, see the 1D case
						if (num_neighbors == maxneighbors && cutoff == 0) {
							zero_dist_ties = true;
							break;
						}
					}

					inn.next();
				}

				// Distance 0 on both axes is in range here (it was reached), so every rectangle at it
				// is a candidate, and the ones kept are the first maxneighbors in their order
				if (zero_dist_ties) {
					if (zero_dist.empty()) {
						vector<ZeroDistNeighbors2D::Item> items;
						for (intervals2d[1].begin_chrom_iter(chromid1, chromid2); !intervals2d[1].isend_chrom(); intervals2d[1].next_in_chrom()) {
							const GInterval2D &interv = intervals2d[1].cur_interval();
							items.push_back({interv.x1, interv.x2, interv.y1, interv.y2, (int64_t)iu.get_orig_interv_idx(interv)});
						}
						zero_dist.build(std::move(items));
					}
					result.erase(result.begin() + first, result.end());
					int64_t id1 = iu.get_orig_interv_idx(*iinterv);
					zero_dist.first(*iinterv, maxneighbors, [&](const ZeroDistNeighbors2D::Item &item) {
						result.push_back(IntervNeighbor2D(id1, item.row, 0, 0, Rectangle(item.start1, item.start2, item.end1, item.end2)));
						iu.verify_max_data_size(result.size());
					});
				}

				// keep the maxneighbors nearest by dist1 + dist2, ties by the smallest coordinates
				if (num_neighbors > maxneighbors) {
					sort(result.begin() + first, result.end());
					result.erase(result.begin() + first + maxneighbors, result.end());
					num_neighbors = maxneighbors;
				}

				if (na_if_notfound && !num_neighbors) {
					result.push_back(IntervNeighbor2D(iu.get_orig_interv_idx(*iinterv), -1, 0, 0));
					iu.verify_max_data_size(result.size());
				}

				if (LOGICAL(_report_progress)[0])
					progress.report(1);
				check_interrupt();
			}
			if (LOGICAL(_report_progress)[0])
				progress.report_last();

			if (result.empty())
				return R_NilValue;

			sort(result.begin(), result.end());

			SEXP answer = iu.create_data_frame(result.size(), Rf_length(_intervs1) + Rf_length(_intervs2) + 2);
			SEXP rdists1, rdists2;
			vector<SEXP> src_cols1;
			vector<SEXP> src_cols2;
			vector<SEXP> tgt_cols;

			iu.define_data_frame_cols(_intervs1, src_cols1, answer, tgt_cols, 0);
			iu.define_data_frame_cols(_intervs2, src_cols2, answer, tgt_cols, Rf_length(_intervs1));
			rprotect(rdists1 = RSaneAllocVector(REALSXP, result.size()));
			rprotect(rdists2 = RSaneAllocVector(REALSXP, result.size()));

			for (uint64_t i = 0; i < result.size(); ++i) {
				const IntervNeighbor2D &r = result[i];
				iu.copy_data_frame_row(src_cols1, r.id1, tgt_cols, i, 0);

				if (r.id2 >= 0) {
					iu.copy_data_frame_row(src_cols2, r.id2, tgt_cols, i, Rf_length(_intervs1));
					REAL(rdists1)[i] = (double)r.dist1;
					REAL(rdists2)[i] = (double)r.dist2;
				} else {
					for (int j = 0; j < Rf_length(_intervs2); ++j)
						iu.set_data_frame_val_nan(tgt_cols, i, Rf_length(_intervs1) + j);
					REAL(rdists1)[i] = NA_REAL;
					REAL(rdists2)[i] = NA_REAL;
				}
			}
			SET_VECTOR_ELT(answer, Rf_length(_intervs1) + Rf_length(_intervs2), rdists1);
			SET_VECTOR_ELT(answer, Rf_length(_intervs1) + Rf_length(_intervs2) + 1, rdists2);

			SEXP colnames = rprotect_ptr(Rf_getAttrib(answer, R_NamesSymbol));
			SET_STRING_ELT(colnames, Rf_length(_intervs1) + Rf_length(_intervs2), Rf_mkChar("dist1"));
			SET_STRING_ELT(colnames, Rf_length(_intervs1) + Rf_length(_intervs2) + 1, Rf_mkChar("dist2"));
			// colnames only - see the 1D branch above.
			runprotect(1);

			return answer;
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
