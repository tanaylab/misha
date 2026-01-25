/*
 * rdbinterval.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <cmath>
#include <cctype>
#include <limits>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <utility>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "TrackExpressionScanner.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "GTrackIntervalsFetcher1D.h"
#include "GTrackIntervalsFetcher2D.h"

using namespace std;

#include "strutil.h"

using namespace rdb;

const char *IntervalPval::COL_NAMES[IntervalPval::NUM_COLS] = { "chrom", "start", "end", "pval" };

const char *ChainInterval::COL_NAMES[ChainInterval::NUM_COLS] = { "chrom", "start", "end", "strand", "chromsrc", "startsrc", "endsrc", "strandsrc", "chain_id", "score" };

IntervUtils::IntervUtils(SEXP envir) :
	m_envir(envir),
	m_num_planned_kids(0),
	m_kid_intervals1d(NULL),
	m_kid_intervals2d(NULL),
	m_config(envir),
	m_validator(m_chrom_key),
	m_chain_converter(*this),
	m_interval_converter(*this)
{
	m_kids_intervals1d.clear();
	m_kids_intervals2d.clear();

    {
        m_allgenome = find_in_misha(m_envir, "ALLGENOME");
    }

	if (Rf_isNull(m_allgenome))
		verror("ALLGENOME variable does not exist");

	if (!Rf_isVector(m_allgenome) || Rf_length(m_allgenome) != 2)
		verror("ALLGENOME variable has invalid type");

	SEXP chroms = VECTOR_ELT(get_rallgenome1d(), GInterval::CHROM);
	SEXP chrom_sizes = VECTOR_ELT(get_rallgenome1d(), GInterval::END);
	SEXP chrom_levels = Rf_getAttrib(chroms, R_LevelsSymbol);
	unsigned num_intervals = (unsigned)Rf_length(chroms);

	for (unsigned i = 0; i < num_intervals; i++) {
		const char *chrom = Rf_isString(chroms) ? CHAR(STRING_ELT(chroms, i)) : CHAR(STRING_ELT(chrom_levels, INTEGER(chroms)[i] - 1));
		double chrom_size = Rf_isReal(chrom_sizes) ? REAL(chrom_sizes)[i] : INTEGER(chrom_sizes)[i];
		try {
			m_chrom_key.add_chrom(chrom, (uint64_t)chrom_size);
		} catch (TGLException &e) {
			verror("Reading ALLGENOME: %s", e.msg());
		}
	}

	// Populate chromosome aliases from R CHROM_ALIAS map
	SEXP chrom_alias = find_in_misha(m_envir, "CHROM_ALIAS");
	if (!Rf_isNull(chrom_alias) && Rf_isVector(chrom_alias)) {
		SEXP alias_names = Rf_getAttrib(chrom_alias, R_NamesSymbol);
		if (!Rf_isNull(alias_names)) {
			int n_aliases = Rf_length(chrom_alias);
			for (int i = 0; i < n_aliases; i++) {
				const char *alias = CHAR(STRING_ELT(alias_names, i));
				const char *canonical = CHAR(STRING_ELT(chrom_alias, i));

				// Find the chromid for the canonical chromosome name
				try {
					int chrom_id = m_chrom_key.chrom2id(canonical);
					m_chrom_key.add_chrom_alias(alias, chrom_id);
				} catch (TGLException &) {
					// Silently skip aliases that point to non-existent chromosomes
				}
			}
		}
	}

    GenomeTrack::set_rnd_func(unif_rand);
}

IntervUtils::~IntervUtils()
{
	for (vector<GIntervalsFetcher1D *>::iterator iinterv = m_kids_intervals1d.begin(); iinterv != m_kids_intervals1d.end(); ++iinterv) 
		delete *iinterv;
	for (vector<GIntervalsFetcher2D *>::iterator iinterv = m_kids_intervals2d.begin(); iinterv != m_kids_intervals2d.end(); ++iinterv) 
		delete *iinterv;
}

bool IntervUtils::track_exists(const char *track_name)
{
	SEXP all_track_names = R_NilValue;
	SEXPCleaner all_track_names_cleaner(all_track_names);

	rprotect(all_track_names = find_in_misha(get_env(), "GTRACKS"));
	if (Rf_isString(all_track_names)) {
		for (int i = 0; i < Rf_length(all_track_names); ++i) {
			if (!strcmp(track_name, CHAR(STRING_ELT(all_track_names, i))))
				return true;
		}
	}
	return false;
}

void IntervUtils::get_all_genome_intervs(GIntervals &intervals) const
{
	intervals.clear();
	convert_rintervs(get_rallgenome1d(), &intervals, NULL);
	intervals.sort();
}

void IntervUtils::get_all_genome_intervs(GIntervals2D &intervals) const
{
	intervals.clear();
	convert_rintervs(get_rallgenome2d(), NULL, &intervals);
	intervals.sort();
}

GIntervalsFetcher1D *IntervUtils::get_kid_intervals1d()
{
	if (!m_kid_intervals1d && !m_kids_intervals1d.empty()) 
		return m_kids_intervals1d[RdbInitializer::get_kid_idx()];

	return NULL;
}

GIntervalsFetcher2D *IntervUtils::get_kid_intervals2d()
{
	if (!m_kid_intervals2d && !m_kids_intervals2d.empty())
		return m_kids_intervals2d[RdbInitializer::get_kid_idx()];

	return NULL;
}



DiagonalBand IntervUtils::convert_band(SEXP rband)
{
	if (Rf_isNull(rband))
		return DiagonalBand();

	if ((!Rf_isReal(rband) && !Rf_isInteger(rband)) || Rf_length(rband) != 2)
		verror("Invalid format of band argument");

	int d1 = Rf_isReal(rband) ? (int)(REAL(rband)[0] > 0 ? REAL(rband)[0] + 0.5 : REAL(rband)[0] - 0.5) : INTEGER(rband)[0];
	int d2 = Rf_isReal(rband) ? (int)(REAL(rband)[1] > 0 ? REAL(rband)[1] + 0.5 : REAL(rband)[1] - 0.5) : INTEGER(rband)[1];

	if (d1 >= d2)
		verror("Invalid band argument: distance1 exceeds distance2");

	return DiagonalBand(d1, d2);
}



bool IntervUtils::is_1d_iterator(SEXP rtrack_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP riterator)
{
	TrackExprScanner scanner(*this);
	TrackExpressionIteratorBase *itr = scanner.create_expr_iterator(rtrack_exprs, scope1d, scope2d, riterator, R_NilValue, false);
	return itr->is_1d();
}

int IntervUtils::prepare4multitasking(SEXP track_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP iterator_policy, SEXP band)
{
	TrackExprScanner scanner(*this);

	scanner.check(track_exprs, scope1d, scope2d, iterator_policy, band);
	if (scanner.get_iterator()->is_1d()) {
		if (scope2d && dynamic_cast<GIntervals2D *>(scope2d))
			((GIntervals2D *)scope2d)->clear();
	} else {
		if (scope1d && dynamic_cast<GIntervals *>(scope1d)) 
			((GIntervals *)scope1d)->clear();
	}

	return prepare4multitasking(scope1d, scope2d);
}

int IntervUtils::prepare4multitasking(GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d)
{
	if (scope1d && !scope1d->size() && scope2d && !scope2d->size()) 
		return 0;

	if ((scope1d && scope1d->size() && scope2d && scope2d->size()) || (!scope1d && !scope2d)) 
		verror("Cannot determine iterator policy");

	int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
	int max_num_pids = min(get_max_processes2core() * num_cores, get_max_processes());
	int num_chroms = scope1d && scope1d->size() ? scope1d->num_chroms() : scope2d->num_chrom_pairs();

	m_kids_intervals1d.clear();
	m_kids_intervals2d.clear();

	if (scope1d && scope1d->size()) {
		int64_t range = scope1d->range();
		int64_t kid_range = 0;
		int num_avail_kids = min(max_num_pids, num_chroms) - 1;
		int num_remaining_chroms = num_chroms - 1;

		// GIntervals::range(chromid) and GIntervals::create_masked_copy() have complexity of O(n) (and only O(1) for GIntervalsBigSet).
		// Calling these functions for each chromosome might be painful - O(m x n) (m = num chromosomes)
		// So we are forced to write a designated code for GIntervals for the sake of efficiency. The total complexity for GIntervals is O(n).
		if (dynamic_cast<GIntervals *>(scope1d)) {
			GIntervals *intervals = (GIntervals *)scope1d;

			m_kids_intervals1d.push_back(new GIntervals());

			for (GIntervals::const_iterator iinterv = intervals->begin(); iinterv != intervals->end(); ++iinterv) {
				if (iinterv != intervals->begin() && iinterv->chromid != (iinterv - 1)->chromid && num_avail_kids) {
					int64_t chrom_range = 0;
					for (GIntervals::const_iterator iinterv2 = iinterv; iinterv2 != intervals->end(); ++iinterv2) {
						if (iinterv->chromid != iinterv2->chromid) 
							break;
						chrom_range += iinterv2->range();
					}

					// should we allocate a new process?
					if ((uint64_t)kid_range > get_min_scope4process() && (uint64_t)range > get_min_scope4process() && (uint64_t)(kid_range + chrom_range) >= (uint64_t)(range / num_avail_kids)) {
						--num_avail_kids;
						kid_range = 0;
						m_kids_intervals1d.push_back(new GIntervals());
					}

					--num_remaining_chroms;
					num_avail_kids = min(num_avail_kids, num_remaining_chroms);
				}

				int64_t interv_range = iinterv->range();
				kid_range += interv_range;
				range -= interv_range;
				((GIntervals *)m_kids_intervals1d.back())->push_back(*iinterv);
			}
		} else {
			set<int> chromids_mask;
			GIntervals all_genome;
			get_all_genome_intervs(all_genome);

			for (GIntervals::const_iterator iinterv = all_genome.begin(); iinterv != all_genome.end(); ++iinterv) {
				int64_t chrom_range = scope1d->range(iinterv->chromid);

				if (!chrom_range) 
					continue;

				if (kid_range && num_avail_kids) {
					// should we allocate a new process?
					if ((uint64_t)kid_range > get_min_scope4process() && (uint64_t)range > get_min_scope4process() && kid_range + chrom_range >= range / num_avail_kids) {
						--num_avail_kids;
						kid_range = 0;
						m_kids_intervals1d.push_back(scope1d->create_masked_copy(chromids_mask));
						chromids_mask.clear();
					}

					--num_remaining_chroms;
					num_avail_kids = min(num_avail_kids, num_remaining_chroms);
				}

				kid_range += chrom_range;
				range -= chrom_range;
				chromids_mask.insert(iinterv->chromid);
			}

			if (!chromids_mask.empty()) 
				m_kids_intervals1d.push_back(scope1d->create_masked_copy(chromids_mask));
		}

//int idx = 0;
//for (vector<GIntervalsFetcher1D *>::iterator i = m_kids_intervals1d.begin(); i != m_kids_intervals1d.end(); ++i){
//GIntervalsFetcher1D *fetcher = *i;
//REprintf("Kid %d:\n", ++idx);
//GIntervals all_genome;
//get_all_genome_intervs(all_genome);
//
//int64_t total_range = 0;
//for (GIntervals::const_iterator iinterv = all_genome.begin(); iinterv != all_genome.end(); ++iinterv) {
//int64_t chrom_range = fetcher->range(iinterv->chromid);
//if (chrom_range) {
//REprintf("\t%s\trange %ld\n", id2chrom(iinterv->chromid).c_str(), chrom_range);
//total_range += chrom_range;
//}
//}
//REprintf("total range: %ld\n\n", total_range);
//}

	} else if (scope2d && scope2d->size()) {
		double surface = scope2d->surface();
		double kid_surface = 0;
		int num_avail_kids = min(max_num_pids, num_chroms) - 1;
		int num_remaining_chroms = num_chroms - 1;

		// GIntervals2D::surface(chromid1, chromid2) and GIntervals2D::create_masked_copy() have complexity of O(n) (and only O(1) for GIntervalsBigSet2D).
		// Calling these functions for each chromosome might be painful - O(m^2 x n) (m = num chromosomes)
		// So we are forced to write a designated code for GIntervals2D for the sake of efficiency. The total complexity for GIntervals is O(n).
		if (dynamic_cast<GIntervals2D *>(scope2d)) {
			GIntervals2D *intervals = (GIntervals2D *)scope2d;

			m_kids_intervals2d.push_back(new GIntervals2D());

			for (GIntervals2D::const_iterator iinterv = intervals->begin(); iinterv != intervals->end(); ++iinterv) {
				if (iinterv != intervals->begin() && !iinterv->is_same_chrom(*(iinterv - 1)) && num_avail_kids) {
					double chrom_surface = 0;
					for (GIntervals2D::const_iterator iinterv2 = iinterv; iinterv2 != intervals->end(); ++iinterv2) {
						if (!iinterv->is_same_chrom(*iinterv2)) 
							break;
						chrom_surface += iinterv2->surface();
					}

					// should we allocate a new process?
					if (kid_surface > get_min_scope4process() && surface > get_min_scope4process() && kid_surface + chrom_surface >= surface / num_avail_kids) {
						--num_avail_kids;
						kid_surface = 0;
						m_kids_intervals2d.push_back(new GIntervals2D());
					}

					--num_remaining_chroms;
					num_avail_kids = min(num_avail_kids, num_remaining_chroms);
				}

				double interv_surface = iinterv->surface();
				kid_surface += interv_surface;
				surface -= interv_surface;
				((GIntervals2D *)m_kids_intervals2d.back())->push_back(*iinterv);
			}
		} else {
			set<ChromPair> chrompairs_mask;
			GIntervals2D all_genome;
			get_all_genome_intervs(all_genome);

			for (GIntervals2D::const_iterator iinterv = all_genome.begin(); iinterv != all_genome.end(); ++iinterv) {
				double chrom_surface = scope2d->surface(iinterv->chromid1(), iinterv->chromid2());

				if (!chrom_surface) 
					continue;

				if (kid_surface && num_avail_kids) {
					// should we allocate a new process?
					if (kid_surface > get_min_scope4process() && kid_surface > get_min_scope4process() && kid_surface + chrom_surface >= surface / num_avail_kids) {
						--num_avail_kids;
						kid_surface = 0;
						m_kids_intervals2d.push_back(scope2d->create_masked_copy(chrompairs_mask));
						chrompairs_mask.clear();
					}

					--num_remaining_chroms;
					num_avail_kids = min(num_avail_kids, num_remaining_chroms);
				}

				kid_surface += chrom_surface;
				surface -= chrom_surface;
				chrompairs_mask.insert(ChromPair(iinterv->chromid1(), iinterv->chromid2()));
			}

			if (!chrompairs_mask.empty())
				m_kids_intervals2d.push_back(scope2d->create_masked_copy(chrompairs_mask));
		}

//int idx = 0;
//for (vector<GIntervalsFetcher2D *>::iterator i = m_kids_intervals2d.begin(); i != m_kids_intervals2d.end(); ++i){
//GIntervalsFetcher2D *fetcher = *i;
//REprintf("Kid %d:\n", ++idx);
//GIntervals2D all_genome;
//get_all_genome_intervs(all_genome);
//
//double total_surface = 0;
//for (GIntervals2D::const_iterator iinterv = all_genome.begin(); iinterv != all_genome.end(); ++iinterv) {
//double chrom_surface = fetcher->surface(iinterv->chromid1(), iinterv->chromid2());
//if (chrom_surface) {
//total_surface += chrom_surface;
//REprintf("\t%s\t%s\tsurface %g\n", id2chrom(iinterv->chromid1()).c_str(), id2chrom(iinterv->chromid2()).c_str(), chrom_surface);
//}
//}
//REprintf("total_surface: %g\n\n", total_surface);
//}
	}

	m_num_planned_kids = max(m_kids_intervals1d.size(), m_kids_intervals2d.size());
	return m_num_planned_kids;
}

static uint64_t estimate_fixed_bin_bins(GIntervalsFetcher1D *intervals1d, int64_t binsize)
{
	if (!intervals1d || binsize <= 0)
		return 0;

	uint64_t bins = 0;
	intervals1d->begin_iter();
	for (; !intervals1d->isend(); intervals1d->next()) {
		const GInterval &interval = intervals1d->cur_interval();
		int64_t start_bin = (int64_t)floor(interval.start / (double)binsize);
		int64_t end_bin = (int64_t)ceil(interval.end / (double)binsize);
		if (end_bin > start_bin)
			bins += (uint64_t)(end_bin - start_bin);
	}

	return bins;
}

static uint64_t estimate_fixed_rect_bins(GIntervalsFetcher2D *intervals2d, int64_t width, int64_t height)
{
	if (!intervals2d || width <= 0 || height <= 0)
		return 0;

	uint64_t bins = 0;
	intervals2d->begin_iter();
	for (; !intervals2d->isend(); intervals2d->next()) {
		const GInterval2D &interval = intervals2d->cur_interval();
		int64_t start_x = (int64_t)floor(interval.start1() / (double)width);
		int64_t end_x = (int64_t)ceil(interval.end1() / (double)width);
		int64_t start_y = (int64_t)floor(interval.start2() / (double)height);
		int64_t end_y = (int64_t)ceil(interval.end2() / (double)height);
		if (end_x > start_x && end_y > start_y) {
			uint64_t xbins = (uint64_t)(end_x - start_x);
			uint64_t ybins = (uint64_t)(end_y - start_y);
			if (xbins && ybins) {
				uint64_t add = xbins * ybins;
				if (add / xbins != ybins || bins > numeric_limits<uint64_t>::max() - add)
					return numeric_limits<uint64_t>::max();
				bins += add;
			}
		}
	}

	return bins;
}

uint64_t IntervUtils::estimate_num_bins(SEXP iterator_policy, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d) const
{
	if (!iterator_policy || Rf_isNull(iterator_policy))
		return 0;

	if ((Rf_isReal(iterator_policy) || Rf_isInteger(iterator_policy)) && Rf_length(iterator_policy) == 1) {
		int64_t binsize = Rf_isReal(iterator_policy) ? (int64_t)REAL(iterator_policy)[0] : INTEGER(iterator_policy)[0];
		return estimate_fixed_bin_bins(scope1d, binsize);
	}

	if ((Rf_isReal(iterator_policy) || Rf_isInteger(iterator_policy)) && Rf_length(iterator_policy) == 2) {
		int64_t width = Rf_isReal(iterator_policy) ? (int64_t)REAL(iterator_policy)[0] : INTEGER(iterator_policy)[0];
		int64_t height = Rf_isReal(iterator_policy) ? (int64_t)REAL(iterator_policy)[1] : INTEGER(iterator_policy)[1];
		return estimate_fixed_rect_bins(scope2d, width, height);
	}

	return 0;
}

// Original 2-parameter version for backward compatibility
bool IntervUtils::distribute_task(uint64_t res_const_size,    // data size in bytes for all the result
								  uint64_t res_record_size)   // size in bytes per datum in the result
{
	return distribute_task(res_const_size, res_record_size, rdb::MT_MODE_MMAP);
}

// New 3-parameter version with explicit mode
bool IntervUtils::distribute_task(uint64_t res_const_size,    // data size in bytes for all the result
								  uint64_t res_record_size,   // size in bytes per datum in the result
								  rdb::MultitaskingMode mode)
{
	return distribute_task(res_const_size, res_record_size, mode, 0);
}

bool IntervUtils::distribute_task(uint64_t res_const_size,    // data size in bytes for all the result
								  uint64_t res_record_size,   // size in bytes per datum in the result
								  rdb::MultitaskingMode mode,
								  uint64_t max_records)
{
	// For now, we only support MMAP mode in this function
	// MT_MODE_SINGLE should be handled by the caller (not calling distribute_task at all)
	if (mode != rdb::MT_MODE_MMAP)
		verror("distribute_task called with unsupported mode %d", mode);

	uint64_t max_records_limit = get_max_data_size();
	if (max_records > 0 && max_records < max_records_limit)
		max_records_limit = max_records;

	uint64_t max_res_size = max_records_limit * res_record_size + m_num_planned_kids * res_const_size;

	rdb::prepare4multitasking(res_const_size, res_record_size, max_res_size, get_max_mem_usage(), m_num_planned_kids);

	for (vector<GIntervalsFetcher1D *>::const_iterator ikid_intervs = m_kids_intervals1d.begin(); ikid_intervs != m_kids_intervals1d.end(); ++ikid_intervs) {
		if (!launch_process())
			return true;
	}

	for (vector<GIntervalsFetcher2D *>::const_iterator ikid_intervs = m_kids_intervals2d.begin(); ikid_intervs != m_kids_intervals2d.end(); ++ikid_intervs) {
		if (!launch_process())
			return true;
	}

	wait_for_kids(*this);
	return false;
}

void ChainIntervals::verify_no_src_overlaps(const GenomeChromKey &chromkey, const vector<string> &src_id2chrom) const
{
	for (const_iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (ChainInterval::SrcCompare()(*iinterv, *(iinterv - 1)))
			TGLError<ChainIntervals>(UNSORTED_INTERVALS, "To verify overlaps chain intervals must be sorted by source");

		if (iinterv->chromid_src == (iinterv - 1)->chromid_src && (iinterv - 1)->end_src > iinterv->start_src)
			TGLError<ChainIntervals>(OVERLAPPING_INTERVAL, "Source of chain intervals %s and %s overlap",
					(iinterv - 1)->tostring(chromkey, src_id2chrom).c_str(), iinterv->tostring(chromkey, src_id2chrom).c_str());
	}
}

void ChainIntervals::verify_no_tgt_overlaps(const GenomeChromKey &chromkey, const vector<string> &src_id2chrom) const
{
	for (const_iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (*iinterv < *(iinterv - 1))
			TGLError<ChainIntervals>(UNSORTED_INTERVALS, "To verify overlaps chain intervals must be sorted by target");

		if (iinterv->chromid == (iinterv - 1)->chromid && (iinterv - 1)->end > iinterv->start)
			TGLError<ChainIntervals>(OVERLAPPING_INTERVAL, "Target of chain intervals %s and %s overlap",
					(iinterv - 1)->tostring(chromkey, src_id2chrom).c_str(), iinterv->tostring(chromkey, src_id2chrom).c_str());
	}
}

void ChainIntervals::handle_src_overlaps(const string &policy, const GenomeChromKey &chromkey, const vector<string> &src_id2chrom)
{
	if (policy == "error") {
		verify_no_src_overlaps(chromkey, src_id2chrom);
	} else if (policy == "keep") {
		// Do nothing - allow overlapping source intervals
	} else if (policy == "discard") {
		// Mark all intervals that have overlapping sources
		vector<bool> to_discard(size(), false);

		for (size_t i = 1; i < size(); ++i) {
			if ((*this)[i-1].chromid_src == (*this)[i].chromid_src &&
			    (*this)[i-1].end_src > (*this)[i].start_src) {
				to_discard[i-1] = true;
				to_discard[i] = true;
			}
		}

		// Remove marked intervals
		iterator new_end = begin();
		for (size_t i = 0; i < size(); ++i) {
			if (!to_discard[i]) {
				if (new_end != begin() + i)
					*new_end = (*this)[i];
				++new_end;
			}
		}
		erase(new_end, end());
	} else {
		TGLError("Invalid source overlap policy: %s. Must be 'error', 'keep', or 'discard'", policy.c_str());
	}
}

void ChainIntervals::handle_tgt_overlaps(const string &policy, const GenomeChromKey &chromkey, const vector<string> &src_id2chrom)
{
	string effective_policy = policy;
	if (effective_policy == "auto")
		effective_policy = "auto_score";

	if (effective_policy == "error") {
		verify_no_tgt_overlaps(chromkey, src_id2chrom);
		return;
	}

	if (empty() || effective_policy == "keep")
		return;

	struct Event {
		int64_t pos;
		bool is_start;
		size_t idx;
		bool operator<(const Event &other) const {
			if (pos != other.pos)
				return pos < other.pos;
			if (is_start != other.is_start)
				return is_start < other.is_start; // closes before opens at same coordinate
			return idx < other.idx;
		}
	};

	auto append_slice = [](vector<ChainInterval> &out,
	                       const ChainInterval &interval,
	                       int64_t seg_start,
	                       int64_t seg_end,
	                       bool allow_merge) {
		if (seg_end <= seg_start)
			return;

		ChainInterval slice = interval;
		const int64_t delta = seg_start - interval.start;
		slice.start = seg_start;
		slice.end = seg_end;
		slice.start_src = interval.start_src + delta;
		slice.end_src = slice.start_src + (seg_end - seg_start);

		if (allow_merge && !out.empty()) {
			ChainInterval &prev = out.back();
			if (prev.chromid == slice.chromid &&
			    prev.chain_id == slice.chain_id &&
			    prev.strand == slice.strand &&
			    prev.chromid_src == slice.chromid_src &&
			    prev.strand_src == slice.strand_src &&
			    prev.end == slice.start &&
			    prev.end_src == slice.start_src) {
				prev.end = slice.end;
				prev.end_src = slice.end_src;
				return;
			}
		}

		out.push_back(slice);
	};

	if (effective_policy == "discard") {
		vector<char> to_discard(size(), 0);

		size_t idx = 0;
		while (idx < size()) {
			const int chromid = (*this)[idx].chromid;
			size_t chrom_end = idx;
			while (chrom_end < size() && (*this)[chrom_end].chromid == chromid)
				++chrom_end;

			vector<Event> events;
			events.reserve((chrom_end - idx) * 2);
			for (size_t k = idx; k < chrom_end; ++k) {
				events.push_back({ (*this)[k].start, true, k });
				events.push_back({ (*this)[k].end, false, k });
			}
			sort(events.begin(), events.end());

			std::set<size_t> active;
			size_t e = 0;
			while (e < events.size()) {
				int64_t pos = events[e].pos;
				while (e < events.size() && events[e].pos == pos && !events[e].is_start) {
					active.erase(events[e].idx);
					++e;
				}
				while (e < events.size() && events[e].pos == pos && events[e].is_start) {
					active.insert(events[e].idx);
					++e;
				}
				if (e >= events.size())
					break;

				int64_t next_pos = events[e].pos;
				if (next_pos <= pos || active.size() <= 1)
					continue;

				for (size_t active_idx : active)
					to_discard[active_idx] = 1;
			}

			idx = chrom_end;
		}

		iterator new_end = begin();
		for (size_t i = 0; i < size(); ++i) {
			if (to_discard[i])
				continue;
			if (new_end != begin() + i)
				*new_end = (*this)[i];
			++new_end;
		}
		erase(new_end, end());
		return;
	}

	auto pick_by_score = [this](const std::set<size_t> &active) {
		size_t best = *active.begin();
		for (size_t idx : active) {
			const ChainInterval &candidate = (*this)[idx];
			const ChainInterval &current = (*this)[best];
			if (candidate.score > current.score)
				best = idx;
			else if (candidate.score == current.score) {
				const int64_t cand_span = candidate.end - candidate.start;
				const int64_t curr_span = current.end - current.start;
				if (cand_span > curr_span)
					best = idx;
				else if (cand_span == curr_span && candidate.chain_id < current.chain_id)
					best = idx;
			}
		}
		return best;
	};

	auto pick_by_length = [this](const std::set<size_t> &active) {
		size_t best = *active.begin();
		for (size_t idx : active) {
			const ChainInterval &candidate = (*this)[idx];
			const ChainInterval &current = (*this)[best];
			const int64_t cand_span = candidate.end - candidate.start;
			const int64_t curr_span = current.end - current.start;
			if (cand_span > curr_span)
				best = idx;
			else if (cand_span == curr_span) {
				if (candidate.score > current.score)
					best = idx;
				else if (candidate.score == current.score && candidate.chain_id < current.chain_id)
					best = idx;
			}
		}
		return best;
	};

	auto pick_first = [this](const std::set<size_t> &active) {
		size_t best = *active.begin();
		for (size_t idx : active) {
			if ((*this)[idx].chain_id < (*this)[best].chain_id)
				best = idx;
		}
		return best;
	};

	const bool is_auto_policy = effective_policy == "auto_first" ||
	                            effective_policy == "auto_longer" ||
	                            effective_policy == "auto_score";

	if (!is_auto_policy && effective_policy != "agg")
		TGLError("Invalid target overlap policy: %s. Must be 'error', 'keep', 'auto_first', 'auto_longer', 'auto_score', 'agg', or 'discard'", effective_policy.c_str());

	vector<ChainInterval> resolved;
	resolved.reserve(size() * 2);

	size_t idx = 0;
	while (idx < size()) {
		const int chromid = (*this)[idx].chromid;
		size_t chrom_end = idx;
		while (chrom_end < size() && (*this)[chrom_end].chromid == chromid)
			++chrom_end;

		vector<Event> events;
		events.reserve((chrom_end - idx) * 2);
		for (size_t k = idx; k < chrom_end; ++k) {
			events.push_back({ (*this)[k].start, true, k });
			events.push_back({ (*this)[k].end, false, k });
		}
		sort(events.begin(), events.end());

		std::set<size_t> active;
		size_t e = 0;
		while (e < events.size()) {
			int64_t pos = events[e].pos;
			while (e < events.size() && events[e].pos == pos && !events[e].is_start) {
				active.erase(events[e].idx);
				++e;
			}
			while (e < events.size() && events[e].pos == pos && events[e].is_start) {
				active.insert(events[e].idx);
				++e;
			}
			if (e >= events.size())
				break;

			int64_t next_pos = events[e].pos;
			if (next_pos <= pos || active.empty())
				continue;

			if (active.size() == 1) {
				bool allow_merge = effective_policy != "agg";
				append_slice(resolved, (*this)[*active.begin()], pos, next_pos, allow_merge);
				continue;
			}

			if (effective_policy == "agg") {
				for (size_t active_idx : active)
					append_slice(resolved, (*this)[active_idx], pos, next_pos, false);
				continue;
			}

			size_t winner = 0;
			if (effective_policy == "auto_score")
				winner = pick_by_score(active);
			else if (effective_policy == "auto_longer")
				winner = pick_by_length(active);
			else
				winner = pick_first(active);

			append_slice(resolved, (*this)[winner], pos, next_pos, true);
		}

		idx = chrom_end;
	}

	assign(resolved.begin(), resolved.end());
}

void ChainIntervals::buildSrcAux()
{
	const size_t n = size();
	m_pmax_end_src.assign(n, 0);

	if (n == 0)
		return;

	// Build chrom slices: first/last indices per chromid_src
	// chromid_src is a small dense int assigned by convert_rchain_intervs()
	// Track max id to size vectors
	m_max_src_chromid = -1;
	for (const auto &ci : *this)
		m_max_src_chromid = std::max(m_max_src_chromid, ci.chromid_src);

	// Initialize with sentinel values
	m_chrom_first.assign((size_t)m_max_src_chromid + 1, (size_t)-1);
	m_chrom_last_excl.assign((size_t)m_max_src_chromid + 1, (size_t)-1);

	// Fill first/last indices and compute prefix-max per chromosome
	size_t i = 0;
	while (i < n) {
		const int chrom = (*this)[i].chromid_src;
		size_t first = i;
		if (m_chrom_first[chrom] == (size_t)-1)
			m_chrom_first[chrom] = first;

		// Walk this chromosome slice and build prefix-max
		int64_t pmax = std::numeric_limits<int64_t>::min();
		for (; i < n && (*this)[i].chromid_src == chrom; ++i) {
			const auto &ci = (*this)[i];
			const int64_t end_src = ci.end_src;
			pmax = std::max(pmax, end_src);
			m_pmax_end_src[i] = pmax;
		}
		m_chrom_last_excl[chrom] = i; // exclusive end
	}
}

ChainIntervals::const_iterator ChainIntervals::map_interval(const GInterval &src_interval, GIntervals &tgt_intervs, ChainIntervals::const_iterator hint, std::vector<ChainMappingMetadata> *metadata)
{
	tgt_intervs.clear();
	if (metadata)
		metadata->clear();

	if (empty())
		return end();

	if (front().chromid_src > src_interval.chromid || (front().chromid_src == src_interval.chromid && front().start_src >= src_interval.end))
		return begin();

	if (back().chromid_src < src_interval.chromid || (back().chromid_src == src_interval.chromid && back().start_src + back().end - back().start <= src_interval.start))
		return end() - 1;

	if (check_first_overlap_src(hint, src_interval))
		return add2tgt(hint, src_interval, tgt_intervs, metadata);

	// Note: hint+1 fast path removed to handle non-consecutive overlapping chains correctly
	// With overlapping sources, hint+1 might overlap but not be the *first* overlap

	// run the binary search
	const_iterator istart_interval = begin();
	const_iterator iend_interval = end();

	while (iend_interval - istart_interval > 1) {
		const_iterator imid_interval = istart_interval + (iend_interval - istart_interval) / 2;

		// If mid overlaps, use prefix-max to find the first overlapping chain
		if (imid_interval->do_overlap_src(src_interval)) {
			// Use prefix-max to find the true first overlap, not just consecutive first
			const int qchrom = src_interval.chromid;
			if (qchrom < 0 || qchrom > m_max_src_chromid || m_chrom_first[qchrom] == (size_t)-1) {
				// Shouldn't happen due to guards above, but fallback to consecutive scan
				const_iterator first_overlap = imid_interval;
				while (first_overlap > begin() && (first_overlap - 1)->do_overlap_src(src_interval))
					--first_overlap;
				return add2tgt(first_overlap, src_interval, tgt_intervs, metadata);
			}

			size_t first = m_chrom_first[qchrom];
			size_t mid_idx = imid_interval - begin();

			// Binary search leftmost L in [first..mid_idx] with pmax_end[L] > query.start
			// We know mid overlaps, so pmax_end[mid_idx] > query.start
			size_t lo = first, hi = mid_idx;
			while (lo < hi) {
				size_t mid = lo + (hi - lo) / 2;
				if (m_pmax_end_src[mid] > src_interval.start)
					hi = mid;
				else
					lo = mid + 1;
			}
			return add2tgt(begin() + lo, src_interval, tgt_intervs, metadata);
		}

		// is mid_interval < interval?
		if (imid_interval->chromid_src < src_interval.chromid || (imid_interval->chromid_src == src_interval.chromid && imid_interval->start_src < src_interval.start))
			istart_interval = imid_interval;
		else
			iend_interval = imid_interval;
	}

	// After binary search, use prefix-max to find the first overlapping chain efficiently
	// Get chromosome slice bounds for this source chromosome
	const int qchrom = src_interval.chromid;

	// Guard: if chromosome doesn't exist in our slices (shouldn't happen due to fast guards above)
	if (qchrom < 0 || qchrom > m_max_src_chromid || m_chrom_first[qchrom] == (size_t)-1) {
		// Fallback to checking iend_interval
		if (iend_interval != end() && iend_interval->do_overlap_src(src_interval))
			return add2tgt(iend_interval, src_interval, tgt_intervs, metadata);
		return istart_interval;
	}

	size_t first = m_chrom_first[qchrom];
	size_t lastEx = m_chrom_last_excl[qchrom];

	// Recompute pos: last index with start_src < query.start within [first, lastEx)
	const_iterator slice_begin = begin() + first;
	const_iterator slice_end = begin() + lastEx;
	auto lb = std::lower_bound(slice_begin, slice_end, src_interval.start,
		[](const ChainInterval &ci, int64_t qstart) {
			return ci.start_src < qstart;
		});

	// lower_bound returns first with start_src >= qstart; we want last < qstart
	ptrdiff_t pos_abs = (lb == slice_begin) ? -1 : (lb - begin() - 1);

	if (pos_abs >= 0) {
		size_t pos = (size_t)pos_abs;
		// Check if there exists an earlier overlap from the left: pmax_end[pos] > query.start
		if (m_pmax_end_src[pos] > src_interval.start) {
			// Binary search leftmost L in [first..pos] with pmax_end[L] > query.start
			size_t lo = first, hi = pos;
			while (lo < hi) {
				size_t mid = lo + (hi - lo) / 2;
				if (m_pmax_end_src[mid] > src_interval.start)
					hi = mid;
				else
					lo = mid + 1;
			}
			// Start enumeration from lo
			const_iterator last_checked = add2tgt(begin() + lo, src_interval, tgt_intervs, metadata);
			// Optimize hint for next query: return pos if it's ahead and in same chrom
			if (pos >= (size_t)(last_checked - begin()) && pos < lastEx)
				return begin() + pos;
			return last_checked;
		}
	}

	// No left overlap exists; check the right neighbor (iend_interval)
	if (iend_interval != end() && iend_interval->do_overlap_src(src_interval))
		return add2tgt(iend_interval, src_interval, tgt_intervs, metadata);

	// Nothing overlaps
	return istart_interval;
}

ChainIntervals::const_iterator ChainIntervals::add2tgt(const_iterator hint, const GInterval &src_interval, GIntervals &tgt_intervs, std::vector<ChainMappingMetadata> *metadata)
{
	const_iterator last_checked = hint;

	// Temporary storage for all mapped intervals with their metadata
	struct MappedInterval {
		GInterval interval;
		double score;
		int64_t chain_id;
		int64_t span;
		const_iterator chain_iter;
		int     chromid_src;
		int64_t start_src_mapped;  // Intersection start (common_start)
		int64_t end_src_mapped;    // Intersection end (common_end)
	};
	std::vector<MappedInterval> all_mapped;

	// Continue until we reach a chain whose start_src >= src_interval.end or different chromid
	// This handles the case where overlapping source regions cause non-consecutive overlapping chains
	while (hint != end() && hint->chromid_src == src_interval.chromid && hint->start_src < src_interval.end) {
		// Quick early exit: if chain ends before query starts, no overlap possible
		if (hint->end_src <= src_interval.start) {
			last_checked = hint;
			++hint;
			continue;
		}

		if (hint->do_overlap_src(src_interval)) {
			int64_t common_start = max(hint->start_src, src_interval.start);
			int64_t common_end = min(hint->end_src, src_interval.end);

			int64_t tgt_start, tgt_end;
			if (hint->strand == 0) {
				// Target is on + strand: offset goes in positive direction
				tgt_start = hint->start + common_start - hint->start_src;
				tgt_end = hint->start + common_end - hint->start_src;
			} else {
				// Target is on - strand: offset goes in negative direction
				tgt_start = hint->end - (common_end - hint->start_src);
				tgt_end = hint->end - (common_start - hint->start_src);
			}

			MappedInterval mapped;
			mapped.interval = GInterval(hint->chromid, tgt_start, tgt_end, hint->strand);
			mapped.score = hint->score;
			mapped.chain_id = hint->chain_id;
			mapped.span = tgt_end - tgt_start;
			mapped.chain_iter = hint;
			mapped.chromid_src = hint->chromid_src;
			mapped.start_src_mapped = common_start;
			mapped.end_src_mapped = common_end;
			all_mapped.push_back(mapped);

			last_checked = hint;
			++hint;
		} else {
			// No overlap, but continue scanning in case there are non-consecutive overlaps
			last_checked = hint;
			++hint;
		}
	}

	// Apply score-based filtering if policy is auto_score
	if (m_tgt_overlap_policy == "auto_score" || m_tgt_overlap_policy == "auto") {
		if (!all_mapped.empty()) {
			// Build list of intervals to keep after score-based filtering
			std::vector<bool> keep(all_mapped.size(), true);

			// For each pair of overlapping target intervals, discard the worse one
			for (size_t i = 0; i < all_mapped.size(); ++i) {
				if (!keep[i]) continue;

				for (size_t j = i + 1; j < all_mapped.size(); ++j) {
					if (!keep[j]) continue;

					// Check if intervals overlap
					const GInterval &int_i = all_mapped[i].interval;
					const GInterval &int_j = all_mapped[j].interval;

					if (int_i.chromid == int_j.chromid && int_i.start < int_j.end && int_j.start < int_i.end) {
						// They overlap - compare and keep the better one
						bool i_wins;

						// Compare by score first
						if (all_mapped[i].score != all_mapped[j].score) {
							i_wins = all_mapped[i].score > all_mapped[j].score;
						}
						// If scores equal, compare by span
						else if (all_mapped[i].span != all_mapped[j].span) {
							i_wins = all_mapped[i].span > all_mapped[j].span;
						}
						// If spans also equal, use chain_id (lower is better for stability)
						else {
							i_wins = all_mapped[i].chain_id < all_mapped[j].chain_id;
						}

						if (i_wins) {
							keep[j] = false;
						} else {
							keep[i] = false;
							break; // i is discarded, no need to compare further
						}
					}
				}
			}

			// Copy filtered intervals to output
			for (size_t i = 0; i < all_mapped.size(); ++i) {
				if (keep[i]) {
					tgt_intervs.push_back(all_mapped[i].interval);
					if (metadata) {
						ChainMappingMetadata meta;
						meta.score = all_mapped[i].score;
						meta.chain_id = all_mapped[i].chain_id;
						meta.chromid_src = all_mapped[i].chromid_src;
						meta.start_src = all_mapped[i].start_src_mapped;
						meta.end_src = all_mapped[i].end_src_mapped;
						metadata->push_back(meta);
					}
				}
			}
		}
	} else {
		// For other policies (keep, etc.), add all intervals
		for (const auto &mapped : all_mapped) {
			tgt_intervs.push_back(mapped.interval);
			if (metadata) {
				ChainMappingMetadata meta;
				meta.score = mapped.score;
				meta.chain_id = mapped.chain_id;
				meta.chromid_src = mapped.chromid_src;
				meta.start_src = mapped.start_src_mapped;
				meta.end_src = mapped.end_src_mapped;
				metadata->push_back(meta);
			}
		}
	}

	return (last_checked != end()) ? last_checked : (hint != begin()) ? hint - 1 : begin();
}
