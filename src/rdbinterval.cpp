/*
 * rdbinterval.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <cmath>
#include <cctype>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

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

const char *ChainInterval::COL_NAMES[ChainInterval::NUM_COLS] = { "chrom", "start", "end", "strand", "chromsrc", "startsrc", "endsrc", "strandsrc" };

IntervUtils::IntervUtils(SEXP envir)
{
	m_envir = envir;
	m_num_planned_kids = 0;
	m_multitasking = -1;
	m_kid_intervals1d = NULL;
	m_kid_intervals2d = NULL;

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

unsigned IntervUtils::get_rintervs_type_mask(SEXP rintervals, const char *error_msg_prefix) const
{
	if (!Rf_isVector(rintervals))
		verror("%sInvalid format of intervals argument", error_msg_prefix);

	if (Rf_length(rintervals) == 2) {
		if (get_rintervs_type_mask(VECTOR_ELT(rintervals, 0), error_msg_prefix) != INTERVS1D || get_rintervs_type_mask(VECTOR_ELT(rintervals, 1), error_msg_prefix) != INTERVS2D)
			verror("%sInvalid format of intervals argument", error_msg_prefix);
		return INTERVS1D | INTERVS2D;
	}

	SEXP colnames = Rf_getAttrib(rintervals, R_NamesSymbol);

	if (!Rf_isString(colnames) || Rf_length(colnames) < GInterval::NUM_COLS)
		verror("%sInvalid format of intervals argument", error_msg_prefix);

	IntervUtils::IntervsType type = INTERVS1D;

	for (unsigned i = 0; i < GInterval::NUM_COLS; i++) {
		if (strcmp(CHAR(STRING_ELT(colnames, i)), GInterval::COL_NAMES[i])) {
			type = INTERVS2D;
			break;
		}
	}

	if (type == INTERVS2D) {
		for (unsigned i = 0; i < GInterval2D::NUM_COLS; i++) {
			if (strcmp(CHAR(STRING_ELT(colnames, i)), GInterval2D::COL_NAMES[i]))
				verror("Invalid format of intervals: column names do not match neither 1d nor 2d intervals");
		}
	}

	if (type == INTERVS1D) {
		SEXP starts = VECTOR_ELT(rintervals, GInterval::START);
		SEXP ends = VECTOR_ELT(rintervals, GInterval::END);
		SEXP strands = R_NilValue;
		SEXP colnames = Rf_getAttrib(rintervals, R_NamesSymbol);

		for (int i = 0; i < Rf_length(rintervals); i++) {
			if (!strcmp(CHAR(STRING_ELT(colnames, i)), "strand")) {
				if (Rf_length(VECTOR_ELT(rintervals, i)) != Rf_length(VECTOR_ELT(rintervals, GInterval::CHROM)))
					verror("%sNumber of rows in column %s differs than the number of rows in column strand", error_msg_prefix, GInterval::COL_NAMES[GInterval::CHROM]);
				break;
			}
		}

		for (unsigned i = 0; i < GInterval::NUM_COLS; i++) {
			if (i != 0 && Rf_length(VECTOR_ELT(rintervals, i)) != Rf_length(VECTOR_ELT(rintervals, i - 1)))
				verror("%sNumber of rows in column %s differs than the number of rows in column %s", error_msg_prefix, GInterval::COL_NAMES[i - 1], GInterval::COL_NAMES[i]);
		}

		if ((!Rf_isReal(starts) && !Rf_isInteger(starts)) || (!Rf_isReal(ends) && !Rf_isInteger(ends)) || (strands != R_NilValue && !Rf_isReal(strands) && !Rf_isInteger(strands)))
			verror("%sInvalid format of intervals argument", error_msg_prefix);

	} else if (type == INTERVS2D) {
		SEXP starts1 = VECTOR_ELT(rintervals, GInterval2D::START1);
		SEXP ends1 = VECTOR_ELT(rintervals, GInterval2D::END1);
		SEXP starts2 = VECTOR_ELT(rintervals, GInterval2D::START2);
		SEXP ends2 = VECTOR_ELT(rintervals, GInterval2D::END2);

		for (unsigned i = 0; i < GInterval2D::NUM_COLS; i++) {
			if (i != 0 && Rf_length(VECTOR_ELT(rintervals, i)) != Rf_length(VECTOR_ELT(rintervals, i - 1)))
				verror("%sNumber of rows in column %s differs than the number of rows in column %s", error_msg_prefix, GInterval2D::COL_NAMES[i - 1], GInterval2D::COL_NAMES[i]);
		}

		if ((!Rf_isReal(starts1) && !Rf_isInteger(starts1)) || (!Rf_isReal(ends1) && !Rf_isInteger(ends1)) || (!Rf_isReal(starts2) && !Rf_isInteger(starts2)) || (!Rf_isReal(ends2) && !Rf_isInteger(ends2)))
			verror("%sInvalid format of intervals argument", error_msg_prefix);
	} else
		verror("%sUnexpected intervals type", error_msg_prefix);

	return type;
}

unsigned IntervUtils::convert_rintervs(SEXP rintervals, GIntervalsFetcher1D **intervals1d, GIntervalsFetcher2D **intervals2d, bool null_if_interv_nonexist,
									   const GenomeChromKey *chromkey, const char *error_msg_prefix, bool verify) const
{
	monitor_memusage();
	if (intervals1d) 
		*intervals1d = NULL;
	if (intervals2d) 
		*intervals2d = NULL;
	try {
		if (Rf_isString(rintervals) && Rf_length(rintervals) == 1) {
			const char *full_interv_name = CHAR(STRING_ELT(rintervals, 0));

			if (GIntervalsBigSet::isbig(full_interv_name, *this)) { // big intervals set?
				SEXP meta = GIntervalsMeta::load_meta(interv2path(get_env(), full_interv_name).c_str());

				if (GIntervalsMeta1D::is1d(meta)) {
					if (intervals1d) {
						*intervals1d = new GIntervalsBigSet1D(full_interv_name, meta, *this);
						if (intervals2d)
							*intervals2d = new GIntervals2D();  // create an empty intervals set for convenience
						return IntervUtils::INTERVS1D;
					}
					verror("%sExpecting 2D intervals while receiving 1D intervals set %s", error_msg_prefix, full_interv_name);
				} else {
					if (intervals2d) {
						*intervals2d = new GIntervalsBigSet2D(full_interv_name, meta, *this);
						if (intervals1d)
							*intervals1d = new GIntervals();  // create an empty intervals set for convenience
						return IntervUtils::INTERVS2D;
					}
					verror("%sExpecting 1D intervals while receiving 2D intervals set %s", error_msg_prefix, full_interv_name);
				}
			} else if (GTrackIntervalsFetcher::isbig(full_interv_name, *this)) {
				string trackpath = track2path(get_env(), full_interv_name);

				if (access((trackpath + "/.meta").c_str(), R_OK) && errno == ENOENT)
					GTrackIntervalsFetcher::create_track_meta(full_interv_name, *this);

				SEXP meta = GIntervalsMeta::load_meta(trackpath.c_str());

				if (GIntervalsMeta1D::is1d(meta)) {
					if (intervals1d) {
						GenomeTrack::Type track_type = GenomeTrack::get_type(trackpath.c_str(), get_chromkey(), true);

						if (track_type == GenomeTrack::FIXED_BIN)
							verror("%s is track of %s type which cannot be used as a replacement for intervals",
								   full_interv_name, GenomeTrack::TYPE_NAMES[track_type]);
						else if (track_type == GenomeTrack::SPARSE)
							*intervals1d = new GTrackIntervalsFetcher1D<GenomeTrackSparse>(full_interv_name, meta, *this);
						else if (track_type == GenomeTrack::ARRAYS)
							*intervals1d = new GTrackIntervalsFetcher1D<GenomeTrackArrays>(full_interv_name, meta, *this);

						if (intervals2d)
							*intervals2d = new GIntervals2D();  // create an empty intervals set for convenience
						return IntervUtils::INTERVS1D;
					}
					verror("%sExpecting 2D intervals while receiving 1D intervals set %s", error_msg_prefix, full_interv_name);
				} else {
					if (intervals2d) {
						GenomeTrack::Type track_type = GenomeTrack::get_type(trackpath.c_str(), get_chromkey(), true);

						if (track_type == GenomeTrack::RECTS)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackRectsRects>(full_interv_name, meta, *this);
						else if (track_type == GenomeTrack::POINTS)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackRectsPoints>(full_interv_name, meta, *this);
						else if (track_type == GenomeTrack::COMPUTED)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackComputed>(full_interv_name, meta, *this);

						if (intervals1d)
							*intervals1d = new GIntervals();  // create an empty intervals set for convenience
						return IntervUtils::INTERVS2D;
					}
					verror("%sExpecting 1D intervals while receiving 2D intervals set %s", error_msg_prefix, full_interv_name);
				}
			}
		}

		if (intervals1d)
			*intervals1d = new GIntervals();
		if (intervals2d)
			*intervals2d = new GIntervals2D();

		unsigned intervs_type_mask;
		convert_rintervs(rintervals, intervals1d ? (GIntervals *)*intervals1d : NULL, intervals2d ? (GIntervals2D *)*intervals2d : NULL,
						 null_if_interv_nonexist, chromkey, error_msg_prefix, &intervs_type_mask, verify);
		return intervs_type_mask;
	} catch (...) {
		if (intervals1d) {
			delete *intervals1d;
			*intervals1d = NULL;
		}
		if (intervals2d) {
			delete *intervals2d;
			*intervals2d = NULL;
		}
		throw;
	}
}

SEXP IntervUtils::convert_rintervs(SEXP rintervals, GIntervals *intervals, GIntervals2D *intervals2d, bool null_if_interv_nonexist,
								   const GenomeChromKey *chromkey, const char *error_msg_prefix, unsigned *pintervs_type_mask, bool verify) const
{
	monitor_memusage();

	if (intervals)
		intervals->clear();

	if (intervals2d)
		intervals2d->clear();

	if (pintervs_type_mask) 
		*pintervs_type_mask = 0;

	bool loaded_from_file = Rf_isString(rintervals);

	// rintervals is the name of the intervals file
	if (Rf_isString(rintervals) && Rf_length(rintervals) == 1) {
		const char *full_interv_name = CHAR(STRING_ELT(rintervals, 0));
		SEXP gintervs;
		bool interv_found = false;

        {
            rprotect(gintervs = find_in_misha(m_envir, "GINTERVS"));
        }
		if (Rf_isString(gintervs)) {
			for (int iinterv = 0; iinterv < Rf_length(gintervs); ++iinterv) {
				const char *interv = CHAR(STRING_ELT(gintervs, iinterv));
				if (!strcmp(full_interv_name, interv)) {
					interv_found = true;
					break;
				}
			}
		}

		if (!interv_found) {
			if (GTrackIntervalsFetcher::isbig(full_interv_name, *this))
				verror("Tracks cannot be used as a substitute for intervals in this function");

			if (null_if_interv_nonexist)
				return R_NilValue;

			verror("%sInterval %s does not exist", error_msg_prefix, full_interv_name);
		}

		string path = interv2path(m_envir, full_interv_name);
		struct stat stat_res;
		if (!stat(path.c_str(), &stat_res) && S_ISDIR(stat_res.st_mode))
			verror("%s is a big intervals set. Big intervals sets are not supported by the function.", full_interv_name);
		rprotect(rintervals = RSaneUnserialize(path.c_str()));
        }

	if (TYPEOF(rintervals) == PROMSXP) {
		if (MISHA_PRENV(rintervals) == R_NilValue)
			rintervals = MISHA_PRVALUE(rintervals);
		else
			rintervals = eval_in_R(MISHA_PREXPR(rintervals), MISHA_PRENV(rintervals));
	}

	SEXP _rintervals = rintervals;
	unsigned intervs_type_mask = get_rintervs_type_mask(rintervals, error_msg_prefix);

	if (pintervs_type_mask) 
		*pintervs_type_mask = intervs_type_mask;

	if (intervs_type_mask == (INTERVS1D | INTERVS2D))
		rintervals = VECTOR_ELT(_rintervals, 0);

	if (intervs_type_mask == INTERVS1D && !intervals)
		verror("%sExpecting 2D intervals while receiving 1D intervals", error_msg_prefix);

	if (intervs_type_mask == INTERVS2D && !intervals2d)
		verror("%sExpecting 1D intervals while receiving 2D intervals", error_msg_prefix);

	if ((intervs_type_mask & INTERVS1D) && intervals) {
		SEXP chroms = VECTOR_ELT(rintervals, GInterval::CHROM);
		SEXP chrom_levels = Rf_getAttrib(chroms, R_LevelsSymbol);
		SEXP starts = VECTOR_ELT(rintervals, GInterval::START);
		SEXP ends = VECTOR_ELT(rintervals, GInterval::END);
		SEXP strands = R_NilValue;
		SEXP colnames = Rf_getAttrib(rintervals, R_NamesSymbol);
		unsigned num_intervals = (unsigned)Rf_length(starts);

		for (int i = 0; i < Rf_length(rintervals); i++) {
			if (!strcmp(CHAR(STRING_ELT(colnames, i)), "strand")) {
				strands = VECTOR_ELT(rintervals, i);
				break;
			}
		}

		for (unsigned i = 0; i < num_intervals; i++) {
			if ((Rf_isFactor(chroms) && INTEGER(chroms)[i] < 0) ||
				(Rf_isReal(starts) && std::isnan(REAL(starts)[i])) || (Rf_isReal(ends) && std::isnan(REAL(ends)[i])) ||
				(strands != R_NilValue && Rf_isReal(strands) && std::isnan(REAL(strands)[i])))
				verror("%sInvalid format of interval at index %d", error_msg_prefix, i + 1);

			const char *chrom = Rf_isString(chroms) ? CHAR(STRING_ELT(chroms, i)) : CHAR(STRING_ELT(chrom_levels, INTEGER(chroms)[i] - 1));
			int chromid = chromkey ? chromkey->chrom2id(chrom) : chrom2id(chrom);
			int64_t start = (int64_t)(Rf_isReal(starts) ? REAL(starts)[i] : INTEGER(starts)[i]);
			int64_t end = (int64_t)(Rf_isReal(ends) ? REAL(ends)[i] : INTEGER(ends)[i]);
			char strand = 0;

			if (strands != R_NilValue)
				strand = (char)(Rf_isReal(strands) ? REAL(strands)[i] : INTEGER(strands)[i]);

			GInterval interval(chromid, start, end, strand, (void *)(intptr_t)i);

			if (loaded_from_file && interval.start == interval.end) {
				interval.end++;
				if (Rf_isReal(ends))
					REAL(ends)[i]++;
				else
					INTEGER(ends)[i]++;
			}

			if (verify) 
				interval.verify(chromkey ? *chromkey : m_chrom_key);
			intervals->push_back(interval);
		}
	}

	if (intervs_type_mask == (INTERVS1D | INTERVS2D))
		rintervals = VECTOR_ELT(_rintervals, 1);

	if ((intervs_type_mask & INTERVS2D) && intervals2d) {
		SEXP chroms1 = VECTOR_ELT(rintervals, GInterval2D::CHROM1);
		SEXP chrom_levels1 = Rf_getAttrib(chroms1, R_LevelsSymbol);
		SEXP starts1 = VECTOR_ELT(rintervals, GInterval2D::START1);
		SEXP ends1 = VECTOR_ELT(rintervals, GInterval2D::END1);
		SEXP chroms2 = VECTOR_ELT(rintervals, GInterval2D::CHROM2);
		SEXP chrom_levels2 = Rf_getAttrib(chroms2, R_LevelsSymbol);
		SEXP starts2 = VECTOR_ELT(rintervals, GInterval2D::START2);
		SEXP ends2 = VECTOR_ELT(rintervals, GInterval2D::END2);
		unsigned num_intervals = (unsigned)Rf_length(starts1);

		for (unsigned i = 0; i < num_intervals; i++) {
			if ((Rf_isFactor(chroms1) && INTEGER(chroms1)[i] < 0) ||
				(Rf_isReal(starts1) && std::isnan(REAL(starts1)[i])) || (Rf_isReal(ends1) && std::isnan(REAL(ends1)[i])) ||
				(Rf_isFactor(chroms2) && INTEGER(chroms2)[i] < 0) ||
				(Rf_isReal(starts2) && std::isnan(REAL(starts2)[i])) || (Rf_isReal(ends2) && std::isnan(REAL(ends2)[i])))
				verror("%sInvalid format of interval at index %d", error_msg_prefix, i + 1);

			const char *chrom1 = Rf_isString(chroms1) ? CHAR(STRING_ELT(chroms1, i)) : CHAR(STRING_ELT(chrom_levels1, INTEGER(chroms1)[i] - 1));
			const char *chrom2 = Rf_isString(chroms2) ? CHAR(STRING_ELT(chroms2, i)) : CHAR(STRING_ELT(chrom_levels2, INTEGER(chroms2)[i] - 1));
			int chromid1 = chromkey ? chromkey->chrom2id(chrom1) : chrom2id(chrom1);
			int chromid2 = chromkey ? chromkey->chrom2id(chrom2) : chrom2id(chrom2);
			int64_t start1 = (int64_t)(Rf_isReal(starts1) ? REAL(starts1)[i] : INTEGER(starts1)[i]);
			int64_t start2 = (int64_t)(Rf_isReal(starts2) ? REAL(starts2)[i] : INTEGER(starts2)[i]);
			int64_t end1 = (int64_t)(Rf_isReal(ends1) ? REAL(ends1)[i] : INTEGER(ends1)[i]);
			int64_t end2 = (int64_t)(Rf_isReal(ends2) ? REAL(ends2)[i] : INTEGER(ends2)[i]);

			GInterval2D interval(chromid1, start1, end1, chromid2, start2, end2, (void *)(intptr_t)i);

			if (verify) 
				interval.verify(chromkey ? *chromkey : m_chrom_key);
			intervals2d->push_back(interval);
		}
	}

	return _rintervals;
}

SEXP IntervUtils::convert_intervs(GIntervalsFetcher1D *intervals, unsigned num_cols, bool null_if_empty, bool use_original_index) const
{
	monitor_memusage();

	if (null_if_empty && !intervals->size())
		return R_NilValue;

	unsigned num_chroms = m_chrom_key.get_num_chroms();

	SEXP answer;
	SEXP chroms, chroms_idx, starts, ends;
	SEXP row_names;
	SEXP col_names;

	rprotect(answer = RSaneAllocVector(VECSXP, num_cols));
    rprotect(chroms_idx = RSaneAllocVector(INTSXP, intervals->size()));
    rprotect(starts = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(ends = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(chroms = RSaneAllocVector(STRSXP, num_chroms));
    rprotect(col_names = RSaneAllocVector(STRSXP, num_cols));
    rprotect(row_names = RSaneAllocVector(INTSXP, intervals->size()));

	for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
		const GInterval &interval = intervals->cur_interval();
		uint64_t index = use_original_index ? get_orig_interv_idx(interval) : intervals->iter_index();

		INTEGER(chroms_idx)[index] = interval.chromid + 1;
		REAL(starts)[index] = interval.start;
		REAL(ends)[index] = interval.end;
		INTEGER(row_names)[index] = index + 1;
	}

	for (unsigned id = 0; id < (unsigned)num_chroms; ++id)
		SET_STRING_ELT(chroms, id, Rf_mkChar(m_chrom_key.id2chrom(id).c_str()));

	for (int i = 0; i < GInterval::NUM_COLS; i++)
		SET_STRING_ELT(col_names, i, Rf_mkChar(GInterval::COL_NAMES[i]));

    Rf_setAttrib(chroms_idx, R_LevelsSymbol, chroms);
    Rf_setAttrib(chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(answer, GInterval::CHROM, chroms_idx);
    SET_VECTOR_ELT(answer, GInterval::START, starts);
    SET_VECTOR_ELT(answer, GInterval::END, ends);

    Rf_setAttrib(answer, R_NamesSymbol, col_names);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, row_names);

	return answer;
}

SEXP IntervUtils::convert_intervs(GIntervalsFetcher2D *intervals, unsigned num_cols, bool null_if_empty, bool use_original_index) const
{
	monitor_memusage();

	if (null_if_empty && !intervals->size())
		return R_NilValue;

	unsigned num_chroms = m_chrom_key.get_num_chroms();

	SEXP answer;
	SEXP chroms1, chroms2, chroms_idx1, chroms_idx2, starts1, starts2, ends1, ends2;
	SEXP row_names;
	SEXP col_names;

	rprotect(answer = RSaneAllocVector(VECSXP, num_cols));
    rprotect(chroms1 = RSaneAllocVector(STRSXP, num_chroms));
    rprotect(starts1 = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(ends1 = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(chroms_idx1 = RSaneAllocVector(INTSXP, intervals->size()));
    rprotect(chroms_idx2 = RSaneAllocVector(INTSXP, intervals->size()));
    rprotect(starts2 = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(ends2 = RSaneAllocVector(REALSXP, intervals->size()));
    rprotect(chroms1 = RSaneAllocVector(STRSXP, num_chroms));
    rprotect(chroms2 = RSaneAllocVector(STRSXP, num_chroms));
    rprotect(col_names = RSaneAllocVector(STRSXP, num_cols));
    rprotect(row_names = RSaneAllocVector(INTSXP, intervals->size()));

	for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
		const GInterval2D &interval = intervals->cur_interval();
		uint64_t index = use_original_index ? get_orig_interv_idx(interval) : intervals->iter_index();

		INTEGER(chroms_idx1)[index] = interval.chromid1() + 1;
		REAL(starts1)[index] = interval.start1();
        REAL(ends1)[index] = interval.end1();
        INTEGER(chroms_idx2)[index] = interval.chromid2() + 1;
        REAL(starts2)[index] = interval.start2();
        REAL(ends2)[index] = interval.end2();
		INTEGER(row_names)[index] = index + 1;
	}

    for (unsigned id = 0; id < (unsigned)num_chroms; ++id) {
        SET_STRING_ELT(chroms1, id, Rf_mkChar(m_chrom_key.id2chrom(id).c_str()));
        SET_STRING_ELT(chroms2, id, Rf_mkChar(m_chrom_key.id2chrom(id).c_str()));
    }

    for (int i = 0; i < GInterval2D::NUM_COLS; i++)
        SET_STRING_ELT(col_names, i, Rf_mkChar(GInterval2D::COL_NAMES[i]));


    Rf_setAttrib(chroms_idx1, R_LevelsSymbol, chroms1);
    Rf_setAttrib(chroms_idx1, R_ClassSymbol, Rf_mkString("factor"));
    Rf_setAttrib(chroms_idx2, R_LevelsSymbol, chroms2);
    Rf_setAttrib(chroms_idx2, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(answer, GInterval2D::CHROM1, chroms_idx1);
    SET_VECTOR_ELT(answer, GInterval2D::START1, starts1);
    SET_VECTOR_ELT(answer, GInterval2D::END1, ends1);
    SET_VECTOR_ELT(answer, GInterval2D::CHROM2, chroms_idx2);
    SET_VECTOR_ELT(answer, GInterval2D::START2, starts2);
    SET_VECTOR_ELT(answer, GInterval2D::END2, ends2);

    Rf_setAttrib(answer, R_NamesSymbol, col_names);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, row_names);

    return answer;
}

void IntervUtils::convert_rchain_intervs(SEXP rchain, ChainIntervals &chain_intervs, vector<string> &src_id2chrom)
{
	if (!Rf_isVector(rchain) || Rf_length(rchain) != ChainInterval::NUM_COLS)
		TGLError("Invalid format of chain argument");

	SEXP colnames = Rf_getAttrib(rchain, R_NamesSymbol);

	if (!Rf_isString(colnames) || Rf_length(colnames) != ChainInterval::NUM_COLS)
		verror("Invalid format of chain argument");

	for (unsigned i = 0; i < ChainInterval::NUM_COLS; i++) {
		if (strcmp(CHAR(STRING_ELT(colnames, i)), ChainInterval::COL_NAMES[i]))
			verror("Invalid format of chain argument");
	}

	// convert the first 3 columns of the data frame
	GIntervals intervs;
	convert_rintervs(rchain, &intervs, NULL);

	// convert the rest of the columns
	SEXP src_chroms = VECTOR_ELT(rchain, ChainInterval::CHROM_SRC);
	SEXP src_chrom_levels = Rf_getAttrib(src_chroms, R_LevelsSymbol);
	SEXP src_starts = VECTOR_ELT(rchain, ChainInterval::START_SRC);
	// Note: src_ends (END_SRC) is not used as it's calculated from start_src + (end - start) in constructor
	SEXP src_strands = VECTOR_ELT(rchain, ChainInterval::STRAND_SRC);
	SEXP tgt_strands = VECTOR_ELT(rchain, ChainInterval::STRAND);

	for (unsigned i = 0; i < ChainInterval::NUM_COLS; i++) {
		if (i != 0 && Rf_length(VECTOR_ELT(rchain, i)) != Rf_length(VECTOR_ELT(rchain, i - 1)))
			verror("Number of rows in column %s differs than the number of rows in column %s", ChainInterval::COL_NAMES[i - 1], ChainInterval::COL_NAMES[i]);
	}

	if (!Rf_isReal(src_starts) && !Rf_isInteger(src_starts))
		verror("Invalid format of intervals argument");

	unordered_map<string, int> src_chrom2id;

	for (unsigned i = 0; i < intervs.size(); i++) {
		if ((Rf_isFactor(src_chroms) && INTEGER(src_chroms)[i] < 0) || (Rf_isReal(src_starts) && std::isnan(REAL(src_starts)[i])))
			verror("Invalid format of interval at index %d", i + 1);

		const char *src_chrom = Rf_isString(src_chroms) ? CHAR(STRING_ELT(src_chroms, i)) : CHAR(STRING_ELT(src_chrom_levels, INTEGER(src_chroms)[i] - 1));
		unordered_map<string, int>::const_iterator isrc_chrom2id = src_chrom2id.find(src_chrom);
		int src_chromid;

		if (isrc_chrom2id == src_chrom2id.end()) {
			src_chromid = src_id2chrom.size();
			src_id2chrom.push_back(src_chrom);
			src_chrom2id[src_chrom] = src_chromid;
		} else
			src_chromid = isrc_chrom2id->second;

		int64_t src_start = (int64_t)(Rf_isReal(src_starts) ? REAL(src_starts)[i] : INTEGER(src_starts)[i]);

		// Read strands and convert +1/-1 to 0/1 for internal storage
		int tgt_strand_val = Rf_isReal(tgt_strands) ? (int)REAL(tgt_strands)[i] : INTEGER(tgt_strands)[i];
		int src_strand_val = Rf_isReal(src_strands) ? (int)REAL(src_strands)[i] : INTEGER(src_strands)[i];
		int tgt_strand = (tgt_strand_val == 1) ? 0 : 1;
		int src_strand = (src_strand_val == 1) ? 0 : 1;

		ChainInterval interval(intervs[i].chromid, intervs[i].start, intervs[i].end, tgt_strand, src_chromid, src_start, src_strand);
		interval.chain_id = i;  // Assign monotonically increasing ID for stable sorting

		// Validate that chain intervals are non-zero length (strict validation)
		if (interval.start >= interval.end)
			verror("Chain file contains zero-length or invalid interval at row %d (start=%lld, end=%lld)",
			       i + 1, (long long)interval.start, (long long)interval.end);

		interval.verify(m_chrom_key, src_id2chrom);
		chain_intervs.push_back(interval);
	}
}

SEXP IntervUtils::convert_chain_intervs(const ChainIntervals &chain_intervs, vector<string> &src_id2chrom)
{
	GIntervals tmp_intervals;
	tmp_intervals.reserve(chain_intervs.size());
	for (ChainIntervals::const_iterator iinterval = chain_intervs.begin(); iinterval != chain_intervs.end(); ++iinterval)
		tmp_intervals.push_back((GInterval)*iinterval);

    SEXP answer = convert_intervs(&tmp_intervals, ChainInterval::NUM_COLS);
	SEXP src_chroms, src_chroms_idx, src_starts, src_ends, src_strands, tgt_strands;
    SEXP col_names = Rf_getAttrib(answer, R_NamesSymbol);
    rprotect(col_names);
	unsigned num_src_chroms = src_id2chrom.size();

    rprotect(src_chroms_idx = RSaneAllocVector(INTSXP, chain_intervs.size()));
    rprotect(src_starts = RSaneAllocVector(REALSXP, chain_intervs.size()));
    rprotect(src_ends = RSaneAllocVector(REALSXP, chain_intervs.size()));
    rprotect(src_chroms = RSaneAllocVector(STRSXP, num_src_chroms));
    rprotect(src_strands = RSaneAllocVector(INTSXP, chain_intervs.size()));
    rprotect(tgt_strands = RSaneAllocVector(INTSXP, chain_intervs.size()));

	for (ChainIntervals::const_iterator iinterval = chain_intervs.begin(); iinterval != chain_intervs.end(); ++iinterval) {
		INTEGER(src_chroms_idx)[iinterval - chain_intervs.begin()] = iinterval->chromid_src + 1;
		REAL(src_starts)[iinterval - chain_intervs.begin()] = iinterval->start_src;
		REAL(src_ends)[iinterval - chain_intervs.begin()] = iinterval->end_src;
		// Convert 0/1 to +1/-1 for R output
		INTEGER(src_strands)[iinterval - chain_intervs.begin()] = (iinterval->strand_src == 0) ? 1 : -1;
		INTEGER(tgt_strands)[iinterval - chain_intervs.begin()] = (iinterval->strand == 0) ? 1 : -1;
	}

	for (unsigned id = 0; id < num_src_chroms; ++id)
		SET_STRING_ELT(src_chroms, id, Rf_mkChar(src_id2chrom[id].c_str()));

    for (int i = 0; i < ChainInterval::NUM_COLS; i++)
		SET_STRING_ELT(col_names, i, Rf_mkChar(ChainInterval::COL_NAMES[i]));

    Rf_setAttrib(src_chroms_idx, R_LevelsSymbol, src_chroms);
    Rf_setAttrib(src_chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(answer, ChainInterval::STRAND, tgt_strands);
    SET_VECTOR_ELT(answer, ChainInterval::CHROM_SRC, src_chroms_idx);
    SET_VECTOR_ELT(answer, ChainInterval::START_SRC, src_starts);
    SET_VECTOR_ELT(answer, ChainInterval::END_SRC, src_ends);
    SET_VECTOR_ELT(answer, ChainInterval::STRAND_SRC, src_strands);

    runprotect(5);
    return answer;
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

SEXP IntervUtils::create_data_frame(int numrows, int numcols, SEXP attrs_src)
{
	SEXP answer, row_names, col_names;

    answer = rprotect_ptr(RSaneAllocVector(VECSXP, numcols));
    col_names = rprotect_ptr(RSaneAllocVector(STRSXP, numcols));
    row_names = rprotect_ptr(RSaneAllocVector(INTSXP, numrows));

	for (int i = 0; i < numrows; ++i)
		INTEGER(row_names)[i] = i + 1;

    if (attrs_src != R_NilValue) 
        Rf_copyMostAttrib(attrs_src, answer);

    Rf_setAttrib(answer, R_NamesSymbol, col_names);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, row_names);

	return answer;
}

void IntervUtils::define_data_frame_cols(SEXP src, vector<SEXP> &src_cols, SEXP tgt, vector<SEXP> &tgt_cols, int tgt_col_offset)
{
	SEXP src_class = Rf_getAttrib(src, R_ClassSymbol);

	if (Rf_isNull(src_class) || !Rf_isString(src_class) ||
		(!(Rf_length(src_class) == 1 && !strcmp(CHAR(STRING_ELT(src_class, 0)), "data.frame")) &&
		!(Rf_length(src_class) == 3 && !strcmp(CHAR(STRING_ELT(src_class, 0)), "tbl_df")  &&
          !strcmp(CHAR(STRING_ELT(src_class, 1)), "tbl") && !strcmp(CHAR(STRING_ELT(src_class, 2)), "data.frame"))))
		verror("Copied object is not a data frame or tibble");

	if (Rf_length(tgt) < Rf_length(src) + tgt_col_offset)
		verror("Attempt to copy data frame columns beyond the valid size");

	int numrows = Rf_length(Rf_getAttrib(tgt, R_RowNamesSymbol));
    SEXP src_colnames = rprotect_ptr(Rf_getAttrib(src, R_NamesSymbol));
    SEXP tgt_colnames = rprotect_ptr(Rf_getAttrib(tgt, R_NamesSymbol));

	if (Rf_isNull(src_colnames) || !Rf_isString(src_colnames))
		verror("Invalid source data frame for a copy");

	src_cols.resize(Rf_length(src));
	if (tgt_cols.size() < (uint64_t)(Rf_length(tgt) + tgt_col_offset)){ 
		tgt_cols.resize(Rf_length(tgt) + tgt_col_offset);
	}

	for (int col = 0; col < Rf_length(src); ++col) {
		SEXP src_col = VECTOR_ELT(src, col);
		SEXP tgt_col;

        tgt_col = rprotect_ptr(RSaneAllocVector(TYPEOF(src_col), numrows));

		if (!Rf_isInteger(src_col) && !Rf_isReal(src_col) && !Rf_isLogical(src_col) && !Rf_isString(src_col) && !Rf_isFactor(src_col))
			verror("Unsupported type found in a data frame: %s", Rf_type2char(TYPEOF(src_col)));

		Rf_copyMostAttrib(src_col, tgt_col);
		SET_STRING_ELT(tgt_colnames, col + tgt_col_offset, STRING_ELT(src_colnames, col));
		src_cols[col] = src_col;
		tgt_cols[col + tgt_col_offset] = tgt_col;

        SET_VECTOR_ELT(tgt, col + tgt_col_offset, tgt_col);
    }
}

void IntervUtils::copy_data_frame_row(const vector<SEXP> &src_cols, int src_row, const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset)
{
	for (uint64_t col = 0; col < src_cols.size(); ++col) {
		SEXP src_col = src_cols[col];
		SEXP tgt_col = tgt_cols[col + tgt_col_offset];

		if (Rf_isInteger(src_col) || Rf_isFactor(src_col))
			INTEGER(tgt_col)[tgt_row] = INTEGER(src_col)[src_row];
		else if (Rf_isReal(src_col))
			REAL(tgt_col)[tgt_row] = REAL(src_col)[src_row];
		else if (Rf_isLogical(src_col))
			LOGICAL(tgt_col)[tgt_row] = LOGICAL(src_col)[src_row];
		else if (Rf_isString(src_col))
			SET_STRING_ELT(tgt_col, tgt_row, Rf_mkChar(CHAR(STRING_ELT(src_col, src_row))));
	}
}

void IntervUtils::copy_data_frame_rows(const vector<SEXP> &src_cols, int src_row, int num_rows, const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset)
{
	for (uint64_t col = 0; col < src_cols.size(); ++col) {
		SEXP src_col = src_cols[col];
		SEXP tgt_col = tgt_cols[col + tgt_col_offset];

		if (Rf_isInteger(src_col) || Rf_isFactor(src_col)) {
			int *src_vals = INTEGER(src_col);
			int *tgt_vals = INTEGER(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isReal(src_col)) {
			double *src_vals = REAL(src_col);
			double *tgt_vals = REAL(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isLogical(src_col)) {
			int *src_vals = LOGICAL(src_col);
			int *tgt_vals = LOGICAL(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isString(src_col)) {
			for (int i = 0; i < num_rows; ++i) 
				SET_STRING_ELT(tgt_col, tgt_row + i, Rf_mkChar(CHAR(STRING_ELT(src_col, src_row + i))));
		}
	}
}

void IntervUtils::set_data_frame_val_nan(const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col)
{
	SEXP rtgt_col = tgt_cols[tgt_col];

	if (Rf_isInteger(rtgt_col) || Rf_isFactor(rtgt_col))
		INTEGER(rtgt_col)[tgt_row] = NA_INTEGER;
	else if (Rf_isReal(rtgt_col))
		REAL(rtgt_col)[tgt_row] = NA_REAL;
	else if (Rf_isLogical(rtgt_col))
		LOGICAL(rtgt_col)[tgt_row] = NA_LOGICAL;
	else if (Rf_isString(rtgt_col))
		SET_STRING_ELT(rtgt_col, tgt_row, NA_STRING);
}

void IntervUtils::restrict_bins(int64_t maxbins, GIntervals &intervals, unsigned binsize) const
{
	for (GIntervals::const_iterator iinterval = intervals.begin(); iinterval != intervals.end(); ++iinterval) {
		int64_t bins = max((int64_t)0, (int64_t)ceil(iinterval->end / binsize) - (int64_t)(iinterval->start / binsize));

		if (bins > maxbins)
			verror("The interval %s [%ld, %ld) covers too wide range of samples that might cause memory allocation failure.\n"
					"(bins covered: %ld, bins limit: %ld)\n", id2chrom(iinterval->chromid).c_str(), iinterval->start, iinterval->end, bins, maxbins);
	}
}

bool IntervUtils::get_multitasking() const
{
	if (m_multitasking < 0) {
        SEXP r_multitasking = Rf_GetOption1(Rf_install("gmultitasking"));

		if (Rf_isLogical(r_multitasking))
			m_multitasking = (int)LOGICAL(r_multitasking)[0];
		else
			m_multitasking = false;
	}
	return (bool)m_multitasking;
}

uint64_t IntervUtils::get_max_processes() const
{
	if (!m_max_processes) {
        SEXP r_max_processes = Rf_GetOption1(Rf_install("gmax.processes"));

		if (Rf_isReal(r_max_processes))
			m_max_processes = (uint64_t)REAL(r_max_processes)[0];
		else if (Rf_isInteger(r_max_processes))
			m_max_processes = INTEGER(r_max_processes)[0];
		else
			m_max_processes = 64;
		if (m_max_processes < 1) 
			m_max_processes = 64;
	}
	return m_max_processes;
}

uint64_t IntervUtils::get_max_processes2core() const
{
	if (!m_max_processes2core) {
        SEXP r_max_processes2core = Rf_GetOption1(Rf_install("gmax.processes2core"));

		if (Rf_isReal(r_max_processes2core))
			m_max_processes2core = (uint64_t)REAL(r_max_processes2core)[0];
		else if (Rf_isInteger(r_max_processes2core))
			m_max_processes2core = INTEGER(r_max_processes2core)[0];
		else
			m_max_processes2core = 4;
		if (m_max_processes2core < 1) 
			m_max_processes2core = 4;
	}
	return m_max_processes2core;
}

uint64_t IntervUtils::get_min_scope4process() const
{
	if (!m_min_scope4process) {
        SEXP r_min_scope4process = Rf_GetOption1(Rf_install("gmin.scope4process"));

		if (Rf_isReal(r_min_scope4process))
			m_min_scope4process = (uint64_t)REAL(r_min_scope4process)[0];
		else if (Rf_isInteger(r_min_scope4process))
			m_min_scope4process = INTEGER(r_min_scope4process)[0];
		else
			m_min_scope4process = 10000;
	}
	return m_min_scope4process;
}

uint64_t IntervUtils::get_min_seqs_work4process() const
{
	if (!m_min_seqs_work4process) {
        SEXP r_min_seqs_work4process = Rf_GetOption1(Rf_install("gmin.seqs.work4process"));

		if (Rf_isReal(r_min_seqs_work4process))
			m_min_seqs_work4process = (uint64_t)REAL(r_min_seqs_work4process)[0];
		else if (Rf_isInteger(r_min_seqs_work4process))
			m_min_seqs_work4process = INTEGER(r_min_seqs_work4process)[0];
		else
			m_min_seqs_work4process = 100000;
	}
	return m_min_seqs_work4process;
}

uint64_t IntervUtils::get_max_data_size() const
{
	if (!m_max_data_size) {
        SEXP r_max_data_size = Rf_GetOption1(Rf_install("gmax.data.size"));

		if (Rf_isReal(r_max_data_size))
			m_max_data_size = (uint64_t)REAL(r_max_data_size)[0];
		else if (Rf_isInteger(r_max_data_size))
			m_max_data_size = INTEGER(r_max_data_size)[0];
		else
			m_max_data_size = numeric_limits<uint64_t>::max();
	}
	return m_max_data_size;
}

uint64_t IntervUtils::get_max_mem_usage() const
{
	if (!m_max_mem_usage) {
        SEXP r_max_mem_usage = Rf_GetOption1(Rf_install("gmax.mem.usage"));

		if (Rf_isReal(r_max_mem_usage))
			m_max_mem_usage = (uint64_t)REAL(r_max_mem_usage)[0] * 1000;
		else if (Rf_isInteger(r_max_mem_usage))
			m_max_mem_usage = INTEGER(r_max_mem_usage)[0] * 1000;
		else
			m_max_mem_usage = numeric_limits<uint64_t>::max();
	}
	return m_max_mem_usage;
}

uint64_t IntervUtils::get_big_intervals_size() const
{
	if (!m_big_intervals_size) {
        SEXP r_big_intervals_size = Rf_GetOption1(Rf_install("gbig.intervals.size"));

		if (Rf_isReal(r_big_intervals_size))
			m_big_intervals_size = (uint64_t)REAL(r_big_intervals_size)[0];
		else if (Rf_isInteger(r_big_intervals_size))
			m_big_intervals_size = INTEGER(r_big_intervals_size)[0];
		else
			m_big_intervals_size = numeric_limits<uint64_t>::max();
		m_big_intervals_size = min(m_big_intervals_size, get_max_data_size());
	}
	return m_big_intervals_size;
}

uint64_t IntervUtils::get_quantile_edge_data_size() const
{
	if (!m_quantile_edge_data_size) {
        SEXP r_quantile_edge_data_size = Rf_GetOption1(Rf_install("gquantile.edge.data.size"));

		if (Rf_isReal(r_quantile_edge_data_size))
			m_quantile_edge_data_size = (uint64_t)REAL(r_quantile_edge_data_size)[0];
		else if (Rf_isInteger(r_quantile_edge_data_size))
			m_quantile_edge_data_size = INTEGER(r_quantile_edge_data_size)[0];
		else
			m_quantile_edge_data_size = 0;
	}
	return m_quantile_edge_data_size;
}

uint64_t IntervUtils::get_track_chunk_size() const
{
	if (!m_track_chunk_size) {
        SEXP r_track_chunk_size = Rf_GetOption1(Rf_install("gtrack.chunk.size"));

		if (Rf_isReal(r_track_chunk_size))
			m_track_chunk_size = (uint64_t)REAL(r_track_chunk_size)[0];
		else if (Rf_isInteger(r_track_chunk_size))
			m_track_chunk_size = INTEGER(r_track_chunk_size)[0];
		else
			m_track_chunk_size = 100000;
	}
	return m_track_chunk_size;
}

uint64_t IntervUtils::get_track_num_chunks() const
{
	if (!m_track_num_chunks) {
        SEXP r_track_num_chunks = Rf_GetOption1(Rf_install("gtrack.num.chunks"));

		if (Rf_isReal(r_track_num_chunks))
			m_track_num_chunks = (uint64_t)REAL(r_track_num_chunks)[0];
		else if (Rf_isInteger(r_track_num_chunks))
			m_track_num_chunks = INTEGER(r_track_num_chunks)[0];
		else
			m_track_num_chunks = 0;
	}
	return m_track_num_chunks;
}

bool IntervUtils::is_1d_iterator(SEXP rtrack_exprs, GIntervalsFetcher1D *scope1d, GIntervalsFetcher2D *scope2d, SEXP riterator)
{
	TrackExprScanner scanner(*this);
	TrackExpressionIteratorBase *itr = scanner.create_expr_iterator(rtrack_exprs, scope1d, scope2d, riterator, R_NilValue, false);
	return itr->is_1d();
}

void IntervUtils::verify_max_data_size(uint64_t data_size, const char *data_name, bool check_all_kids)
{
	if (data_size > get_max_data_size())
		verror("%s size exceeded the maximal allowed (%ld).\n"
				"Try to bound the scope of the function.\n"
				"Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).",
			   data_name, get_max_data_size());

	if (check_all_kids) 
		update_res_data_size(data_size);
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

bool IntervUtils::distribute_task(uint64_t res_const_size,    // data size in bytes for all the result
								  uint64_t res_record_size)   // size in bytes per datum in the result
{
	uint64_t max_res_size = get_max_data_size() * res_record_size + m_num_planned_kids * res_const_size;

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

		if (iinterv->chromid_src == (iinterv - 1)->chromid_src && (iinterv - 1)->start_src + (iinterv - 1)->end - (iinterv - 1)->start > iinterv->start_src)
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
			    (*this)[i-1].start_src + (*this)[i-1].end - (*this)[i-1].start > (*this)[i].start_src) {
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
	if (policy == "error") {
		verify_no_tgt_overlaps(chromkey, src_id2chrom);
	} else if (policy == "auto") {
		// Auto-resolve by truncating/splitting intervals
		if (empty())
			return;

		set<ChainInterval, ChainInterval::SetCompare> sorted_intervs;
		for (iterator iinterv = begin(); iinterv != end(); ++iinterv)
			sorted_intervs.insert(*iinterv);

		set<ChainInterval>::iterator iinterv2 = sorted_intervs.begin();
		set<ChainInterval>::iterator iinterv1 = iinterv2++;

		while (iinterv2 != sorted_intervs.end()) {
			if (iinterv1->chromid == iinterv2->chromid && iinterv1->end > iinterv2->start) {
				// Overlapping intervals detected - truncate/split
				int64_t tgt_end1 = iinterv1->end;

				// Truncate interv1
				((ChainInterval &)*iinterv1).end = iinterv2->start;

				// Adjust or split interv2
				if (tgt_end1 < iinterv2->end) {
					// The two intervals intersect - create non-overlapping interv2
					ChainInterval interv(iinterv2->chromid, tgt_end1, iinterv2->end, iinterv2->strand,
							iinterv2->chromid_src, iinterv2->start_src + tgt_end1 - iinterv2->start, iinterv2->strand_src);
					interv.chain_id = iinterv2->chain_id;  // Preserve chain_id for deterministic sorting
					sorted_intervs.erase(iinterv2);

					if (interv.start != interv.end)
						sorted_intervs.insert(interv);
				} else {
					// interval1 contains interval2 => split interval1
					ChainInterval interv(iinterv1->chromid, iinterv1->end + iinterv2->end - iinterv2->start, tgt_end1, iinterv1->strand,
							iinterv1->chromid_src, iinterv1->start_src + iinterv2->end - iinterv1->start, iinterv1->strand_src);
					interv.chain_id = iinterv1->chain_id;  // Preserve chain_id for deterministic sorting
					sorted_intervs.erase(iinterv2);
					if (interv.start != interv.end)
						sorted_intervs.insert(interv);
				}

				// Remove zero-length intervals
				if (iinterv1->start == iinterv1->end) {
					iinterv2 = iinterv1;
					--iinterv2;
					sorted_intervs.erase(iinterv1);
				} else
					iinterv2 = iinterv1;
			}
			iinterv1 = iinterv2;
			++iinterv2;
		}

		// Copy back the resolved intervals
		clear();
		for (set<ChainInterval>::const_iterator iinterval = sorted_intervs.begin(); iinterval != sorted_intervs.end(); ++iinterval)
			push_back(*iinterval);
	} else if (policy == "discard") {
		// Mark all intervals that have overlapping targets
		vector<bool> to_discard(size(), false);

		for (size_t i = 1; i < size(); ++i) {
			if ((*this)[i-1].chromid == (*this)[i].chromid &&
			    (*this)[i-1].end > (*this)[i].start) {
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
	} else if (policy == "keep") {
		// Do nothing - allow overlapping target intervals
	} else {
		TGLError("Invalid target overlap policy: %s. Must be 'error', 'keep', 'auto', or 'discard'", policy.c_str());
	}
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
			const int64_t end_src = ci.start_src + (ci.end - ci.start);
			pmax = std::max(pmax, end_src);
			m_pmax_end_src[i] = pmax;
		}
		m_chrom_last_excl[chrom] = i; // exclusive end
	}
}

ChainIntervals::const_iterator ChainIntervals::map_interval(const GInterval &src_interval, GIntervals &tgt_intervs, ChainIntervals::const_iterator hint)
{
	tgt_intervs.clear();

	if (empty())
		return end();

	if (front().chromid_src > src_interval.chromid || (front().chromid_src == src_interval.chromid && front().start_src >= src_interval.end))
		return begin();

	if (back().chromid_src < src_interval.chromid || (back().chromid_src == src_interval.chromid && back().start_src + back().end - back().start <= src_interval.start))
		return end() - 1;

	if (check_first_overlap_src(hint, src_interval))
		return add2tgt(hint, src_interval, tgt_intervs);

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
				return add2tgt(first_overlap, src_interval, tgt_intervs);
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
			return add2tgt(begin() + lo, src_interval, tgt_intervs);
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
			return add2tgt(iend_interval, src_interval, tgt_intervs);
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
			const_iterator last_checked = add2tgt(begin() + lo, src_interval, tgt_intervs);
			// Optimize hint for next query: return pos if it's ahead and in same chrom
			if (pos >= (size_t)(last_checked - begin()) && pos < lastEx)
				return begin() + pos;
			return last_checked;
		}
	}

	// No left overlap exists; check the right neighbor (iend_interval)
	if (iend_interval != end() && iend_interval->do_overlap_src(src_interval))
		return add2tgt(iend_interval, src_interval, tgt_intervs);

	// Nothing overlaps
	return istart_interval;
}

ChainIntervals::const_iterator ChainIntervals::add2tgt(const_iterator hint, const GInterval &src_interval, GIntervals &tgt_intervs)
{
	const_iterator last_checked = hint;

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

			// Debug assertions to catch coordinate corruption
			assert(tgt_start < tgt_end && "Target interval has zero or negative length");
			assert(tgt_start >= hint->start && tgt_end <= hint->end && "Target coordinates out of chain bounds");

			tgt_intervs.push_back(GInterval(hint->chromid, tgt_start, tgt_end, 0));
		}

		last_checked = hint;
		++hint;
	}
	return last_checked;
}
