#include "IntervalConverter.h"
#include "rdbinterval.h" // For IntervUtils, ChainIntervals, GIntervals, GIntervals2D, etc.
#include "rdbutils.h" // For verror, rprotect, runprotect, RSaneAllocVector, monitor_memusage, find_in_misha, interv2path, track2path, eval_in_R, RSaneUnserialize
#include "GIntervals.h"
#include "GIntervals2D.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "GTrackIntervalsFetcher1D.h"
#include "GTrackIntervalsFetcher2D.h"
#include "GenomeTrack.h"
#include "TGLException.h"
#include <cstring> // For strcmp
#include <cmath> // For isnan
#include <sys/stat.h> // For stat, S_ISDIR
#include <unistd.h> // For access, errno, ENOENT
#include <string>

using namespace std;
using namespace rdb;

IntervalConverter::IntervalConverter(rdb::IntervUtils &iu) :
	m_iu(iu)
{
}

unsigned IntervalConverter::get_rintervs_type_mask(SEXP rintervals, const char *error_msg_prefix) const
{
	if (!Rf_isVector(rintervals))
		verror("%sInvalid format of intervals argument", error_msg_prefix);

	if (Rf_length(rintervals) == 2) {
		if (get_rintervs_type_mask(VECTOR_ELT(rintervals, 0), error_msg_prefix) != rdb::IntervUtils::INTERVS1D || get_rintervs_type_mask(VECTOR_ELT(rintervals, 1), error_msg_prefix) != rdb::IntervUtils::INTERVS2D)
			verror("%sInvalid format of intervals argument", error_msg_prefix);
		return rdb::IntervUtils::INTERVS1D | rdb::IntervUtils::INTERVS2D;
	}

	SEXP colnames = Rf_getAttrib(rintervals, R_NamesSymbol);

	if (!Rf_isString(colnames) || Rf_length(colnames) < GInterval::NUM_COLS)
		verror("%sInvalid format of intervals argument", error_msg_prefix);

	rdb::IntervUtils::IntervsType type = rdb::IntervUtils::INTERVS1D;

	for (unsigned i = 0; i < GInterval::NUM_COLS; i++) {
		if (strcmp(CHAR(STRING_ELT(colnames, i)), GInterval::COL_NAMES[i])) {
			type = rdb::IntervUtils::INTERVS2D;
			break;
		}
	}

	if (type == rdb::IntervUtils::INTERVS2D) {
		for (unsigned i = 0; i < GInterval2D::NUM_COLS; i++) {
			if (strcmp(CHAR(STRING_ELT(colnames, i)), GInterval2D::COL_NAMES[i]))
				verror("Invalid format of intervals: column names do not match neither 1d nor 2d intervals");
		}
	}

	if (type == rdb::IntervUtils::INTERVS1D) {
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

	} else if (type == rdb::IntervUtils::INTERVS2D) {
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

unsigned IntervalConverter::convert_rintervs(SEXP rintervals, GIntervalsFetcher1D **intervals1d, GIntervalsFetcher2D **intervals2d, bool null_if_interv_nonexist,
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

			if (GIntervalsBigSet::isbig(full_interv_name, m_iu)) { // big intervals set?
				SEXP meta = GIntervalsMeta::load_meta(interv2path(m_iu.get_env(), full_interv_name).c_str());

				if (GIntervalsMeta1D::is1d(meta)) {
					if (intervals1d) {
						*intervals1d = new GIntervalsBigSet1D(full_interv_name, meta, m_iu);
						if (intervals2d)
							*intervals2d = new GIntervals2D();  // create an empty intervals set for convenience
						return rdb::IntervUtils::INTERVS1D;
					}
					verror("%sExpecting 2D intervals while receiving 1D intervals set %s", error_msg_prefix, full_interv_name);
				} else {
					if (intervals2d) {
						*intervals2d = new GIntervalsBigSet2D(full_interv_name, meta, m_iu);
						if (intervals1d)
							*intervals1d = new GIntervals();  // create an empty intervals set for convenience
						return rdb::IntervUtils::INTERVS2D;
					}
					verror("%sExpecting 1D intervals while receiving 2D intervals set %s", error_msg_prefix, full_interv_name);
				}
			} else if (GTrackIntervalsFetcher::isbig(full_interv_name, m_iu)) {
				string trackpath = track2path(m_iu.get_env(), full_interv_name);

				if (access((trackpath + "/.meta").c_str(), R_OK) && errno == ENOENT)
					GTrackIntervalsFetcher::create_track_meta(full_interv_name, m_iu);

				SEXP meta = GIntervalsMeta::load_meta(trackpath.c_str());

				if (GIntervalsMeta1D::is1d(meta)) {
					if (intervals1d) {
						GenomeTrack::Type track_type = GenomeTrack::get_type(trackpath.c_str(), m_iu.get_chromkey(), true);

						if (track_type == GenomeTrack::FIXED_BIN)
							verror("%s is track of %s type which cannot be used as a replacement for intervals",
								   full_interv_name, GenomeTrack::TYPE_NAMES[track_type]);
						else if (track_type == GenomeTrack::SPARSE)
							*intervals1d = new GTrackIntervalsFetcher1D<GenomeTrackSparse>(full_interv_name, meta, m_iu);
						else if (track_type == GenomeTrack::ARRAYS)
							*intervals1d = new GTrackIntervalsFetcher1D<GenomeTrackArrays>(full_interv_name, meta, m_iu);

						if (intervals2d)
							*intervals2d = new GIntervals2D();  // create an empty intervals set for convenience
						return rdb::IntervUtils::INTERVS1D;
					}
					verror("%sExpecting 2D intervals while receiving 1D intervals set %s", error_msg_prefix, full_interv_name);
				} else {
					if (intervals2d) {
						GenomeTrack::Type track_type = GenomeTrack::get_type(trackpath.c_str(), m_iu.get_chromkey(), true);

						if (track_type == GenomeTrack::RECTS)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackRectsRects>(full_interv_name, meta, m_iu);
						else if (track_type == GenomeTrack::POINTS)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackRectsPoints>(full_interv_name, meta, m_iu);
						else if (track_type == GenomeTrack::COMPUTED)
							*intervals2d = new GTrackIntervalsFetcher2D<GenomeTrackComputed>(full_interv_name, meta, m_iu);

						if (intervals1d)
							*intervals1d = new GIntervals();  // create an empty intervals set for convenience
						return rdb::IntervUtils::INTERVS2D;
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

SEXP IntervalConverter::convert_rintervs(SEXP rintervals, GIntervals *intervals, GIntervals2D *intervals2d, bool null_if_interv_nonexist,
										  const GenomeChromKey *chromkey, const char *error_msg_prefix, unsigned *pintervs_type_mask, 
										  bool verify, bool skip_missing_chroms) const
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
            rprotect(gintervs = find_in_misha(m_iu.get_env(), "GINTERVS"));
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
			if (GTrackIntervalsFetcher::isbig(full_interv_name, m_iu))
				verror("Tracks cannot be used as a substitute for intervals in this function");

			if (null_if_interv_nonexist)
				return R_NilValue;

			verror("%sInterval %s does not exist", error_msg_prefix, full_interv_name);
		}

		string path = interv2path(m_iu.get_env(), full_interv_name);
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

	// Handle empty intervals - determine type based on what the caller expects
	// For empty intervals, the caller's GIntervals/GIntervals2D objects are already
	// default-constructed as empty, so we just need to return the appropriate type mask
	if (Rf_isVector(rintervals) && Rf_length(rintervals) == 0) {
		unsigned empty_intervs_type_mask;
		if (intervals && !intervals2d) {
			empty_intervs_type_mask = rdb::IntervUtils::INTERVS1D;
		} else if (intervals2d && !intervals) {
			empty_intervs_type_mask = rdb::IntervUtils::INTERVS2D;
		} else {
			// Both or neither provided - default to 1D
			empty_intervs_type_mask = rdb::IntervUtils::INTERVS1D;
		}

		if (pintervs_type_mask)
			*pintervs_type_mask = empty_intervs_type_mask;

		// Return early - the caller's intervals objects are already empty by default
		return _rintervals;
	}

	unsigned intervs_type_mask = get_rintervs_type_mask(rintervals, error_msg_prefix);

	if (pintervs_type_mask)
		*pintervs_type_mask = intervs_type_mask;

	if (intervs_type_mask == (rdb::IntervUtils::INTERVS1D | rdb::IntervUtils::INTERVS2D))
		rintervals = VECTOR_ELT(_rintervals, 0);

	if (intervs_type_mask == rdb::IntervUtils::INTERVS1D && !intervals)
		verror("%sExpecting 2D intervals while receiving 1D intervals", error_msg_prefix);

	if (intervs_type_mask == rdb::IntervUtils::INTERVS2D && !intervals2d)
		verror("%sExpecting 1D intervals while receiving 2D intervals", error_msg_prefix);

	if ((intervs_type_mask & rdb::IntervUtils::INTERVS1D) && intervals) {
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
			int chromid;
			try {
				chromid = chromkey ? chromkey->chrom2id(chrom) : m_iu.chrom2id(chrom);
			} catch (const TGLException &) {
				// If skip_missing_chroms is enabled, skip intervals with chromosomes not in the key
				if (skip_missing_chroms && chromkey)
					continue;
				throw;
			}
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
				interval.verify(chromkey ? *chromkey : m_iu.get_chromkey());
			intervals->push_back(interval);
		}
	}

	if (intervs_type_mask == (rdb::IntervUtils::INTERVS1D | rdb::IntervUtils::INTERVS2D))
		rintervals = VECTOR_ELT(_rintervals, 1);

	if ((intervs_type_mask & rdb::IntervUtils::INTERVS2D) && intervals2d) {
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
			int chromid1, chromid2;
			try {
				chromid1 = chromkey ? chromkey->chrom2id(chrom1) : m_iu.chrom2id(chrom1);
				chromid2 = chromkey ? chromkey->chrom2id(chrom2) : m_iu.chrom2id(chrom2);
			} catch (const TGLException &) {
				// If skip_missing_chroms is enabled, skip intervals with chromosomes not in the key
				if (skip_missing_chroms && chromkey)
					continue;
				throw;
			}
			int64_t start1 = (int64_t)(Rf_isReal(starts1) ? REAL(starts1)[i] : INTEGER(starts1)[i]);
			int64_t start2 = (int64_t)(Rf_isReal(starts2) ? REAL(starts2)[i] : INTEGER(starts2)[i]);
			int64_t end1 = (int64_t)(Rf_isReal(ends1) ? REAL(ends1)[i] : INTEGER(ends1)[i]);
			int64_t end2 = (int64_t)(Rf_isReal(ends2) ? REAL(ends2)[i] : INTEGER(ends2)[i]);

			GInterval2D interval(chromid1, start1, end1, chromid2, start2, end2, (void *)(intptr_t)i);

			if (verify)
				interval.verify(chromkey ? *chromkey : m_iu.get_chromkey());
			intervals2d->push_back(interval);
		}
	}

	return _rintervals;
}

SEXP IntervalConverter::convert_intervs(GIntervalsFetcher1D *intervals, unsigned num_cols, bool null_if_empty, bool use_original_index) const
{
	monitor_memusage();

	if (null_if_empty && !intervals->size())
		return R_NilValue;

	unsigned num_chroms = m_iu.get_chromkey().get_num_chroms();

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
		uint64_t index = use_original_index ? m_iu.get_orig_interv_idx(interval) : intervals->iter_index();

		INTEGER(chroms_idx)[index] = interval.chromid + 1;
		REAL(starts)[index] = interval.start;
		REAL(ends)[index] = interval.end;
		INTEGER(row_names)[index] = index + 1;
	}

	for (unsigned id = 0; id < (unsigned)num_chroms; ++id)
		SET_STRING_ELT(chroms, id, Rf_mkChar(m_iu.get_chromkey().id2chrom(id).c_str()));

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

SEXP IntervalConverter::convert_intervs(GIntervalsFetcher2D *intervals, unsigned num_cols, bool null_if_empty, bool use_original_index) const
{
	monitor_memusage();

	if (null_if_empty && !intervals->size())
		return R_NilValue;

	unsigned num_chroms = m_iu.get_chromkey().get_num_chroms();

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
		uint64_t index = use_original_index ? m_iu.get_orig_interv_idx(interval) : intervals->iter_index();

		INTEGER(chroms_idx1)[index] = interval.chromid1() + 1;
		REAL(starts1)[index] = interval.start1();
        REAL(ends1)[index] = interval.end1();
        INTEGER(chroms_idx2)[index] = interval.chromid2() + 1;
        REAL(starts2)[index] = interval.start2();
        REAL(ends2)[index] = interval.end2();
		INTEGER(row_names)[index] = index + 1;
	}

    for (unsigned id = 0; id < (unsigned)num_chroms; ++id) {
        SET_STRING_ELT(chroms1, id, Rf_mkChar(m_iu.get_chromkey().id2chrom(id).c_str()));
        SET_STRING_ELT(chroms2, id, Rf_mkChar(m_iu.get_chromkey().id2chrom(id).c_str()));
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

