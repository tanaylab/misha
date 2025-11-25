#include <cstdint>
#include <cmath>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>

#include "rdbutils.h"

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "GenomeTrack.h"
#include "GenomeTrackArrays.h"
#include "GenomeTrackComputed.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "TrackExpressionVars.h"
#include "TrackExpressionParams.h"
#include "GenomeSeqFetch.h"

const char *TrackExpressionVars::Track_var::FUNC_NAMES[TrackExpressionVars::Track_var::NUM_FUNCS] = {
	"avg", "min", "max", "nearest", "stddev", "sum", "quantile",
	"global.percentile", "global.percentile.min", "global.percentile.max",
	"weighted.sum", "area", "pwm", "pwm.max", "pwm.max.pos", "pwm.count", "kmer.count", "kmer.frac",
    "max.pos.abs", "max.pos.relative", "min.pos.abs", "min.pos.relative",
    "exists", "size", "sample", "sample.pos.abs", "sample.pos.relative",
    "first", "first.pos.abs", "first.pos.relative", "last", "last.pos.abs", "last.pos.relative"};

const char *TrackExpressionVars::Interv_var::FUNC_NAMES[TrackExpressionVars::Interv_var::NUM_FUNCS] = { "distance", "distance.center", "coverage", "neighbor.count" };

const char *TrackExpressionVars::Value_var::FUNC_NAMES[TrackExpressionVars::Value_var::NUM_FUNCS] = {
	"avg", "min", "max", "stddev", "sum", "quantile",
	"nearest",
	"exists", "size",
	"first", "last", "sample",
	"first.pos.abs", "first.pos.relative",
	"last.pos.abs", "last.pos.relative",
	"sample.pos.abs", "sample.pos.relative",
	"min.pos.abs", "min.pos.relative",
	"max.pos.abs", "max.pos.relative"
};

using namespace rdb;


TrackExpressionVars::TrackExpressionVars(rdb::IntervUtils &iu) :
	m_iu(iu)
{
	m_imdfs1d.reserve(10000);
	m_imdfs2d.reserve(10000);
	m_track_n_imdfs.reserve(10000);
	m_groot = get_groot(m_iu.get_env());
	m_shared_seqfetch.set_seqdir(m_groot + "/seq");
}

TrackExpressionVars::~TrackExpressionVars()
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ivar++)
		runprotect(ivar->rvar);
	for (Interv_vars::iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ivar++)
		runprotect(ivar->rvar);
	for (Value_vars::iterator ivar = m_value_vars.begin(); ivar != m_value_vars.end(); ivar++)
		runprotect(ivar->rvar);
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
		delete itrack_n_imdf->track;
}

void TrackExpressionVars::parse_exprs(const vector<string> &track_exprs)
{
	enum { TRACK, VTRACK, NUM_VAR_TYPES };

	SEXP rtracknames[NUM_VAR_TYPES] = { R_NilValue, R_NilValue };
	SEXP gvtracks = R_NilValue;
	SEXP vtracks = R_NilValue;

    // retrieve track names (split nested calls and protect environment)
    rtracknames[TRACK] = rprotect_ptr(find_in_misha(m_iu.get_env(), "GTRACKS"));

	// retrieve virtual track names (it's more complex since virtual track names are burried in a list of lists)
    rtracknames[VTRACK] = R_NilValue;
    gvtracks = rprotect_ptr(find_in_misha(m_iu.get_env(), "GVTRACKS"));

    if (!Rf_isNull(gvtracks) && !Rf_isSymbol(gvtracks)) { 
        SEXP gwds = rprotect_ptr(Rf_getAttrib(gvtracks, R_NamesSymbol));

		if (!Rf_isVector(gvtracks) || (Rf_length(gvtracks) && !Rf_isString(gwds)) || Rf_length(gwds) != Rf_length(gvtracks))
			verror("Invalid format of GVTRACKS variable.\n"
				"To continue working with virtual tracks please remove this variable from the environment.");

        const char *gwd = get_gwd(m_iu.get_env());

		for (int i = 0; i < Rf_length(gwds); ++i) {
			if (!strcmp(gwd, CHAR(STRING_ELT(gwds, i)))) {
                vtracks = VECTOR_ELT(gvtracks, i);
                SEXP vtracknames = Rf_getAttrib(vtracks, R_NamesSymbol);

				if (!Rf_isVector(vtracks) || (Rf_length(vtracks) && !Rf_isString(vtracknames)) || Rf_length(vtracknames) != Rf_length(vtracks))
					verror("Invalid format of GVTRACKS variable.\n"
						"To continue working with virtual tracks please remove this variable from the environment.");

                rtracknames[VTRACK] = vtracknames;
			}
		}
	}
    // Unprotect: rtracknames[TRACK], gvtracks (if non-null), gwds if we protected it
    if (!Rf_isNull(gvtracks) && !Rf_isSymbol(gvtracks)) {
        runprotect(2); // rtracknames[TRACK], gvtracks
    } else {
        runprotect(1); // rtracknames[TRACK]
    }

	for (vector<string>::const_iterator iexpr = track_exprs.begin(); iexpr != track_exprs.end(); ++iexpr) {
		for (int var_type = 0; var_type < NUM_VAR_TYPES; ++var_type) {
			if (!Rf_isString(rtracknames[var_type]))
				continue;

			for (int itrack = 0; itrack < Rf_length(rtracknames[var_type]); ++itrack) {
				string track = CHAR(STRING_ELT(rtracknames[var_type], itrack));
				uint64_t pos = 0;

				while ((pos = iexpr->find(track, pos)) != string::npos) {
					if (is_var(*iexpr, pos, pos + track.size())) {
						if (var_type == TRACK)
							add_track_var(track);
						else
                            add_vtrack_var(track, VECTOR_ELT(vtracks, itrack));
						break;
					}
					pos += track.size();
				}
			}
		}
	}

	for (Interv_vars::iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ++ivar) {
		ivar->siinterv = ivar->sintervs.begin();
		if (ivar->val_func == Interv_var::DIST || ivar->val_func == Interv_var::NEIGHBOR_COUNT)
			ivar->eiinterv = ivar->eintervs.begin();
	}

	// Value vars don't need initialization - they use track objects
}

void TrackExpressionVars::parse_imdf(SEXP rvtrack, const string &vtrack, Iterator_modifier1D *imdf1d, Iterator_modifier2D *imdf2d)
{
	SEXP rimdf = get_rvector_col(rvtrack, "itr", vtrack.c_str(), false);

	if (Rf_isNull(rimdf))
		return;

	string vname = vtrack + "$itr";

	SEXP rtype = get_rvector_col(rimdf, "type", vname.c_str(), true);

	if (!Rf_isString(rtype) || Rf_length(rtype) != 1)
		verror("Invalid format of virtual track %s", vtrack.c_str());

	string type(CHAR(STRING_ELT(rtype, 0)));
	transform(type.begin(), type.end(), type.begin(), ::tolower);

	if (type == "1d") {
		if (!imdf1d)
			verror("Virtual track %s: 1D iterator modifier cannot be used with source that supports only 2D iterators", vtrack.c_str());

		SEXP rdim = get_rvector_col(rimdf, "dim", vname.c_str(), false);

		if (Rf_isNull(rdim))
			imdf1d->dim = Iterator_modifier1D::DIM_NONE;
		else {
			if (!(Rf_isReal(rdim) || Rf_isInteger(rdim)) || Rf_length(rdim) != 1)
				verror("Virtual track %s: invalid dimension projection of iterator modifier", vtrack.c_str());

			double dim = Rf_isReal(rdim) ? REAL(rdim)[0] : INTEGER(rdim)[0];

			if (!dim)
				imdf1d->dim = Iterator_modifier1D::DIM_NONE;
			else if (dim == 1)
				imdf1d->dim = Iterator_modifier1D::DIM1;
			else if (dim == 2)
				imdf1d->dim = Iterator_modifier1D::DIM2;
			else
				verror("Virtual track %s: invalid dimension projection of iterator modifier", vtrack.c_str());
		}

		enum { SSHIFT, ESHIFT, NUM_SHIFTS };
		static const char *shift_col_names[NUM_SHIFTS] = { "sshift", "eshift" };

		for (int i = 0; i < NUM_SHIFTS; i++) {
			SEXP rshift = get_rvector_col(rimdf, shift_col_names[i], vname.c_str(), false);

			if (!(Rf_isReal(rshift) || Rf_isInteger(rshift)) || Rf_length(rshift) != 1)
				verror("Virtual track %s: %s must be an integer", vtrack.c_str(), shift_col_names[i]);

			int64_t shift = Rf_isReal(rshift) ? (int64_t)REAL(rshift)[0] : INTEGER(rshift)[0];

			if (i == SSHIFT)
				imdf1d->sshift = shift;
			else if (i == ESHIFT)
				imdf1d->eshift = shift;
		}
	} else if (type == "2d") {
		if (!imdf2d)
			verror("Virtual track %s: 2D iterator modifier cannot be used with source that supports only 1D iterators", vtrack.c_str());

		enum { SSHIFT1, ESHIFT1, SSHIFT2, ESHIFT2, NUM_SHIFTS };
		static const char *shift_col_names[NUM_SHIFTS] = { "sshift1", "eshift1", "sshift2", "eshift2" };

		for (int i = 0; i < NUM_SHIFTS; i++) {
			SEXP rshift = get_rvector_col(rimdf, shift_col_names[i], vname.c_str(), false);

			if (!(Rf_isReal(rshift) || Rf_isInteger(rshift)) || Rf_length(rshift) != 1)
				verror("Virtual track %s: %s must be an integer", vtrack.c_str(), shift_col_names[i]);

			int64_t shift = Rf_isReal(rshift) ? (int64_t)REAL(rshift)[0] : INTEGER(rshift)[0];

			if (i == SSHIFT1)
				imdf2d->sshift1 = shift;
			else if (i == ESHIFT1)
				imdf2d->eshift1 = shift;
			else if (i == SSHIFT2)
				imdf2d->sshift2 = shift;
			else if (i == ESHIFT2)
				imdf2d->eshift2 = shift;
		}
	} else
		verror("Virtual track %s: invalid type of iterator modifier", vtrack.c_str());
}

TrackExpressionVars::Iterator_modifier1D *TrackExpressionVars::add_imdf(const Iterator_modifier1D &imdf)
{
	if (imdf == Iterator_modifier1D())
		return NULL;

	for (Iterator_modifiers1D::iterator iimdf = m_imdfs1d.begin(); iimdf != m_imdfs1d.end(); ++iimdf) {
		if (*iimdf == imdf)
			return &*iimdf;
	}

	if (m_imdfs1d.size() == m_imdfs1d.capacity())
		verror("Reached the limit of maximal number of tracks");

	m_imdfs1d.push_back(imdf);
	return &m_imdfs1d.back();
}

TrackExpressionVars::Iterator_modifier2D *TrackExpressionVars::add_imdf(const Iterator_modifier2D &imdf)
{
	if (imdf == Iterator_modifier2D())
		return NULL;

	for (Iterator_modifiers2D::iterator iimdf = m_imdfs2d.begin(); iimdf != m_imdfs2d.end(); ++iimdf) {
		if (*iimdf == imdf)
			return &*iimdf;
	}

	if (m_imdfs2d.size() == m_imdfs2d.capacity())
		verror("Reached the limit of maximal number of tracks");

	m_imdfs2d.push_back(imdf);
	return &m_imdfs2d.back();
}

TrackExpressionVars::Track_n_imdf &TrackExpressionVars::add_track_n_imdf(const string &track, GenomeTrack::Type track_type,
																		 const vector<unsigned> &slice, GenomeTrackArrays::SliceFunctions slice_func, double slice_percentile,
																		 const Iterator_modifier1D &imdf1d, const Iterator_modifier2D &imdf2d)
{
	Iterator_modifier1D *pimdf1d = add_imdf(imdf1d);
	Iterator_modifier2D *pimdf2d = add_imdf(imdf2d);

	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		if (itrack_n_imdf->name == track &&
			itrack_n_imdf->slice == slice && itrack_n_imdf->slice_func == slice_func && itrack_n_imdf->slice_percentile == slice_percentile &&
			itrack_n_imdf->imdf1d == pimdf1d && itrack_n_imdf->imdf2d == pimdf2d)
			return *itrack_n_imdf;
	}

	if (m_track_n_imdfs.size() == m_track_n_imdfs.capacity())
		verror("Reached the limit of maximal number of tracks");

	m_track_n_imdfs.push_back(Track_n_imdf());
	Track_n_imdf &track_n_imdf = m_track_n_imdfs.back();
	track_n_imdf.name = track;
	track_n_imdf.track = NULL;
	track_n_imdf.type = track_type;
	track_n_imdf.slice = slice;
	track_n_imdf.slice_func = slice_func;
	track_n_imdf.slice_percentile = slice_percentile;
	track_n_imdf.imdf1d = pimdf1d;
	track_n_imdf.imdf2d = pimdf2d;
	return track_n_imdf;
}

void TrackExpressionVars::attach_filter_to_var(SEXP rvtrack, const string &vtrack, Track_var &var)
{
	// Check for filter field in rvtrack
	SEXP rfilter = get_rvector_col(rvtrack, "filter", vtrack.c_str(), false);

	if (Rf_isNull(rfilter) || !Rf_isString(rfilter) || Rf_length(rfilter) != 1) {
		var.filter = nullptr;
		return;
	}

	// Get filter key
	const char *filter_key = CHAR(STRING_ELT(rfilter, 0));

	// Look up filter in registry
	var.filter = FilterRegistry::instance().get(filter_key);

	if (var.filter == nullptr) {
		verror("Virtual track %s: filter with key '%s' not found in registry", vtrack.c_str(), filter_key);
	}
}

void TrackExpressionVars::attach_filter_to_var(SEXP rvtrack, const string &vtrack, Interv_var &var)
{
	// Check for filter field in rvtrack
	SEXP rfilter = get_rvector_col(rvtrack, "filter", vtrack.c_str(), false);

	if (Rf_isNull(rfilter) || !Rf_isString(rfilter) || Rf_length(rfilter) != 1) {
		var.filter = nullptr;
		return;
	}

	// Get filter key
	const char *filter_key = CHAR(STRING_ELT(rfilter, 0));

	// Look up filter in registry
	var.filter = FilterRegistry::instance().get(filter_key);

	if (var.filter == nullptr) {
		verror("Virtual track %s: filter with key '%s' not found in registry", vtrack.c_str(), filter_key);
	}
}

void TrackExpressionVars::attach_filter_to_var(SEXP rvtrack, const string &vtrack, Value_var &var)
{
	// Check for filter field in rvtrack
	SEXP rfilter = get_rvector_col(rvtrack, "filter", vtrack.c_str(), false);

	if (Rf_isNull(rfilter) || !Rf_isString(rfilter) || Rf_length(rfilter) != 1) {
		var.filter = nullptr;
		return;
	}

	// Get filter key
	const char *filter_key = CHAR(STRING_ELT(rfilter, 0));

	// Look up filter in registry
	var.filter = FilterRegistry::instance().get(filter_key);

	if (var.filter == nullptr) {
		verror("Virtual track %s: filter with key '%s' not found in registry", vtrack.c_str(), filter_key);
	}
}

void TrackExpressionVars::add_vtrack_var(const string &vtrack, SEXP rvtrack)
{
	// Check for existing track vars
	for (Track_vars::const_iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar)
	{
		if (ivar->var_name == vtrack)
			return;
	}

	// Check for existing interval vars
	for (Interv_vars::const_iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ++ivar)
	{
		if (ivar->var_name == vtrack)
			return;
	}

	// Check for existing value vars
	for (Value_vars::const_iterator ivar = m_value_vars.begin(); ivar != m_value_vars.end(); ++ivar)
	{
		if (ivar->var_name == vtrack)
			return;
	}

	// First check if this is a PWM function
	SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
    if (!Rf_isNull(rfunc) && Rf_isString(rfunc)) {
        string func = CHAR(STRING_ELT(rfunc, 0));
        transform(func.begin(), func.end(), func.begin(), ::tolower);
        
        if (func == "pwm" || func == "pwm.max" || func == "pwm.max.pos" || func == "pwm.count") {
            // Create the Track_var without a Track_n_imdf
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            var.val_func = (func == "pwm" ? Track_var::PWM :
                           func == "pwm.max" ? Track_var::PWM_MAX :
                           func == "pwm.max.pos" ? Track_var::PWM_MAX_POS :
                           Track_var::PWM_COUNT);
            var.track_n_imdf = nullptr;  // No track needed for PWM
            var.seq_imdf1d = nullptr;
            
            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);

			// Parse PWM parameters using helper struct
			TrackExprParams::PWMParams pwm_params = TrackExprParams::PWMParams::parse(rparams, vtrack);

			// Construct scorer with shared sequence fetcher for caching
			var.pwm_scorer = std::make_unique<PWMScorer>(
				pwm_params.core.pssm,
				&m_shared_seqfetch,
				pwm_params.extend_flag,
				func == "pwm" ? PWMScorer::TOTAL_LIKELIHOOD :
				func == "pwm.max" ? PWMScorer::MAX_LIKELIHOOD :
				func == "pwm.max.pos" ? PWMScorer::MAX_LIKELIHOOD_POS :
				PWMScorer::MOTIF_COUNT,
				static_cast<char>(pwm_params.core.strand_mode),
				pwm_params.core.spat_factor,
				pwm_params.core.spat_bin_size,
				static_cast<float>(pwm_params.core.score_thresh)
			);

            // Parse optional iterator modifier (sshift/eshift) for sequence-based vtracks
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);

            var.percentile = numeric_limits<double>::quiet_NaN();
            var.requires_pv = false;

            // Attach filter if present
            attach_filter_to_var(rvtrack, vtrack, var);
            return;
        } else if (func == "kmer.count" || func == "kmer.frac")	{
			// Create the Track_var without a Track_n_imdf
			m_track_vars.push_back(Track_var());
			Track_var &var = m_track_vars.back();
			var.var_name = vtrack;
			var.val_func = func == "kmer.count" ? Track_var::KMER_COUNT : Track_var::KMER_FRAC;
			var.track_n_imdf = nullptr; // No track needed for kmer
            var.seq_imdf1d = nullptr;
			var.percentile = numeric_limits<double>::quiet_NaN();
			SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);

			// Parse KMER parameters using helper struct
			TrackExprParams::KmerParams kmer_params = TrackExprParams::KmerParams::parse(rparams, vtrack);

			KmerCounter::CountMode mode = func == "kmer.count" ? KmerCounter::SUM : KmerCounter::FRACTION;
			var.kmer_counter = std::make_unique<KmerCounter>(kmer_params.kmer.c_str(), &m_shared_seqfetch, mode,
			                                                   kmer_params.extend, kmer_params.strand);

            // Parse optional iterator modifier (sshift/eshift) for sequence-based vtracks
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);

            var.requires_pv = false;

            // Attach filter if present
            attach_filter_to_var(rvtrack, vtrack, var);
			return;
		}
	}

	// Only check source if not a PWM function
	SEXP rsrc = get_rvector_col(rvtrack, "src", vtrack.c_str(), false);
	if (Rf_isNull(rsrc))
	{
		verror("Source must be specified for non-PWM virtual tracks");
	}

    if (Rf_isString(rsrc) && Rf_length(rsrc) == 1)
    {
        string track(CHAR(STRING_ELT(rsrc, 0)));

        SEXP gtracks = find_in_misha(m_iu.get_env(), "GTRACKS");
        if (Rf_isString(gtracks))
        {
            for (int itrack = 0; itrack < Rf_length(gtracks); itrack++)
            {
                if (!strcmp(CHAR(STRING_ELT(gtracks, itrack)), track.c_str()))
                {
                    add_vtrack_var_src_track(rvtrack, vtrack, track);
                    return;
                }
            }
        }
    }

	// Check if rsrc is a data.frame with interval columns and value column (value-based vtrack)
	if (TYPEOF(rsrc) == VECSXP) {  // data.frame is a VECSXP (list)
		SEXP names = Rf_getAttrib(rsrc, R_NamesSymbol);
		if (names != R_NilValue && TYPEOF(names) == STRSXP) {
			// Check if it has chrom, start, end columns
			bool has_chrom = false, has_start = false, has_end = false;
			int value_idx = -1;

			for (int i = 0; i < Rf_length(names); i++) {
				const char *name = CHAR(STRING_ELT(names, i));
				if (!strcmp(name, "chrom")) { has_chrom = true; }
				else if (!strcmp(name, "start")) { has_start = true; }
				else if (!strcmp(name, "end")) { has_end = true; }
				// Look for numeric columns (not chrom/start/end/strand) as value columns
				else if (strcmp(name, "strand") != 0 && strcmp(name, "intervalID") != 0) {
					SEXP col = VECTOR_ELT(rsrc, i);
					if (TYPEOF(col) == REALSXP || TYPEOF(col) == INTSXP) {
						if (value_idx == -1) {  // Take first numeric column as value
							value_idx = i;
						}
					}
				}
			}

			// Check what function is requested
			SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
			string func;
			if (Rf_isNull(rfunc))
				func = Value_var::FUNC_NAMES[Value_var::AVG];
			else {
				func = CHAR(STRING_ELT(rfunc, 0));
				transform(func.begin(), func.end(), func.begin(), ::tolower);
			}

			// Check if this is an interval-based function (distance, distance.center, coverage, neighbor.count)
			bool is_interval_func = (!strcmp(func.c_str(), "distance") ||
			                         !strcmp(func.c_str(), "distance.center") ||
			                         !strcmp(func.c_str(), "coverage") ||
			                         !strcmp(func.c_str(), "neighbor.count"));

			// If we have interval columns AND a value column AND it's NOT an interval function,
			// treat as value-based vtrack
			if (has_chrom && has_start && has_end && value_idx >= 0 && !is_interval_func) {
				GIntervals intervs;
				vector<float> vals;

				// Convert to intervals and extract values
				m_iu.convert_rintervs(rsrc, &intervs, NULL);

				// Extract values
				SEXP value_col = VECTOR_ELT(rsrc, value_idx);
				int n = Rf_length(value_col);
				vals.reserve(n);

				if (TYPEOF(value_col) == REALSXP) {
					for (int i = 0; i < n; i++) {
						vals.push_back(static_cast<float>(REAL(value_col)[i]));
					}
				} else if (TYPEOF(value_col) == INTSXP) {
					for (int i = 0; i < n; i++) {
						int ival = INTEGER(value_col)[i];
						if (ival == NA_INTEGER)
							vals.push_back(numeric_limits<float>::quiet_NaN());
						else
							vals.push_back(static_cast<float>(ival));
					}
				}

				add_vtrack_var_src_value(rvtrack, vtrack, intervs, vals);
				return;
			}
		}
	}

	GIntervals intervs1d;
	GIntervals2D intervs2d;

	try
	{
		m_iu.convert_rintervs(rsrc, &intervs1d, &intervs2d);
	} catch (TGLException &e) {
		verror("Source of virtual track %s was not recognized neither as a track nor as intervals. (\"%s\")", vtrack.c_str(), e.msg());
	}
	add_vtrack_var_src_interv(rvtrack, vtrack, intervs1d, intervs2d);
}

TrackExpressionVars::Track_var &TrackExpressionVars::add_track_var(const string &track)
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		if (ivar->var_name == track)
			return *ivar;
	}

	GenomeTrack::Type track_type = GenomeTrack::get_type(track2path(m_iu.get_env(), track).c_str(), m_iu.get_chromkey());
	Track_n_imdf &track_n_imdf = add_track_n_imdf(track, track_type, vector<unsigned>(), GenomeTrackArrays::S_AVG, 0, Iterator_modifier1D(), Iterator_modifier2D());

	m_track_vars.push_back(Track_var());
	Track_var &var = m_track_vars.back();
	var.var_name = track;
	var.val_func = Track_var::REG;
	var.percentile = numeric_limits<double>::quiet_NaN();
	var.requires_pv = false;
	var.track_n_imdf = &track_n_imdf;
	return var;
}

TrackExpressionVars::Track_var &TrackExpressionVars::add_vtrack_var_src_track(SEXP rvtrack, const string &vtrack, const string &track)
{
	GenomeTrack::Type track_type = GenomeTrack::get_type(track2path(m_iu.get_env(), track).c_str(), m_iu.get_chromkey());
	Iterator_modifier1D imdf1d;
	Iterator_modifier2D imdf2d;

	if (GenomeTrack::is_1d(track_type))
		parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
	else
		parse_imdf(rvtrack, vtrack, NULL, &imdf2d);

	SEXP rslice = get_rvector_col(rvtrack, "slice", vtrack.c_str(), false);
	vector<unsigned> slice;
	GenomeTrackArrays::SliceFunctions slice_func = GenomeTrackArrays::S_AVG;
	double slice_percentile = 0;

	if (!Rf_isNull(rslice)) {
		if (track_type != GenomeTrack::ARRAYS)
			verror("Slices are not supported by %s tracks", GenomeTrack::TYPE_NAMES[track_type]);

		SEXP rslice_func = get_rvector_col(rslice, "func", vtrack.c_str(), false);

    	if (!Rf_isNull(rslice_func)) {
    		string slice_func_str;

    		if (track_type != GenomeTrack::ARRAYS)
    			verror("Slices are not supported by %s tracks", GenomeTrack::TYPE_NAMES[track_type]);

    		if (!Rf_isString(rslice_func))
    			verror("slice function argument must be a string");

    		SEXP rslice_params = get_rvector_col(rslice, "params", vtrack.c_str(), false);

    		slice_func_str = CHAR(STRING_ELT(rslice_func, 0));
    		transform(slice_func_str.begin(), slice_func_str.end(), slice_func_str.begin(), ::tolower);

    		int ifunc;
    		for (ifunc = 0; ifunc < GenomeTrackArrays::NUM_S_FUNCS; ++ifunc) {
    			if (!strcmp(slice_func_str.c_str(), GenomeTrackArrays::SLICE_FUNCTION_NAMES[ifunc])) {
    				slice_func = (GenomeTrackArrays::SliceFunctions)ifunc;
					if (slice_func == GenomeTrackArrays::S_QUANTILE) {
						if (Rf_isNull(rslice_params))
							verror("Virtual track %s: slice function %s requires an additional parameter (percentile) to be specified", vtrack.c_str(), slice_func_str.c_str());
						if (!Rf_isReal(rslice_params) || Rf_length(rslice_params) != 1)
							verror("Virtual track %s: invalid parameters used for function %s", vtrack.c_str(), slice_func_str.c_str());
						slice_percentile = REAL(rslice_params)[0];
						if (slice_percentile < 0 || slice_percentile > 1)
							verror("Virtual track %s: parameter (percentile) used for function %s is out of range", vtrack.c_str(), slice_func_str.c_str());
					} else if (!Rf_isNull(rslice_params))
    					verror("Virtual track %s: slice function %s does not accept any parameters", vtrack.c_str(), slice_func_str.c_str());
    				break;
    			}
    		}

    		if (ifunc >= GenomeTrackArrays::NUM_S_FUNCS)
    			verror("Virtual track %s: invalid function %s used with a track", vtrack.c_str(), slice_func_str.c_str());
    	}

		SEXP rslice_idx = get_rvector_col(rslice, "slice", vtrack.c_str(), false);
		if (!Rf_isNull(rslice_idx)) {
			if (!Rf_isReal(rslice_idx) && !Rf_isInteger(rslice_idx)) 
				verror("Virtual track %s: invalid slice parameters", vtrack.c_str());

			for (int i = 0; i < Rf_length(rslice_idx); ++i) {
				double idx = Rf_isReal(rslice_idx) ? REAL(rslice_idx)[i] : INTEGER(rslice_idx)[i];
				if (idx < 1 || idx != (double)(int)idx)
					verror("Virtual track %s: slice indices must be positive integers", vtrack.c_str());
				slice.push_back((unsigned)(idx - 1));
			}

			sort(slice.begin(), slice.end());
			vector<unsigned>::iterator new_end = unique(slice.begin(), slice.end());
			slice.resize(new_end - slice.begin());
		}
	}

	Track_n_imdf &track_n_imdf = add_track_n_imdf(track, track_type, slice, slice_func, slice_percentile, imdf1d, imdf2d);

	m_track_vars.push_back(Track_var());
	Track_var &var = m_track_vars.back();
	var.var_name = vtrack;
	var.track_n_imdf = &track_n_imdf;

	SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
	SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);
	string func;

	if (Rf_isNull(rfunc))
		func = Track_var::FUNC_NAMES[Track_var::REG];
	else {
		if (!Rf_isString(rfunc))
			verror("Function argument must be a string");

		func = CHAR(STRING_ELT(rfunc, 0));
		transform(func.begin(), func.end(), func.begin(), ::tolower);
	}

	int ifunc;
	for (ifunc = 0; ifunc < Track_var::NUM_FUNCS; ++ifunc) {
		if (!strcmp(func.c_str(), Track_var::FUNC_NAMES[ifunc])) {
			if ((GenomeTrack::is_1d(track_type) && (ifunc == Track_var::WEIGHTED_SUM || ifunc == Track_var::OCCUPIED_AREA)) ||
					(GenomeTrack::is_2d(track_type) && (ifunc == Track_var::REG_NEAREST || ifunc == Track_var::STDDEV || ifunc == Track_var::SUM || ifunc == Track_var::QUANTILE)) ||
					(track_type != GenomeTrack::FIXED_BIN && (ifunc == Track_var::PV || ifunc == Track_var::PV_MIN || ifunc == Track_var::PV_MAX)))
				verror("Virtual track %s: function %s is not supported by %s tracks", vtrack.c_str(), func.c_str(), GenomeTrack::TYPE_NAMES[track_type]);
	
			if (ifunc == Track_var::QUANTILE) {
				if (Rf_isNull(rparams))
					verror("Virtual track %s: function %s requires an additional parameter (percentile) to be specified", vtrack.c_str(), func.c_str());
				if (!Rf_isReal(rparams) || Rf_length(rparams) != 1)
					verror("Virtual track %s: invalid parameters used for function %s", vtrack.c_str(), func.c_str());
				var.percentile = REAL(rparams)[0];
				if (var.percentile < 0 || var.percentile > 1)
					verror("Virtual track %s: parameter (percentile) used for function %s is out of range", vtrack.c_str(), func.c_str());
				
			} else	if (ifunc == Track_var::PWM || ifunc == Track_var::PWM_MAX || ifunc == Track_var::PWM_MAX_POS || ifunc == Track_var::PWM_COUNT) {
				var.percentile = numeric_limits<double>::quiet_NaN();

				TrackExprParams::PWMParams pwm_params = TrackExprParams::PWMParams::parse(rparams, vtrack);

				// Read extend parameter from vtrack (overrides default)
				SEXP rextend = get_rvector_col(rvtrack, "extend", vtrack.c_str(), false);
				bool extend = true; // default value (matching constructor default)
				if (!Rf_isNull(rextend)) {
					if (!Rf_isLogical(rextend))
						verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
					extend = LOGICAL(rextend)[0];
				}
				pwm_params.extend_flag = extend;
				pwm_params.core.extend = extend ? std::max(0, (int)pwm_params.core.pssm.size() - 1) : 0;

				var.pwm_scorer = std::make_unique<PWMScorer>(
					pwm_params.core.pssm,
					&m_shared_seqfetch,
					pwm_params.extend_flag,
					ifunc == Track_var::PWM ? PWMScorer::TOTAL_LIKELIHOOD :
					ifunc == Track_var::PWM_MAX ? PWMScorer::MAX_LIKELIHOOD :
					ifunc == Track_var::PWM_MAX_POS ? PWMScorer::MAX_LIKELIHOOD_POS :
					PWMScorer::MOTIF_COUNT,
					static_cast<char>(pwm_params.core.strand_mode),
					pwm_params.core.spat_factor,
					pwm_params.core.spat_bin_size,
					static_cast<float>(pwm_params.core.score_thresh));
			} else if(ifunc == Track_var::KMER_COUNT || ifunc == Track_var::KMER_FRAC) {
				if (Rf_isNull(rparams)){
					verror("Virtual track %s: function %s requires a parameter (kmer string)", vtrack.c_str(), func.c_str());
				}

				if (!Rf_isString(rparams) || Rf_length(rparams) != 1){
					verror("Virtual track %s: function %s requires a string parameter", vtrack.c_str(), func.c_str());
				}
				var.percentile = numeric_limits<double>::quiet_NaN();

				// Read extend parameter from vtrack
				SEXP rextend = get_rvector_col(rvtrack, "extend", vtrack.c_str(), false);
				bool extend = true; // default value (matching constructor default)
				if (!Rf_isNull(rextend)) {
					if (!Rf_isLogical(rextend))
						verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
					extend = LOGICAL(rextend)[0];
				}

				// Read strand parameter from vtrack
				SEXP rstrand = get_rvector_col(rvtrack, "strand", vtrack.c_str(), false);
				char strand = 0; // default value (both strands)
				if (!Rf_isNull(rstrand)) {
					if (!Rf_isReal(rstrand) || Rf_length(rstrand) != 1)
						verror("Virtual track %s: strand parameter must be numeric", vtrack.c_str());
					strand = (char)REAL(rstrand)[0];
				}

				// Get the kmer string
				const char *kmer = CHAR(STRING_ELT(rparams, 0));
				KmerCounter::CountMode mode = ifunc == Track_var::KMER_COUNT ? KmerCounter::SUM : KmerCounter::FRACTION;

				var.kmer_counter = std::make_unique<KmerCounter>(kmer, &m_shared_seqfetch, mode, extend, strand);
				break;

			} else{
				var.percentile = numeric_limits<double>::quiet_NaN();
				if (!Rf_isNull(rparams))
					verror("Virtual track %s: function %s does not accept any parameters", vtrack.c_str(), func.c_str());
			}

			var.val_func = (Track_var::Val_func)ifunc;
			var.requires_pv = ifunc == Track_var::PV || ifunc == Track_var::PV_MIN || ifunc == Track_var::PV_MAX;
			break;
		}
	}

	if (ifunc >= Track_var::NUM_FUNCS)
        verror("Virtual track %s: invalid function %s used for a virtual track", vtrack.c_str(), func.c_str());

	// Attach filter if present
	// Check if filter is being applied to 2D track
	SEXP rfilter = get_rvector_col(rvtrack, "filter", vtrack.c_str(), false);
	if (!Rf_isNull(rfilter) && GenomeTrack::is_2d(track_type)) {
		verror("Virtual track %s: Filters are not yet supported for 2D tracks", vtrack.c_str());
	}

	attach_filter_to_var(rvtrack, vtrack, var);

	return var;
}

TrackExpressionVars::Interv_var &TrackExpressionVars::add_vtrack_var_src_interv(SEXP rvtrack, const string &vtrack, GIntervals &intervs1d, GIntervals2D &intervs2d)
{
	if (intervs1d.empty() && !intervs2d.empty())
		verror("Virtual track %s: virtual tracks do not support 2D intervals as a source", vtrack.c_str());

	m_interv_vars.push_back(Interv_var());
	Interv_var &var = m_interv_vars.back();

	var.var_name = vtrack;
	Iterator_modifier1D imdf1d;
	parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
	var.imdf1d = add_imdf(imdf1d);

	SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
	SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);
	string func;

	if (Rf_isNull(rfunc))
		func = Interv_var::FUNC_NAMES[Interv_var::DIST];
	else {
		func = CHAR(STRING_ELT(rfunc, 0));
		transform(func.begin(), func.end(), func.begin(), ::tolower);
	}

	if (!strcmp(func.c_str(), Interv_var::FUNC_NAMES[Interv_var::DIST])) {
		var.val_func = Interv_var::DIST;

		double dist_margin = 0;

		if (!Rf_isNull(rparams)) {
			if (Rf_isReal(rparams) && Rf_length(rparams) == 1)
				dist_margin = REAL(rparams)[0];
			else if (Rf_isInteger(rparams) && Rf_length(rparams) == 1)
				dist_margin = INTEGER(rparams)[0];
			else
				verror("Virtual track %s: invalid parameters used for function %s", vtrack.c_str(), func.c_str());
		}
		if (dist_margin < 0)
			verror("Virtual track %s, function %s: margin cannot be a negative number", vtrack.c_str(), func.c_str());

		var.dist_margin = dist_margin * 0.5;
		var.sintervs.swap(intervs1d);
		var.sintervs.sort();
		var.eintervs = var.sintervs;
		var.eintervs.sort(GIntervals::compare_by_end_coord);
		// cannot not set var.siinterv and var.eiinterv now because these iterators will be invalidated one a new var is added to m_interv_vars

		if (var.dist_margin) {
			// if dist_margin is not zero => we are measuring distances from the centers of the interval, therefore the intervals cannot overlap
			for (GIntervals::const_iterator iinterv = var.sintervs.begin() + 1; iinterv < var.sintervs.end(); ++iinterv) {
				if (iinterv->do_touch(*(iinterv - 1)))
					verror("Virtual track %s: intervals are overlapping and hence incompatible with %s function having non-zero (%g) margin", vtrack.c_str(), func.c_str(), var.dist_margin);
			}
		}
	} else if (!strcmp(func.c_str(), Interv_var::FUNC_NAMES[Interv_var::DIST_CENTER])) {
		var.val_func = Interv_var::DIST_CENTER;

		if (!Rf_isNull(rparams))
			verror("Virtual track %s: function %s does not accept any parameters", vtrack.c_str(), func.c_str());

		var.dist_margin = 0.;
		var.sintervs.swap(intervs1d);
		var.sintervs.sort();
		// cannot not set var.siinterv and var.eiinterv now because these iterators will be invalidated one a new var is added to m_interv_vars

		for (GIntervals::const_iterator iinterv = var.sintervs.begin() + 1; iinterv < var.sintervs.end(); ++iinterv) {
			if (iinterv->do_touch(*(iinterv - 1)))
				verror("Virtual track %s: intervals are overlapping and hence incompatible with %s function", vtrack.c_str(), func.c_str());
		}
	} else if (!strcmp(func.c_str(), Interv_var::FUNC_NAMES[Interv_var::COVERAGE])) {
        var.val_func = Interv_var::COVERAGE;

        if (!Rf_isNull(rparams))
            verror("Virtual track %s: function %s does not accept any parameters", vtrack.c_str(), func.c_str());

        var.dist_margin = 0.;
        var.sintervs.swap(intervs1d);
        var.sintervs.sort();
        var.sintervs.unify_overlaps(); // Unify overlaps since we want total coverage
	
	} else if (!strcmp(func.c_str(), Interv_var::FUNC_NAMES[Interv_var::NEIGHBOR_COUNT])) {
		var.val_func = Interv_var::NEIGHBOR_COUNT;

		double dist_margin = 0;

		if (!Rf_isNull(rparams)) {
			if (Rf_isReal(rparams) && Rf_length(rparams) == 1)
				dist_margin = REAL(rparams)[0];
			else if (Rf_isInteger(rparams) && Rf_length(rparams) == 1)
				dist_margin = INTEGER(rparams)[0];
			else
				verror("Virtual track %s: invalid parameters used for function %s", vtrack.c_str(), func.c_str());
		}
		if (dist_margin < 0)
			verror("Virtual track %s, function %s: distance cannot be a negative number", vtrack.c_str(), func.c_str());

		var.dist_margin = dist_margin;
		var.sintervs.swap(intervs1d);
		var.sintervs.sort();

		const GenomeChromKey &chromkey = m_iu.get_chromkey();

		// Create expanded intervals but do NOT unify them
		var.eintervs.clear();
		var.eintervs.reserve(var.sintervs.size());
		for (const auto &interv : var.sintervs) {
			double expanded_start = static_cast<double>(interv.start) - var.dist_margin;
			double expanded_end = static_cast<double>(interv.end) + var.dist_margin;

			int64_t new_start = expanded_start < 0.0 ? 0 : static_cast<int64_t>(std::floor(expanded_start));
			uint64_t chrom_size = chromkey.get_chrom_size(interv.chromid);
			double chrom_size_d = static_cast<double>(chrom_size);
			if (expanded_end > chrom_size_d)
				expanded_end = chrom_size_d;

			int64_t new_end = static_cast<int64_t>(std::ceil(expanded_end));
			if (new_end < new_start)
				new_end = new_start;

			var.eintervs.push_back(GInterval(interv.chromid, new_start, new_end, interv.strand));
		}
		var.eintervs.sort();
	} else {
		verror("Virtual track %s: invalid function %s used with intervals", vtrack.c_str(), func.c_str());
	}

	// Attach filter if present
	attach_filter_to_var(rvtrack, vtrack, var);

	return var;
}

TrackExpressionVars::Value_var &TrackExpressionVars::add_vtrack_var_src_value(SEXP rvtrack, const string &vtrack, GIntervals &intervs, vector<float> &vals)
{
	if (intervs.size() != vals.size())
		verror("Virtual track %s: number of intervals (%zu) does not match number of values (%zu)", vtrack.c_str(), intervs.size(), vals.size());

	m_value_vars.push_back(Value_var());
	Value_var &var = m_value_vars.back();

	var.var_name = vtrack;
	Iterator_modifier1D imdf1d;
	parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
	var.imdf1d = add_imdf(imdf1d);

	SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
	SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);
	string func;

	if (Rf_isNull(rfunc))
		func = Value_var::FUNC_NAMES[Value_var::AVG];
	else {
		func = CHAR(STRING_ELT(rfunc, 0));
		transform(func.begin(), func.end(), func.begin(), ::tolower);
	}

	// Parse function name
	bool func_found = false;
	for (int i = 0; i < Value_var::NUM_FUNCS; i++) {
		if (!strcmp(func.c_str(), Value_var::FUNC_NAMES[i])) {
			var.val_func = static_cast<Value_var::Val_func>(i);
			func_found = true;
			break;
		}
	}

	if (!func_found)
		verror("Virtual track %s: invalid function %s used with value-based intervals", vtrack.c_str(), func.c_str());

	// Parse percentile parameter for quantile function
	var.percentile = 0;
	if (var.val_func == Value_var::QUANTILE) {
		if (Rf_isNull(rparams))
			verror("Virtual track %s: function %s requires percentile parameter", vtrack.c_str(), func.c_str());

		if (Rf_isReal(rparams) && Rf_length(rparams) == 1)
			var.percentile = REAL(rparams)[0];
		else if (Rf_isInteger(rparams) && Rf_length(rparams) == 1)
			var.percentile = INTEGER(rparams)[0];
		else
			verror("Virtual track %s: invalid percentile parameter for function %s", vtrack.c_str(), func.c_str());

		if (var.percentile < 0 || var.percentile > 1)
			verror("Virtual track %s: percentile must be between 0 and 1, got %g", vtrack.c_str(), var.percentile);
	} else {
		if (!Rf_isNull(rparams))
			verror("Virtual track %s: function %s does not accept any parameters", vtrack.c_str(), func.c_str());
	}

	// Sort intervals and values together
	// Create pairs, sort, then separate
	vector<pair<GInterval, float>> paired;
	paired.reserve(intervs.size());
	for (size_t i = 0; i < intervs.size(); i++) {
		paired.push_back(make_pair(intervs[i], vals[i]));
	}

	sort(paired.begin(), paired.end(), [](const pair<GInterval, float> &a, const pair<GInterval, float> &b) {
		return a.first < b.first;
	});

	// Separate back into sorted intervals and values
	GIntervals sorted_intervals;
	vector<float> sorted_values;
	sorted_intervals.reserve(paired.size());
	sorted_values.reserve(paired.size());
	for (size_t i = 0; i < paired.size(); i++) {
		sorted_intervals.push_back(paired[i].first);
		sorted_values.push_back(paired[i].second);
	}

	// Check for overlaps (Phase 1: error only)
	for (size_t i = 1; i < sorted_intervals.size(); i++) {
		if (sorted_intervals[i-1].chromid == sorted_intervals[i].chromid &&
		    sorted_intervals[i-1].end > sorted_intervals[i].start) {
			verror("Virtual track %s: overlapping intervals detected. Interval at position %zu [%s:%ld-%ld] overlaps with previous interval [%s:%ld-%ld]. "
			       "Value-based virtual tracks do not support overlapping intervals.",
			       vtrack.c_str(), i,
			       m_iu.get_chromkey().id2chrom(sorted_intervals[i].chromid).c_str(),
			       sorted_intervals[i].start, sorted_intervals[i].end,
			       m_iu.get_chromkey().id2chrom(sorted_intervals[i-1].chromid).c_str(),
			       sorted_intervals[i-1].start, sorted_intervals[i-1].end);
		}
	}

	// Create the track object
	var.track = make_shared<GenomeTrackInMemory>();

	// Determine chromid (assuming all intervals are on same chrom, or 0 if empty)
	int chromid = sorted_intervals.empty() ? 0 : sorted_intervals.front().chromid;

	// Initialize track with data
	var.track->init_from_data(sorted_intervals, sorted_values, chromid);

	// Register required function
	switch (var.val_func) {
		case Value_var::AVG:
			var.track->register_function(GenomeTrack1D::AVG);
			var.track->register_function(GenomeTrack1D::SIZE);
			break;
		case Value_var::MIN:
			var.track->register_function(GenomeTrack1D::MIN);
			break;
		case Value_var::MAX:
			var.track->register_function(GenomeTrack1D::MAX);
			break;
		case Value_var::SUM:
			var.track->register_function(GenomeTrack1D::SUM);
			break;
		case Value_var::STDDEV:
			var.track->register_function(GenomeTrack1D::STDDEV);
			var.track->register_function(GenomeTrack1D::SIZE);
			break;
		case Value_var::QUANTILE:
			var.track->register_quantile(10000, 1000, 1000);
			break;
		case Value_var::NEAREST:
			var.track->register_function(GenomeTrack1D::NEAREST);
			break;
		case Value_var::EXISTS:
			var.track->register_function(GenomeTrack1D::EXISTS);
			break;
		case Value_var::SIZE:
			var.track->register_function(GenomeTrack1D::SIZE);
			break;
		case Value_var::FIRST:
			var.track->register_function(GenomeTrack1D::FIRST);
			break;
		case Value_var::LAST:
			var.track->register_function(GenomeTrack1D::LAST);
			break;
		case Value_var::SAMPLE:
			var.track->register_function(GenomeTrack1D::SAMPLE);
			break;
		case Value_var::FIRST_POS_ABS:
		case Value_var::FIRST_POS_REL:
			var.track->register_function(GenomeTrack1D::FIRST_POS);
			break;
		case Value_var::LAST_POS_ABS:
		case Value_var::LAST_POS_REL:
			var.track->register_function(GenomeTrack1D::LAST_POS);
			break;
		case Value_var::SAMPLE_POS_ABS:
		case Value_var::SAMPLE_POS_REL:
			var.track->register_function(GenomeTrack1D::SAMPLE_POS);
			break;
		case Value_var::MIN_POS_ABS:
		case Value_var::MIN_POS_REL:
			var.track->register_function(GenomeTrack1D::MIN_POS);
			break;
		case Value_var::MAX_POS_ABS:
		case Value_var::MAX_POS_REL:
			var.track->register_function(GenomeTrack1D::MAX_POS);
			break;
		case Value_var::NUM_FUNCS:
		default:
			// NUM_FUNCS is a sentinel value, not a valid function type
			break;
	}

	// Attach filter if present
	attach_filter_to_var(rvtrack, vtrack, var);

	return var;
}

void TrackExpressionVars::register_track_functions()
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		// Skip sequence-based variables since they don't have associated tracks
        if (TrackExpressionVars::is_sequence_based_function(ivar->val_func)) {
            continue;
        }
		GenomeTrack1D *track1d = GenomeTrack::is_1d(ivar->track_n_imdf->type) ? (GenomeTrack1D *)ivar->track_n_imdf->track : NULL;
		GenomeTrack2D *track2d = GenomeTrack::is_2d(ivar->track_n_imdf->type) ? (GenomeTrack2D *)ivar->track_n_imdf->track : NULL;

		switch (ivar->val_func) {
		case Track_var::REG:
		case Track_var::PV:
			if (track1d)
				track1d->register_function(GenomeTrack1D::AVG);
			else
				track2d->register_function(GenomeTrack2D::AVG);
			break;
		case Track_var::REG_MIN:
		case Track_var::PV_MIN:
			if (track1d)
				track1d->register_function(GenomeTrack1D::MIN);
			else
				track2d->register_function(GenomeTrack2D::MIN);
			break;
		case Track_var::REG_MAX:
		case Track_var::PV_MAX:
			if (track1d)
				track1d->register_function(GenomeTrack1D::MAX);
			else
				track2d->register_function(GenomeTrack2D::MAX);
			break;
		case Track_var::MAX_POS_ABS:
		case Track_var::MAX_POS_REL:
			if (!track1d)
				verror("vtrack functions 'max.pos.abs' and 'max.pos.relative' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::MAX_POS);
			break;
		case Track_var::MIN_POS_ABS:
		case Track_var::MIN_POS_REL:
			if (!track1d)
				verror("vtrack functions 'min.pos.abs' and 'min.pos.relative' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::MIN_POS);
			break;
		case Track_var::REG_NEAREST:
			track1d->register_function(GenomeTrack1D::NEAREST);
			break;
		case Track_var::STDDEV:
			track1d->register_function(GenomeTrack1D::STDDEV);
			break;
		case Track_var::SUM:
			track1d->register_function(GenomeTrack1D::SUM);
			break;
		case Track_var::QUANTILE:
			track1d->register_quantile(m_iu.get_max_data_size(), m_iu.get_quantile_edge_data_size(), m_iu.get_quantile_edge_data_size());
			break;
		case Track_var::WEIGHTED_SUM:
			track2d->register_function(GenomeTrack2D::WEIGHTED_SUM);
			break;
		case Track_var::OCCUPIED_AREA:
			track2d->register_function(GenomeTrack2D::OCCUPIED_AREA);
			break;
		case Track_var::EXISTS:
			if (!track1d)
				verror("vtrack function 'exists' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::EXISTS);
			break;
		case Track_var::SIZE:
			if (!track1d)
				verror("vtrack function 'size' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::SIZE);
			break;
		case Track_var::SAMPLE:
			if (!track1d)
				verror("vtrack function 'sample' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::SAMPLE);
			break;
		case Track_var::SAMPLE_POS_ABS:
		case Track_var::SAMPLE_POS_REL:
			if (!track1d)
				verror("vtrack functions 'sample.pos.abs' and 'sample.pos.relative' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::SAMPLE_POS);
			break;
		case Track_var::FIRST:
			if (!track1d)
				verror("vtrack function 'first' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::FIRST);
			break;
		case Track_var::FIRST_POS_ABS:
		case Track_var::FIRST_POS_REL:
			if (!track1d)
				verror("vtrack functions 'first.pos.abs' and 'first.pos.relative' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::FIRST_POS);
			break;
		case Track_var::LAST:
			if (!track1d)
				verror("vtrack function 'last' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::LAST);
			break;
		case Track_var::LAST_POS_ABS:
		case Track_var::LAST_POS_REL:
			if (!track1d)
				verror("vtrack functions 'last.pos.abs' and 'last.pos.relative' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::LAST_POS);
			break;
		// Sequence-based functions work directly on sequences, no need to register track functions
		default:
			if (!TrackExpressionVars::is_sequence_based_function((Track_var::Val_func)ivar->val_func))
				verror("Unrecognized virtual track function");
			break;
		}

		if (ivar->track_n_imdf->type == GenomeTrack::ARRAYS) {
			GenomeTrackArrays *track = (GenomeTrackArrays *)ivar->track_n_imdf->track;
			if (ivar->track_n_imdf->slice_func == GenomeTrackArrays::S_QUANTILE) 
				track->set_slice_quantile(ivar->track_n_imdf->slice_percentile, m_iu.get_max_data_size(), m_iu.get_quantile_edge_data_size(), m_iu.get_quantile_edge_data_size(),
										  ivar->track_n_imdf->slice);
			else 
				track->set_slice_function(ivar->track_n_imdf->slice_func, ivar->track_n_imdf->slice);
		}
	}
}

void TrackExpressionVars::init(const TrackExpressionIteratorBase &expr_itr)
{
	// First validate iterator compatibility
	for (Track_vars::const_iterator itrack_var = m_track_vars.begin(); itrack_var != m_track_vars.end(); ++itrack_var)
	{
	    // Skip iterator validation for sequence-based variables since they don't have tracks or imdf
        if (TrackExpressionVars::is_sequence_based_function(itrack_var->val_func)) {
            continue;
        }

		if (expr_itr.is_1d())
		{
			if (GenomeTrack::is_1d(itrack_var->track_n_imdf->type) && itrack_var->track_n_imdf->imdf1d && itrack_var->track_n_imdf->imdf1d->dim != Iterator_modifier1D::DIM_NONE)
				verror("Virtual track %s: 1D iterator is used for a virtual track that specifies dimension projection and hence expects 2D iterators", itrack_var->var_name.c_str());
			if (GenomeTrack::is_2d(itrack_var->track_n_imdf->type))
			{
				if (itrack_var->var_name == itrack_var->track_n_imdf->name)
					verror("1D iterator is applied to a 2D track %s", itrack_var->track_n_imdf->name.c_str());
				verror("Virtual track %s: 1D iterator is applied to a 2D track %s", itrack_var->var_name.c_str(), itrack_var->track_n_imdf->name.c_str());
			}
		}
		else if (expr_itr.is_2d() && GenomeTrack::is_1d(itrack_var->track_n_imdf->type) && (!itrack_var->track_n_imdf->imdf1d || itrack_var->track_n_imdf->imdf1d->dim == Iterator_modifier1D::DIM_NONE))
		{
			if (itrack_var->var_name == itrack_var->track_n_imdf->name)
				verror("2D iterator is applied to a 1D track %s without explicit dimension projection", itrack_var->track_n_imdf->name.c_str());
			verror("Virtual track %s: 2D iterator is applied to a 1D track %s without explicit dimension projection", itrack_var->var_name.c_str(), itrack_var->track_n_imdf->name.c_str());
		}
	}

	for (Interv_vars::const_iterator iinterv_var = m_interv_vars.begin(); iinterv_var != m_interv_vars.end(); ++iinterv_var)
	{
		if (expr_itr.is_1d() && iinterv_var->imdf1d && iinterv_var->imdf1d->dim != Iterator_modifier1D::DIM_NONE)
			verror("Virtual track %s: 1D iterator is used for a virtual track that specifies dimension projection and hence expects 2D iterators", iinterv_var->var_name.c_str());
		else if (expr_itr.is_2d() && (!iinterv_var->imdf1d || iinterv_var->imdf1d->dim == Iterator_modifier1D::DIM_NONE))
			verror("Virtual track %s: 2D iterator is used without explicit dimension projection", iinterv_var->var_name.c_str());
	}

	// Collect non-PWM tracks that need file checking
	vector<string> track_names;
	vector<GenomeTrack::Type> track_types;

	for (Track_vars::const_iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar)
	{
		// Skip sequence-based tracks
		if (TrackExpressionVars::is_sequence_based_function(ivar->val_func))
		{
			continue;
		}
		track_names.push_back(ivar->track_n_imdf->name);
		track_types.push_back(ivar->track_n_imdf->type);
	}

	// Do track file validation only for non-PWM tracks
	if (!track_names.empty())
	{
		try
		{
			for (vector<string>::const_iterator itrack_name = track_names.begin(); itrack_name != track_names.end(); ++itrack_name)
			{
				string trackpath(track2path(m_iu.get_env(), *itrack_name));
				GenomeTrack::Type track_type = track_types[itrack_name - track_names.begin()];
				vector<string> filenames;
				unsigned binsize = 0;

				// read the list of chrom files
				get_chrom_files(trackpath.c_str(), filenames);
				sort(filenames.begin(), filenames.end());

				if (GenomeTrack::is_1d(track_type))
				{
					set<int> chromids;

					for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename)
					{
						int chromid = -1;
						GenomeTrackFixedBin gtrack_fbin;

						try
						{
							chromid = GenomeTrack::get_chromid_1d(m_iu.get_chromkey(), *ifilename);
						}
						catch (TGLException &e)
						{
							verror("Track %s: %s\n", itrack_name->c_str(), e.msg());
						}

						chromids.insert(chromid);
						if (track_type == GenomeTrack::FIXED_BIN)
						{
							gtrack_fbin.init_read((trackpath + "/" + *ifilename).c_str(), chromid);
							if (ifilename == filenames.begin())
							{
								binsize = gtrack_fbin.get_bin_size();
							}
							else if (binsize != gtrack_fbin.get_bin_size())
								verror("Track %s: bin size of chroms %s and %s differ (%d and %d respectively)",
									   itrack_name->c_str(), m_iu.id2chrom(GenomeTrack::get_chromid_1d(m_iu.get_chromkey(), *(ifilename - 1))).c_str(),
									   m_iu.id2chrom(chromid).c_str(), binsize, gtrack_fbin.get_bin_size());
						}
					}
				}
				else if (GenomeTrack::is_2d(track_type))
				{
					for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename)
					{
						try
						{
							GenomeTrack::get_chromid_2d(m_iu.get_chromkey(), *ifilename);
						}
						catch (TGLException &e)
						{
							verror("Track %s: %s\n", itrack_name->c_str(), e.msg());
						}
					}
				}
			}
		}
		catch (TGLException &e)
		{
			verror("%s\n", e.msg());
		}
	}

	// Handle percentile variables initialization
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar)
	{
		if (ivar->requires_pv)
		{
			// load binned_pvals
			SEXP val;
			Binned_pv &pv = ivar->pv_binned;
			string pv_fname(track2path(m_iu.get_env(), ivar->track_n_imdf->name) + "/vars/pv.percentiles");

			if (access(pv_fname.c_str(), R_OK) < 0 && errno == ENOENT)
			{
                char command[1000];

				REprintf("Preparing track %s for percentiles queries\n", ivar->track_n_imdf->name.c_str());
				snprintf(command, sizeof(command),
						 "{ "
						 "    .ginteractive = getOption(\".ginteractive\")\n"
						 "    tryCatch({\n"
						 "            options(.ginteractive = FALSE)\n"
						 "            misha:::.gtrack.prepare.pvals(\"%s\")\n"
						 "        },\n"
						 "        finally = { options(.ginteractive = .ginteractive) })"
						 " }",
						 ivar->track_n_imdf->name.c_str());
                run_in_R(command, m_iu.get_env());
			}

			rprotect(val = RSaneUnserialize(pv_fname.c_str()));
			SEXPCleaner val_cleaner(val);
			SEXP breaks = Rf_getAttrib(val, Rf_install("breaks"));

			if (breaks == R_NilValue || !Rf_isReal(breaks) || Rf_length(breaks) != Rf_length(val))
				verror("File %s is in invalid format.", pv_fname.c_str());

			pv.bins.assign(REAL(val), REAL(val) + Rf_length(val));
			pv.binfinder.init(REAL(breaks), Rf_length(breaks));
		}
	}
}

void TrackExpressionVars::define_r_vars(unsigned size)
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ivar++) {
		rprotect(ivar->rvar = RSaneAllocVector(REALSXP, size));
		Rf_defineVar(Rf_install(ivar->var_name.c_str()), ivar->rvar, m_iu.get_env());
		ivar->var = REAL(ivar->rvar);
	}
	for (Interv_vars::iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ivar++) {
		rprotect(ivar->rvar = RSaneAllocVector(REALSXP, size));
		Rf_defineVar(Rf_install(ivar->var_name.c_str()), ivar->rvar, m_iu.get_env());
		ivar->var = REAL(ivar->rvar);
	}
	for (Value_vars::iterator ivar = m_value_vars.begin(); ivar != m_value_vars.end(); ivar++) {
		rprotect(ivar->rvar = RSaneAllocVector(REALSXP, size));
		Rf_defineVar(Rf_install(ivar->var_name.c_str()), ivar->rvar, m_iu.get_env());
		ivar->var = REAL(ivar->rvar);
	}
}

void TrackExpressionVars::start_chrom(const GInterval &interval)
{
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		// Skip track initialization for sequence-based tracks
		bool is_sequence_track = false;
		for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar)
		{
			if (ivar->track_n_imdf == &(*itrack_n_imdf) &&
				TrackExpressionVars::is_sequence_based_function(ivar->val_func))
			{
				is_sequence_track = true;
				break;
			}
		}
		if (is_sequence_track)
			continue;

		try
		{
			string filename(track2path(m_iu.get_env(), itrack_n_imdf->name) + "/" + GenomeTrack::get_1d_filename(m_iu.get_chromkey(), interval.chromid));

			delete itrack_n_imdf->track;
			if (itrack_n_imdf->type == GenomeTrack::FIXED_BIN)
			{
				itrack_n_imdf->track = new GenomeTrackFixedBin;
				((GenomeTrackFixedBin *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid);
			}
			else if (itrack_n_imdf->type == GenomeTrack::SPARSE)
			{
				itrack_n_imdf->track = new GenomeTrackSparse;
				((GenomeTrackSparse *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid);
			}
			else if (itrack_n_imdf->type == GenomeTrack::ARRAYS)
			{
				itrack_n_imdf->track = new GenomeTrackArrays;
				((GenomeTrackArrays *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid);
			}
			else
			{
				verror("Internal error: track %s of type %s is not supported by 1D iterators", itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
			}
		}
		catch (TGLException &e)
		{
			verror("%s\n", e.msg());
		}
	}
	register_track_functions();

	// Invalidate PWM scorer caches on chromosome change
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		if (ivar->pwm_scorer) {
			ivar->pwm_scorer->invalidate_cache();
		}
	}
}

void TrackExpressionVars::start_chrom(const GInterval2D &interval)
{
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf) {
		try {
			if (GenomeTrack::is_1d(itrack_n_imdf->type)) {
				int chromid = 0;

				if (itrack_n_imdf->imdf1d->dim == Iterator_modifier1D::DIM1)
					chromid = interval.chromid1();
				else if (itrack_n_imdf->imdf1d->dim == Iterator_modifier1D::DIM2)
					chromid = interval.chromid2();
				else
					verror("Internal error: no 2D to 1D conversion for track %s", itrack_n_imdf->name.c_str());

				if (chromid != itrack_n_imdf->imdf1d->interval.chromid) {
					string filename(track2path(m_iu.get_env(), itrack_n_imdf->name) + "/" + GenomeTrack::get_1d_filename(m_iu.get_chromkey(), chromid));

					delete itrack_n_imdf->track;
					if (itrack_n_imdf->type == GenomeTrack::FIXED_BIN) {
						itrack_n_imdf->track = new GenomeTrackFixedBin;
						((GenomeTrackFixedBin *)itrack_n_imdf->track)->init_read(filename.c_str(), chromid);
					} else if (itrack_n_imdf->type == GenomeTrack::SPARSE) {
						itrack_n_imdf->track = new GenomeTrackSparse;
						((GenomeTrackSparse *)itrack_n_imdf->track)->init_read(filename.c_str(), chromid);
					} else if (itrack_n_imdf->type == GenomeTrack::ARRAYS) {
						itrack_n_imdf->track = new GenomeTrackArrays;
						((GenomeTrackArrays *)itrack_n_imdf->track)->init_read(filename.c_str(), chromid);
					} else
						verror("Internal error: track %s of type %s is not supported by 1D iterators (projected from 2D)", itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
				}
			} else if (!m_interval2d.is_same_chrom(interval)) {
				string filename(track2path(m_iu.get_env(), itrack_n_imdf->name) + "/" + GenomeTrack::get_2d_filename(m_iu.get_chromkey(), interval.chromid1(), interval.chromid2()));

				delete itrack_n_imdf->track;
				if (itrack_n_imdf->type == GenomeTrack::RECTS) {
					itrack_n_imdf->track = new GenomeTrackRectsRects(m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					((GenomeTrackRectsRects *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
				} else if (itrack_n_imdf->type == GenomeTrack::POINTS) {
					itrack_n_imdf->track = new GenomeTrackRectsPoints(m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					((GenomeTrackRectsPoints *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
				} else if (itrack_n_imdf->type == GenomeTrack::COMPUTED) {
					itrack_n_imdf->track = new GenomeTrackComputed(m_groot, m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					((GenomeTrackComputed *)itrack_n_imdf->track)->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
				} else
					verror("Internal error: track %s of type %s is not supported by 2D iterators", itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
			}
		} catch (TGLException &e) {
			verror("%s\n", e.msg());
		}
	}
	register_track_functions();

	// Invalidate PWM scorer caches on chromosome change
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		if (ivar->pwm_scorer) {
			ivar->pwm_scorer->invalidate_cache();
		}
	}
}

void TrackExpressionVars::set_vars(const GInterval &interval, unsigned idx)
{
    if (m_interval1d.chromid != interval.chromid) {
        start_chrom(interval);
	}

	m_interval1d = interval;

    for (Iterator_modifiers1D::iterator iimdf = m_imdfs1d.begin(); iimdf != m_imdfs1d.end(); ++iimdf)
        iimdf->transform(interval, m_iu.get_chromkey());

    set_vars(idx);
}

void TrackExpressionVars::set_vars(const GInterval2D &interval, const DiagonalBand &band, unsigned idx)
{
	if (!m_interval2d.is_same_chrom(interval))
		start_chrom(interval);

	m_interval2d = interval;
	m_band = band;

    for (Iterator_modifiers1D::iterator iimdf = m_imdfs1d.begin(); iimdf != m_imdfs1d.end(); ++iimdf) {
        iimdf->transform(interval, m_iu.get_chromkey());
    }

    for (Iterator_modifiers2D::iterator iimdf = m_imdfs2d.begin(); iimdf != m_imdfs2d.end(); ++iimdf) {
        iimdf->transform(interval, m_iu.get_chromkey());
    }

	set_vars(idx);
}

void TrackExpressionVars::set_vars(unsigned idx)
{
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		// Skip track setup for sequence-based tracks
		Track_vars::iterator seq_var = m_track_vars.begin();
		for (; seq_var != m_track_vars.end(); ++seq_var)
		{
			if (seq_var->track_n_imdf == &(*itrack_n_imdf) &&
				TrackExpressionVars::is_sequence_based_function(seq_var->val_func))
				break;
		}
		if (seq_var != m_track_vars.end())
			continue;
		try {
			if (GenomeTrack::is_2d(itrack_n_imdf->type)) {
				if (itrack_n_imdf->imdf2d) {
					if (!itrack_n_imdf->imdf2d->out_of_range)
						((GenomeTrack2D *)itrack_n_imdf->track)->read_interval(itrack_n_imdf->imdf2d->interval, m_band);
				} else
					((GenomeTrack2D *)itrack_n_imdf->track)->read_interval(m_interval2d, m_band);
			} else {
				if (itrack_n_imdf->imdf1d) {
					if (!itrack_n_imdf->imdf1d->out_of_range)
						((GenomeTrack1D *)itrack_n_imdf->track)->read_interval(itrack_n_imdf->imdf1d->interval);
				} else
					((GenomeTrack1D *)itrack_n_imdf->track)->read_interval(m_interval1d);
			}
		} catch (TGLException &e) {
			verror("%s", e.msg());
		}
		}

    // Collect sequence-based vtracks for potential batch processing
    vector<Track_var*> kmer_vtracks;
    vector<Track_var*> pwm_vtracks;

    for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
        if (TrackExpressionVars::is_pwm_function(ivar->val_func)) {
            pwm_vtracks.push_back(&*ivar);
        } else if (TrackExpressionVars::is_kmer_function(ivar->val_func)) {
            kmer_vtracks.push_back(&*ivar);
        }
    }

    // Check if any sequence vtrack has a filter
    bool any_filtered = false;
    for (Track_var* ivar : pwm_vtracks) {
        if (ivar->filter) {
            any_filtered = true;
            break;
        }
    }
    if (!any_filtered) {
        for (Track_var* ivar : kmer_vtracks) {
            if (ivar->filter) {
                any_filtered = true;
                break;
            }
        }
    }

    // Batch process if we have multiple sequence vtracks (threshold: 4+) AND no filters
    // (batch processing with filters is complex due to variable-length segments)
    if (!any_filtered && kmer_vtracks.size() + pwm_vtracks.size() >= 4) {
        batch_process_sequence_vtracks(kmer_vtracks, pwm_vtracks, m_interval1d, idx);
    } else {
        // Process individually if too few
        for (Track_var* ivar : pwm_vtracks) {
            const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : m_interval1d;

            // Check if filter applies
            if (ivar->filter) {
                std::vector<GInterval> unmasked_parts;
                ivar->filter->subtract(seq_interval, unmasked_parts);

                if (unmasked_parts.empty()) {
                    // Completely masked
                    ivar->var[idx] = std::numeric_limits<double>::quiet_NaN();
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

                        if (ivar->val_func == Track_var::PWM_MAX) {
                            if (first || part_score > result) {
                                result = part_score;
                            }
                        } else if (ivar->val_func == Track_var::PWM_MAX_POS) {
                            // For position, take the position with max score
                            if (first || part_score > result) {
                                result = part_score;
                            }
                        } else if (ivar->val_func == Track_var::PWM_COUNT) {
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
        for (Track_var* ivar : kmer_vtracks) {
            const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : m_interval1d;

            if (ivar->kmer_counter) {
                // Check if filter applies
                if (ivar->filter) {
                    std::vector<GInterval> unmasked_parts;
                    ivar->filter->subtract(seq_interval, unmasked_parts);

                    if (unmasked_parts.empty()) {
                        // Completely masked
                        ivar->var[idx] = std::numeric_limits<double>::quiet_NaN();
                    } else if (unmasked_parts.size() == 1) {
                        // Single unmasked part
                        ivar->var[idx] = ivar->kmer_counter->score_interval(unmasked_parts[0], m_iu.get_chromkey());
                    } else {
                        // Multiple unmasked parts - score each and aggregate
                        double total_count = 0.0;
                        int64_t total_bases = 0;

                        for (const auto& part : unmasked_parts) {
                            double part_score = ivar->kmer_counter->score_interval(part, m_iu.get_chromkey());

                            if (ivar->val_func == Track_var::KMER_COUNT) {
                                // Sum counts
                                total_count += part_score;
                            } else {
                                // KMER_FRAC - need to weight by length
                                total_count += part_score * part.range();
                                total_bases += part.range();
                            }
                        }

                        if (ivar->val_func == Track_var::KMER_FRAC && total_bases > 0) {
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
                ivar->var[idx] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    // Process regular track vtracks
    for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
        // Skip sequence-based vtracks (already processed above)
        if (TrackExpressionVars::is_sequence_based_function(ivar->val_func)) {
            continue;
        }

		if (GenomeTrack::is_1d(ivar->track_n_imdf->type)) {
			GenomeTrack1D &track = *(GenomeTrack1D *)ivar->track_n_imdf->track;
			const GInterval &base_interval = ivar->track_n_imdf->imdf1d ?
				ivar->track_n_imdf->imdf1d->interval : m_interval1d;
			const int64_t base_start = base_interval.start;

			if (ivar->track_n_imdf->imdf1d && ivar->track_n_imdf->imdf1d->out_of_range)
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			else {
				// Check if filter applies and get unmasked parts
				std::vector<GInterval> unmasked_parts;
				bool has_filter = false;

				if (ivar->filter) {
					const GInterval &eval_interval = ivar->track_n_imdf->imdf1d ?
						ivar->track_n_imdf->imdf1d->interval : m_interval1d;

					ivar->filter->subtract(eval_interval, unmasked_parts);
					has_filter = true;

					if (unmasked_parts.empty()) {
						// Completely masked - return NaN
						ivar->var[idx] = numeric_limits<double>::quiet_NaN();
						continue;
					}
				}

				// If filter exists and resulted in multiple unmasked parts, aggregate across them
				// If only one part (or filter but no filtering needed), use normal path below
				if (has_filter && unmasked_parts.size() > 1) {
					// Aggregate over unmasked parts using helper methods
					double result = std::numeric_limits<double>::quiet_NaN();

					switch (ivar->val_func) {
					case Track_var::REG:
					case Track_var::PV:
						result = aggregate_avg_with_filter(track, unmasked_parts);
						break;
					case Track_var::STDDEV:
						result = aggregate_stddev_with_filter(track, unmasked_parts);
						break;
					case Track_var::SUM:
						result = aggregate_sum_with_filter(track, unmasked_parts);
						break;
					case Track_var::REG_MIN:
					case Track_var::PV_MIN:
						result = aggregate_min_with_filter(track, unmasked_parts);
						break;
					case Track_var::REG_MAX:
					case Track_var::PV_MAX:
						result = aggregate_max_with_filter(track, unmasked_parts);
						break;
					case Track_var::MAX_POS_ABS:
						result = aggregate_max_pos_abs_with_filter(track, unmasked_parts);
						break;
					case Track_var::MAX_POS_REL:
						result = aggregate_max_pos_rel_with_filter(track, unmasked_parts, base_start);
						break;
					case Track_var::MIN_POS_ABS:
						result = aggregate_min_pos_abs_with_filter(track, unmasked_parts);
						break;
					case Track_var::MIN_POS_REL:
						result = aggregate_min_pos_rel_with_filter(track, unmasked_parts, base_start);
						break;
					case Track_var::REG_NEAREST:
						// For nearest, use the first unmasked part (closest to start)
						track.read_interval(unmasked_parts[0]);
						result = track.last_nearest();
						break;
					case Track_var::QUANTILE:
						result = aggregate_quantile_with_filter(track, unmasked_parts, ivar->percentile);
						break;
					case Track_var::EXISTS:
						result = aggregate_exists_with_filter(track, unmasked_parts);
						break;
					case Track_var::SIZE:
						result = aggregate_size_with_filter(track, unmasked_parts);
						break;
					case Track_var::SAMPLE:
						result = aggregate_sample_with_filter(track, unmasked_parts);
						break;
					case Track_var::SAMPLE_POS_ABS:
						result = aggregate_sample_pos_abs_with_filter(track, unmasked_parts);
						break;
					case Track_var::SAMPLE_POS_REL:
						result = aggregate_sample_pos_rel_with_filter(track, unmasked_parts, base_start);
						break;
					case Track_var::FIRST:
						// For first, use the first unmasked part
						track.read_interval(unmasked_parts[0]);
						result = track.last_first();
						break;
					case Track_var::FIRST_POS_ABS:
						// For first position, use the first unmasked part
						track.read_interval(unmasked_parts[0]);
						result = track.last_first_pos();
						break;
					case Track_var::FIRST_POS_REL:
						// For first position relative, use the first unmasked part
						track.read_interval(unmasked_parts[0]);
						result = track.last_first_pos() - base_start;
						break;
					case Track_var::LAST:
						// For last, use the last unmasked part
						track.read_interval(unmasked_parts.back());
						result = track.last_last();
						break;
					case Track_var::LAST_POS_ABS:
						// For last position, use the last unmasked part
						track.read_interval(unmasked_parts.back());
						result = track.last_last_pos();
						break;
					case Track_var::LAST_POS_REL:
						// For last position relative, use the last unmasked part
						track.read_interval(unmasked_parts.back());
						result = track.last_last_pos() - base_start;
						break;
					// Sequence-based functions are already handled above
					default:
						if (!TrackExpressionVars::is_sequence_based_function(ivar->val_func))
							verror("Internal error: unsupported function %d", ivar->val_func);
						break;
					}

					ivar->var[idx] = result;

					// Restore track state after filter processing
					// The aggregate functions call track.read_interval() which modifies the shared track state.
					// Other vtracks sharing this Track_n_imdf need the original interval restored.
					track.read_interval(base_interval);
				} else if (has_filter && unmasked_parts.size() == 1) {
					// Single unmasked part - read just that part
					track.read_interval(unmasked_parts[0]);

					switch (ivar->val_func) {
					case Track_var::REG:
					case Track_var::PV:
						ivar->var[idx] = track.last_avg();
						break;
					case Track_var::REG_MIN:
					case Track_var::PV_MIN:
						ivar->var[idx] = track.last_min();
						break;
					case Track_var::REG_MAX:
					case Track_var::PV_MAX:
						ivar->var[idx] = track.last_max();
						break;
					case Track_var::MAX_POS_ABS:
						ivar->var[idx] = track.last_max_pos();
						break;
					case Track_var::MAX_POS_REL:
						ivar->var[idx] = track.last_max_pos() - base_start;
						break;
					case Track_var::MIN_POS_ABS:
						ivar->var[idx] = track.last_min_pos();
						break;
					case Track_var::MIN_POS_REL:
						ivar->var[idx] = track.last_min_pos() - base_start;
						break;
					case Track_var::REG_NEAREST:
						ivar->var[idx] = track.last_nearest();
						break;
					case Track_var::STDDEV:
						ivar->var[idx] = track.last_stddev();
						break;
					case Track_var::SUM:
						ivar->var[idx] = track.last_sum();
						break;
					case Track_var::QUANTILE:
						ivar->var[idx] = track.last_quantile(ivar->percentile);
						break;
					case Track_var::EXISTS:
						ivar->var[idx] = track.last_exists();
						break;
					case Track_var::SIZE:
						ivar->var[idx] = track.last_size();
						break;
					case Track_var::SAMPLE:
						ivar->var[idx] = track.last_sample();
						break;
					case Track_var::SAMPLE_POS_ABS:
						ivar->var[idx] = track.last_sample_pos();
						break;
					case Track_var::SAMPLE_POS_REL:
						ivar->var[idx] = track.last_sample_pos() - base_start;
						break;
					case Track_var::FIRST:
						ivar->var[idx] = track.last_first();
						break;
					case Track_var::FIRST_POS_ABS:
						ivar->var[idx] = track.last_first_pos();
						break;
					case Track_var::FIRST_POS_REL:
						ivar->var[idx] = track.last_first_pos() - base_start;
						break;
					case Track_var::LAST:
						ivar->var[idx] = track.last_last();
						break;
					case Track_var::LAST_POS_ABS:
						ivar->var[idx] = track.last_last_pos();
						break;
					case Track_var::LAST_POS_REL:
						ivar->var[idx] = track.last_last_pos() - base_start;
						break;
					// Sequence-based functions are already handled above
					default:
						if (!TrackExpressionVars::is_sequence_based_function(ivar->val_func))
							verror("Internal error: unsupported function %d", ivar->val_func);
						break;
					}

					// Restore track state after processing single filtered part
					track.read_interval(base_interval);
				} else {
					// No filter or single unmasked part - use normal path
					switch (ivar->val_func) {
				case Track_var::REG:
				case Track_var::PV:
					ivar->var[idx] = track.last_avg();
					break;
				case Track_var::REG_MIN:
				case Track_var::PV_MIN:
					ivar->var[idx] = track.last_min();
					break;
				case Track_var::REG_MAX:
				case Track_var::PV_MAX:
					ivar->var[idx] = track.last_max();
					break;
				case Track_var::MAX_POS_ABS:
					ivar->var[idx] = track.last_max_pos();
					break;
				case Track_var::MAX_POS_REL:
					ivar->var[idx] = track.last_max_pos() - base_start;
					break;
				case Track_var::MIN_POS_ABS:
					ivar->var[idx] = track.last_min_pos();
					break;
				case Track_var::MIN_POS_REL:
					ivar->var[idx] = track.last_min_pos() - base_start;
					break;
				case Track_var::REG_NEAREST:
					ivar->var[idx] = track.last_nearest();
					break;
				case Track_var::STDDEV:
					ivar->var[idx] = track.last_stddev();
					break;
				case Track_var::SUM:
					ivar->var[idx] = track.last_sum();
					break;
				case Track_var::QUANTILE:
					ivar->var[idx] = track.last_quantile(ivar->percentile);
					break;
				case Track_var::EXISTS:
					ivar->var[idx] = track.last_exists();
					break;
				case Track_var::SIZE:
					ivar->var[idx] = track.last_size();
					break;
				case Track_var::SAMPLE:
					ivar->var[idx] = track.last_sample();
					break;
				case Track_var::SAMPLE_POS_ABS:
					ivar->var[idx] = track.last_sample_pos();
					break;
				case Track_var::SAMPLE_POS_REL:
					ivar->var[idx] = track.last_sample_pos() - base_start;
					break;
				case Track_var::FIRST:
					ivar->var[idx] = track.last_first();
					break;
				case Track_var::FIRST_POS_ABS:
					ivar->var[idx] = track.last_first_pos();
					break;
				case Track_var::FIRST_POS_REL:
					ivar->var[idx] = track.last_first_pos() - base_start;
					break;
				case Track_var::LAST:
					ivar->var[idx] = track.last_last();
					break;
				case Track_var::LAST_POS_ABS:
					ivar->var[idx] = track.last_last_pos();
					break;
				case Track_var::LAST_POS_REL:
					ivar->var[idx] = track.last_last_pos() - base_start;
					break;
				// Sequence-based functions are already handled above
				default:
					if (!TrackExpressionVars::is_sequence_based_function(ivar->val_func))
						verror("Internal error: unsupported function %d", ivar->val_func);
					break;
				}

					if (ivar->requires_pv) {
						double val = ivar->var[idx];
						if (!std::isnan(val)) {
							int bin = ivar->pv_binned.binfinder.val2bin(val);
							if (bin < 0) {
								if (val <= ivar->pv_binned.binfinder.get_breaks().front())
									ivar->var[idx] = ivar->pv_binned.bins[0];
								else
									ivar->var[idx] = 1.;
							} else
								ivar->var[idx] = ivar->pv_binned.bins[bin];
						}
					}
				}
			}
		} else {
			GenomeTrack2D &track = *(GenomeTrack2D *)ivar->track_n_imdf->track;

			if (ivar->track_n_imdf->imdf2d && ivar->track_n_imdf->imdf2d->out_of_range)
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			else {
				switch (ivar->val_func) {
				case Track_var::REG:
					ivar->var[idx] = track.last_avg();
					break;
				case Track_var::REG_MIN:
					ivar->var[idx] = track.last_min();
					break;
				case Track_var::REG_MAX:
					ivar->var[idx] = track.last_max();
					break;
				case Track_var::WEIGHTED_SUM:
					ivar->var[idx] = track.last_weighted_sum();
					break;
				case Track_var::OCCUPIED_AREA:
					ivar->var[idx] = track.last_occupied_area();
					break;
				default:
					verror("Internal error: unsupported function %d", ivar->val_func);
				}
			}
		}
	}

	// set intervals variables
	for (Interv_vars::iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ++ivar) {
		if (ivar->val_func == Interv_var::DIST) {
			// if iterator modifier exists, iterator intervals might not come sorted => perform a binary search
			if (ivar->imdf1d) {
				const GInterval &interval = ivar->imdf1d->interval;
				double min_dist = numeric_limits<double>::max();
				double dist;
				int64_t coord = (interval.start + interval.end) / 2;
				GIntervals::const_iterator iinterv = lower_bound(ivar->sintervs.begin(), ivar->sintervs.end(), interval, GIntervals::compare_by_start_coord);

				if (iinterv != ivar->sintervs.end() && iinterv->chromid == interval.chromid)
					min_dist = iinterv->dist2coord(coord, ivar->dist_margin);

				if (iinterv != ivar->sintervs.begin() && (iinterv - 1)->chromid == interval.chromid) {
					dist = (iinterv - 1)->dist2coord(coord, ivar->dist_margin);
					if (fabs(min_dist) > fabs(dist))
						min_dist = dist;
				}

				// if min_dist == double_max then we haven't found an interval with the same chromosome as the iterator interval =>
				// we can skip the second binary search
				if (min_dist == numeric_limits<double>::max())
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				else {
					iinterv = lower_bound(ivar->eintervs.begin(), ivar->eintervs.end(), interval, GIntervals::compare_by_end_coord);

					if (iinterv != ivar->eintervs.end() && iinterv->chromid == interval.chromid) {
						dist = iinterv->dist2coord(coord, ivar->dist_margin);
						if (fabs(min_dist) > fabs(dist))
							min_dist = dist;
					}

					if (iinterv != ivar->eintervs.begin() && (iinterv - 1)->chromid == interval.chromid) {
						dist = (iinterv - 1)->dist2coord(coord, ivar->dist_margin);
						if (fabs(min_dist) > fabs(dist))
							min_dist = dist;
					}

					ivar->var[idx] = min_dist;
				}
			} else {
				const GIntervals *pintervs[2] = { &ivar->sintervs, &ivar->eintervs };
				GIntervals::const_iterator *piinterv[2] = { &ivar->siinterv, &ivar->eiinterv };
				double dist[2] = { 0, 0 };
				const GInterval &interval = m_interval1d;

				for (int i = 0; i < 2; ++i) {
					const GIntervals &intervs = *pintervs[i];
					GIntervals::const_iterator &iinterv = *piinterv[i];

					while (iinterv != intervs.end() && iinterv->chromid < interval.chromid)
						++iinterv;

					if (iinterv == intervs.end() || iinterv->chromid != interval.chromid)
						dist[i] = numeric_limits<double>::quiet_NaN();
					else {
						int64_t coord = (interval.start + interval.end) / 2;
						dist[i] = (double)iinterv->dist2coord(coord, ivar->dist_margin);
						GIntervals::const_iterator iinterv_next = iinterv + 1;

						while (iinterv_next != intervs.end() && iinterv_next->chromid == interval.chromid) {
							double dist_next = iinterv_next->dist2coord(coord, ivar->dist_margin);

							if (fabs(dist[i]) < fabs(dist_next))
								break;

							iinterv = iinterv_next;
							dist[i] = dist_next;
							++iinterv_next;
						}
					}
				}

				ivar->var[idx] = fabs(dist[0]) < fabs(dist[1]) ? dist[0] : dist[1];
			}
		} else if (ivar->val_func == Interv_var::DIST_CENTER) {
			// if iterator modifier exists, iterator intervals might not come sorted => perform a binary search
			if (ivar->imdf1d) {
				int64_t coord = (ivar->imdf1d->interval.start + ivar->imdf1d->interval.end) / 2;
				GInterval interval(ivar->imdf1d->interval.chromid, coord, coord + 1, 0);
				GIntervals::const_iterator iinterv = lower_bound(ivar->sintervs.begin(), ivar->sintervs.end(), interval, GIntervals::compare_by_start_coord);
				double dist = numeric_limits<double>::quiet_NaN();

				ivar->var[idx] = numeric_limits<double>::quiet_NaN();

				if (iinterv != ivar->sintervs.end() && iinterv->chromid == interval.chromid)
					dist = iinterv->dist2center(coord);

				if (dist != numeric_limits<double>::quiet_NaN() && iinterv != ivar->sintervs.begin() && (iinterv - 1)->chromid == interval.chromid)
					dist = (iinterv - 1)->dist2center(coord);

				ivar->var[idx] = dist;
			} else {
				int64_t coord = (m_interval1d.start + m_interval1d.end) / 2;
				GIntervals::const_iterator &iinterv = ivar->siinterv;
				double dist = numeric_limits<double>::quiet_NaN();

				while (iinterv != ivar->sintervs.end() && ivar->siinterv->chromid < m_interval1d.chromid)
					++iinterv;

				while (iinterv != ivar->sintervs.end() && iinterv->chromid == m_interval1d.chromid && iinterv->start <= coord) {
					if (iinterv->end > coord)
						dist = iinterv->dist2center(coord);
					++iinterv;
				}

				ivar->var[idx] = dist;
			}
		} else if (ivar->val_func == Interv_var::NEIGHBOR_COUNT) {
			const GInterval &interval = ivar->imdf1d ? ivar->imdf1d->interval : m_interval1d;

			if (ivar->imdf1d && ivar->imdf1d->out_of_range) {
				ivar->var[idx] = 0;
				continue;
			}

			std::vector<GInterval> eval_intervals;
			if (ivar->filter) {
				ivar->filter->subtract(interval, eval_intervals);
				if (eval_intervals.empty()) {
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
					continue;
				}
			} else {
				eval_intervals.push_back(interval);
			}

			size_t neighbor_count = 0;
			if (!ivar->imdf1d && !ivar->filter) {
				GIntervals::const_iterator &eiter = ivar->eiinterv;
				const GIntervals &expanded = ivar->eintervs;

				while (eiter != expanded.end() && (eiter->chromid < interval.chromid ||
					   (eiter->chromid == interval.chromid && eiter->end <= interval.start)))
					++eiter;

				GIntervals::const_iterator scan = eiter;
				while (scan != expanded.end() && scan->chromid == interval.chromid && scan->start < interval.end) {
					if (scan->end > interval.start)
						++neighbor_count;
					++scan;
				}
			} else {
				const GIntervals &expanded = ivar->eintervs;
				std::unordered_set<size_t> counted;
				counted.reserve(eval_intervals.size() * 2);

				for (const auto &eval_interval : eval_intervals) {
					auto it = lower_bound(expanded.begin(), expanded.end(), eval_interval, GIntervals::compare_by_start_coord);

					// Walk backward to the first interval on this chromosome
					// Since intervals are only sorted by start (not end), we can't make assumptions
					// about whether earlier intervals might have large spans that overlap the query
					while (it != expanded.begin()) {
						auto prev = it - 1;
						if (prev->chromid != eval_interval.chromid)
							break;
						--it;
					}

					// Scan forward checking all intervals on this chromosome for overlaps
					for (; it != expanded.end() && it->chromid == eval_interval.chromid; ++it) {
						if (it->end <= eval_interval.start)
							continue;
						if (it->start >= eval_interval.end)
							break;

						size_t expanded_idx = it - expanded.begin();
						if (counted.insert(expanded_idx).second)
							++neighbor_count;
					}
				}
			}

			ivar->var[idx] = static_cast<double>(neighbor_count);
		} else if (ivar->val_func == Interv_var::COVERAGE)
		{
			const GInterval &interval = ivar->imdf1d ? ivar->imdf1d->interval : m_interval1d;

			if (ivar->imdf1d && ivar->imdf1d->out_of_range)
			{
				ivar->var[idx] = 0;
				continue;
			}

			// Apply filter if present
			std::vector<GInterval> eval_intervals;
			if (ivar->filter) {
				ivar->filter->subtract(interval, eval_intervals);
				if (eval_intervals.empty()) {
					ivar->var[idx] = numeric_limits<double>::quiet_NaN();
					continue;
				}
			} else {
				eval_intervals.push_back(interval);
			}

			int64_t total_overlap = 0;
			int64_t total_unmasked_length = 0;

			// Calculate coverage for each unmasked sub-interval
			for (const auto &eval_interval : eval_intervals) {
				total_unmasked_length += eval_interval.range();
				int64_t sub_overlap = 0;
				GIntervals::const_iterator iinterv;

				// For non-sequential access or first access
				if (ivar->imdf1d || ivar->siinterv == ivar->sintervs.end())
				{
					iinterv = lower_bound(ivar->sintervs.begin(), ivar->sintervs.end(), eval_interval,
										  GIntervals::compare_by_start_coord);

					// Check previous interval too
					if (iinterv != ivar->sintervs.begin())
					{
						auto prev = iinterv - 1;
						if (prev->chromid == eval_interval.chromid && prev->end > eval_interval.start)
						{
							int64_t overlap_start = max(eval_interval.start, prev->start);
							int64_t overlap_end = min(eval_interval.end, prev->end);
							sub_overlap += overlap_end - overlap_start;
						}
					}
				}
				else
				{
					// For sequential access, start from last position
					iinterv = ivar->siinterv;

					// Skip past intervals from previous chromosomes
					while (iinterv != ivar->sintervs.end() && iinterv->chromid < eval_interval.chromid)
						++iinterv;

					// But check if we need to back up
					while (iinterv != ivar->sintervs.begin() &&
						   (iinterv - 1)->chromid == eval_interval.chromid &&
						   (iinterv - 1)->end > eval_interval.start)
					{
						--iinterv;
					}
				}

				// Check forward intervals
				while (iinterv != ivar->sintervs.end() &&
					   iinterv->chromid == eval_interval.chromid &&
					   iinterv->start < eval_interval.end)
				{
					if (iinterv->end > eval_interval.start)
					{
						int64_t overlap_start = max(eval_interval.start, iinterv->start);
						int64_t overlap_end = min(eval_interval.end, iinterv->end);
						sub_overlap += overlap_end - overlap_start;
					}
					++iinterv;
				}

				if (!ivar->imdf1d)
				{
					ivar->siinterv = iinterv;
				}

				total_overlap += sub_overlap;
			}

			ivar->var[idx] = (double)total_overlap / total_unmasked_length;
		}
	}

	// set value variables - now using GenomeTrackInMemory for DRY principle
	for (Value_vars::iterator ivar = m_value_vars.begin(); ivar != m_value_vars.end(); ++ivar) {
		const GInterval &interval = ivar->imdf1d ? ivar->imdf1d->interval : m_interval1d;

		if (ivar->imdf1d && ivar->imdf1d->out_of_range) {
			ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			continue;
		}

		// Apply filter if present
		if (ivar->filter) {
			std::vector<GInterval> eval_intervals;
			ivar->filter->subtract(interval, eval_intervals);

			if (eval_intervals.empty()) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
				continue;
			}

			// For filtered intervals, we need to combine results across multiple sub-intervals
			// This matches the pattern used for track variables with filters
			struct FilteredValueAggregate {
				double total_sum = 0;
				double total_sum_sq = 0;
				int64_t total_size = 0;
				double min_val = numeric_limits<double>::quiet_NaN();
				double max_val = numeric_limits<double>::quiet_NaN();
				double first_val = numeric_limits<double>::quiet_NaN();
				double last_val = numeric_limits<double>::quiet_NaN();
				double sample_val = numeric_limits<double>::quiet_NaN();
				double min_pos = numeric_limits<double>::quiet_NaN();
				double max_pos = numeric_limits<double>::quiet_NaN();
				double first_pos = numeric_limits<double>::quiet_NaN();
				double last_pos = numeric_limits<double>::quiet_NaN();
				double sample_pos = numeric_limits<double>::quiet_NaN();
				double nearest_val = numeric_limits<double>::quiet_NaN();
				bool exists = false;
				bool has_value = false;

				void add_interval(GenomeTrackInMemory &track) {
					double part_size = track.last_size();
					if (part_size <= 0 || std::isnan(part_size))
						return;

					double part_avg = track.last_avg();
					double part_sum = track.last_sum();
					double part_stddev = track.last_stddev();
					if (std::isnan(part_avg) || std::isnan(part_sum))
						return;
					double part_min = track.last_min();
					double part_max = track.last_max();
					double part_first = track.last_first();
					double part_last = track.last_last();
					double part_sample = track.last_sample();
					double part_min_pos = track.last_min_pos();
					double part_max_pos = track.last_max_pos();
					double part_first_pos = track.last_first_pos();
					double part_last_pos = track.last_last_pos();
					double part_sample_pos = track.last_sample_pos();
					double part_nearest = track.last_nearest();
					double part_exists = track.last_exists();

					total_sum += part_sum;
					total_size += static_cast<int64_t>(part_size);

					if (part_size > 1 && !std::isnan(part_stddev)) {
						double var_term = part_stddev * part_stddev * (part_size - 1);
						total_sum_sq += var_term + part_avg * part_avg * part_size;
					} else
						total_sum_sq += part_avg * part_avg * part_size;

					if (!std::isnan(part_min)) {
						if (std::isnan(min_val) || part_min < min_val)
							min_val = part_min;
					}
					if (!std::isnan(part_max)) {
						if (std::isnan(max_val) || part_max > max_val)
							max_val = part_max;
					}

					if (!std::isnan(part_first) && std::isnan(first_val))
						first_val = part_first;
					if (!std::isnan(part_last))
						last_val = part_last;

					if (!std::isnan(part_sample) && std::isnan(sample_val))
						sample_val = part_sample;

					if (!std::isnan(part_min_pos)) {
						if (std::isnan(min_pos) || part_min_pos < min_pos)
							min_pos = part_min_pos;
					}
					if (!std::isnan(part_max_pos)) {
						if (std::isnan(max_pos) || part_max_pos > max_pos)
							max_pos = part_max_pos;
					}
					if (!std::isnan(part_first_pos) && std::isnan(first_pos))
						first_pos = part_first_pos;
					if (!std::isnan(part_last_pos))
						last_pos = part_last_pos;

					if (!std::isnan(part_sample_pos) && std::isnan(sample_pos))
						sample_pos = part_sample_pos;

					if (!std::isnan(part_nearest) && std::isnan(nearest_val))
						nearest_val = part_nearest;

					exists = exists || (!std::isnan(part_exists) && part_exists != 0);
					has_value = true;
				}
			} agg;

			for (const auto &eval_interval : eval_intervals) {
				ivar->track->read_interval(eval_interval);
				agg.add_interval(*ivar->track);
			}

			// Finalize result for filtered case
			if (!agg.has_value) {
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			} else {
				switch (ivar->val_func) {
					case Value_var::AVG:
						ivar->var[idx] = agg.total_size > 0 ? agg.total_sum / agg.total_size : numeric_limits<double>::quiet_NaN();
						break;
					case Value_var::STDDEV:
						if (agg.total_size > 1) {
							double mean = agg.total_sum / agg.total_size;
							double variance = agg.total_sum_sq / (agg.total_size - 1) - mean * mean * ((double)agg.total_size / (agg.total_size - 1));
							ivar->var[idx] = sqrt(std::max(0.0, variance));
						} else
							ivar->var[idx] = numeric_limits<double>::quiet_NaN();
						break;
					case Value_var::SUM:
						ivar->var[idx] = agg.total_sum;
						break;
					case Value_var::MIN:
						ivar->var[idx] = agg.min_val;
						break;
					case Value_var::MAX:
						ivar->var[idx] = agg.max_val;
						break;
					case Value_var::EXISTS:
						ivar->var[idx] = agg.exists ? 1.0 : 0.0;
						break;
					case Value_var::SIZE:
						ivar->var[idx] = agg.total_size;
						break;
					case Value_var::FIRST:
						ivar->var[idx] = agg.first_val;
						break;
					case Value_var::LAST:
						ivar->var[idx] = agg.last_val;
						break;
					case Value_var::SAMPLE:
						ivar->var[idx] = agg.sample_val;
						break;
					case Value_var::FIRST_POS_ABS:
						ivar->var[idx] = agg.first_pos;
						break;
					case Value_var::FIRST_POS_REL:
						ivar->var[idx] = std::isnan(agg.first_pos) ? numeric_limits<double>::quiet_NaN() : agg.first_pos - interval.start;
						break;
					case Value_var::LAST_POS_ABS:
						ivar->var[idx] = agg.last_pos;
						break;
					case Value_var::LAST_POS_REL:
						ivar->var[idx] = std::isnan(agg.last_pos) ? numeric_limits<double>::quiet_NaN() : agg.last_pos - interval.start;
						break;
					case Value_var::SAMPLE_POS_ABS:
						ivar->var[idx] = agg.sample_pos;
						break;
					case Value_var::SAMPLE_POS_REL:
						ivar->var[idx] = std::isnan(agg.sample_pos) ? numeric_limits<double>::quiet_NaN() : agg.sample_pos - interval.start;
						break;
					case Value_var::MIN_POS_ABS:
						ivar->var[idx] = agg.min_pos;
						break;
					case Value_var::MIN_POS_REL:
						ivar->var[idx] = std::isnan(agg.min_pos) ? numeric_limits<double>::quiet_NaN() : agg.min_pos - interval.start;
						break;
					case Value_var::MAX_POS_ABS:
						ivar->var[idx] = agg.max_pos;
						break;
					case Value_var::MAX_POS_REL:
						ivar->var[idx] = std::isnan(agg.max_pos) ? numeric_limits<double>::quiet_NaN() : agg.max_pos - interval.start;
						break;
					case Value_var::NEAREST:
						ivar->var[idx] = agg.nearest_val;
						break;
					case Value_var::QUANTILE: {
						// For quantile with filter, just use first non-empty result
						float val = numeric_limits<float>::quiet_NaN();
						for (const auto &eval_interval : eval_intervals) {
							ivar->track->read_interval(eval_interval);
							val = ivar->track->last_quantile(ivar->percentile);
							if (!std::isnan(val))
								break;
						}
						ivar->var[idx] = val;
						break;
					}
					default:
						ivar->var[idx] = numeric_limits<double>::quiet_NaN();
						break;
				}
			}
		} else {
			// No filter - simple case, use track interface directly
			ivar->track->read_interval(interval);

			// Extract result based on function
			switch (ivar->val_func) {
				case Value_var::AVG:
					ivar->var[idx] = ivar->track->last_avg();
					break;
				case Value_var::MIN:
					ivar->var[idx] = ivar->track->last_min();
					break;
				case Value_var::MAX:
					ivar->var[idx] = ivar->track->last_max();
					break;
				case Value_var::SUM:
					ivar->var[idx] = ivar->track->last_sum();
					break;
				case Value_var::STDDEV:
					ivar->var[idx] = ivar->track->last_stddev();
					break;
				case Value_var::QUANTILE:
					ivar->var[idx] = ivar->track->last_quantile(ivar->percentile);
					break;
				case Value_var::NEAREST:
					ivar->var[idx] = ivar->track->last_nearest();
					break;
				case Value_var::EXISTS:
					ivar->var[idx] = ivar->track->last_exists();
					break;
				case Value_var::SIZE:
					ivar->var[idx] = ivar->track->last_size();
					break;
				case Value_var::FIRST:
					ivar->var[idx] = ivar->track->last_first();
					break;
				case Value_var::LAST:
					ivar->var[idx] = ivar->track->last_last();
					break;
				case Value_var::SAMPLE:
					ivar->var[idx] = ivar->track->last_sample();
					break;
				case Value_var::FIRST_POS_ABS:
					ivar->var[idx] = ivar->track->last_first_pos();
					break;
				case Value_var::FIRST_POS_REL:
					ivar->var[idx] = ivar->track->last_first_pos() - interval.start;
					break;
				case Value_var::LAST_POS_ABS:
					ivar->var[idx] = ivar->track->last_last_pos();
					break;
				case Value_var::LAST_POS_REL:
					ivar->var[idx] = ivar->track->last_last_pos() - interval.start;
					break;
				case Value_var::SAMPLE_POS_ABS:
					ivar->var[idx] = ivar->track->last_sample_pos();
					break;
				case Value_var::SAMPLE_POS_REL:
					ivar->var[idx] = ivar->track->last_sample_pos() - interval.start;
					break;
				case Value_var::MIN_POS_ABS:
					ivar->var[idx] = ivar->track->last_min_pos();
					break;
				case Value_var::MIN_POS_REL:
					ivar->var[idx] = ivar->track->last_min_pos() - interval.start;
					break;
				case Value_var::MAX_POS_ABS:
					ivar->var[idx] = ivar->track->last_max_pos();
					break;
				case Value_var::MAX_POS_REL:
					ivar->var[idx] = ivar->track->last_max_pos() - interval.start;
					break;
				default:
					verror("Internal error: unsupported value variable function %d", ivar->val_func);
			}
		}
	}
}

void TrackExpressionVars::batch_process_sequence_vtracks(vector<Track_var*> &kmer_vtracks,
                                                         vector<Track_var*> &pwm_vtracks,
                                                         const GInterval &interval,
                                                         unsigned idx)
{
	// Group kmers by: (interval, extend_params, strand)
	// Key: "chromid:start-end:max_extension:strand"
	std::unordered_map<std::string, vector<Track_var*>> kmer_groups;

	for (Track_var* var : kmer_vtracks) {
		if (!var->kmer_counter) continue;

		const GInterval &base_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
		char strand = var->kmer_counter->get_strand();
		int64_t extension = var->kmer_counter->get_extend() ? (int64_t)(var->kmer_counter->get_kmer().length() - 1) : 0;

		std::string key = std::to_string(base_interval.chromid) + ":" +
		                  std::to_string(base_interval.start) + "-" +
		                  std::to_string(base_interval.end) + ":" +
		                  std::to_string(extension) + ":" +
		                  std::to_string((int)strand);
		kmer_groups[key].push_back(var);
	}

	// Process each kmer group with HashMap optimization
	for (const auto& [key, group] : kmer_groups) {
		if (group.empty()) continue;

		Track_var* first = group[0];
		const GInterval &base_interval = first->seq_imdf1d ? first->seq_imdf1d->interval : interval;
		char strand = first->kmer_counter->get_strand();
		int64_t extension = first->kmer_counter->get_extend() ? (int64_t)(first->kmer_counter->get_kmer().length() - 1) : 0;

		// Calculate fetch interval with extension
		GInterval fetch_interval = base_interval;
		if (extension > 0) {
			fetch_interval.start = std::max((int64_t)0, base_interval.start - extension);
			fetch_interval.end = std::min((int64_t)m_iu.get_chromkey().get_chrom_size(base_interval.chromid),
			                              base_interval.end + extension);
		}

		// Build HashMap for forward and/or reverse strands
		std::unordered_map<std::string, std::vector<size_t>> fwd_positions;
		std::unordered_map<std::string, std::vector<size_t>> rev_positions;

		// Get max kmer length for this group
		size_t max_kmer_len = 0;
		for (Track_var* var : group) {
			max_kmer_len = std::max(max_kmer_len, var->kmer_counter->get_kmer().length());
		}

		// Fetch and process forward strand if needed
		if (strand == 0 || strand == 1) {
			GInterval fwd_interval = fetch_interval;
			fwd_interval.strand = 1;
			std::vector<char> seq;
			m_shared_seqfetch.read_interval(fwd_interval, m_iu.get_chromkey(), seq);
			std::string sequence(seq.begin(), seq.end());
			std::transform(sequence.begin(), sequence.end(), sequence.begin(),
			               [](unsigned char c) { return std::toupper(c); });

			// Build HashMap: scan once, record all kmer positions
			for (size_t pos = 0; pos + max_kmer_len <= sequence.length(); pos++) {
				for (size_t k = 1; k <= max_kmer_len && pos + k <= sequence.length(); k++) {
					std::string kmer = sequence.substr(pos, k);
					fwd_positions[kmer].push_back(pos);
				}
			}
		}

		// Fetch and process reverse strand if needed
		if (strand == 0 || strand == -1) {
			GInterval rev_interval = fetch_interval;
			rev_interval.strand = -1;
			std::vector<char> seq;
			m_shared_seqfetch.read_interval(rev_interval, m_iu.get_chromkey(), seq);
			std::string sequence(seq.begin(), seq.end());
			std::transform(sequence.begin(), sequence.end(), sequence.begin(),
			               [](unsigned char c) { return std::toupper(c); });

			// Build HashMap for reverse strand
			for (size_t pos = 0; pos + max_kmer_len <= sequence.length(); pos++) {
				for (size_t k = 1; k <= max_kmer_len && pos + k <= sequence.length(); k++) {
					std::string kmer = sequence.substr(pos, k);
					rev_positions[kmer].push_back(pos);
				}
			}
		}

		// Distribute to all kmer counters in this group
		for (Track_var* var : group) {
			const GInterval &var_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
			const std::string &target_kmer = var->kmer_counter->get_kmer();

			var->var[idx] = var->kmer_counter->score_from_positions(
				fwd_positions[target_kmer],
				rev_positions[target_kmer],
				var_interval,
				fetch_interval);
		}
	}

	// Process PWM vtracks using score_interval
	// The shared cache prevents redundant disk I/O
	for (Track_var* var : pwm_vtracks) {
		const GInterval &seq_interval = var->seq_imdf1d ? var->seq_imdf1d->interval : interval;
		var->var[idx] = var->pwm_scorer->score_interval(seq_interval, m_iu.get_chromkey());
	}
}

// Filter aggregation helper methods
double TrackExpressionVars::aggregate_avg_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	long double total_weight = 0;
	long double total_weighted_sum = 0;

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_avg = track.last_avg();
		if (!std::isnan(part_avg)) {
			double part_len = part.end - part.start;
			total_weighted_sum += part_avg * part_len;
			total_weight += part_len;
		}
	}

	return (total_weight > 0) ? (total_weighted_sum / total_weight) : numeric_limits<double>::quiet_NaN();
}

double TrackExpressionVars::aggregate_sum_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double total_sum = 0;
	bool has_value = false;

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_sum = track.last_sum();
		if (!std::isnan(part_sum)) {
			total_sum += part_sum;
			has_value = true;
		}
	}

	return has_value ? total_sum : numeric_limits<double>::quiet_NaN();
}

double TrackExpressionVars::aggregate_min_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double min_val = numeric_limits<double>::infinity();

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_min = track.last_min();
		if (!std::isnan(part_min)) {
			min_val = std::min(min_val, part_min);
		}
	}

	return (min_val != numeric_limits<double>::infinity()) ? min_val : numeric_limits<double>::quiet_NaN();
}

double TrackExpressionVars::aggregate_max_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double max_val = -numeric_limits<double>::infinity();

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_max = track.last_max();
		if (!std::isnan(part_max)) {
			max_val = std::max(max_val, part_max);
		}
	}

	return (max_val != -numeric_limits<double>::infinity()) ? max_val : numeric_limits<double>::quiet_NaN();
}

bool TrackExpressionVars::find_best_max_pos_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double &best_pos)
{
	double best_val = -numeric_limits<double>::infinity();
	bool has_value = false;

	for (const auto &part : parts) {
		track.read_interval(part);
		double part_max = track.last_max();
		double part_pos = track.last_max_pos();

		if (std::isnan(part_max) || std::isnan(part_pos))
			continue;

		if (!has_value || part_max > best_val || (part_max == best_val && part_pos < best_pos)) {
			best_val = part_max;
			best_pos = part_pos;
			has_value = true;
		}
	}

	return has_value;
}

double TrackExpressionVars::aggregate_max_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_max_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos;
}

double TrackExpressionVars::aggregate_max_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_max_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos - base_start;
}

bool TrackExpressionVars::find_best_min_pos_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double &best_pos)
{
	double best_val = numeric_limits<double>::infinity();
	bool has_value = false;

	for (const auto &part : parts) {
		track.read_interval(part);
		double part_min = track.last_min();
		double part_pos = track.last_min_pos();

		if (std::isnan(part_min) || std::isnan(part_pos))
			continue;

		if (!has_value || part_min < best_val || (part_min == best_val && part_pos < best_pos)) {
			best_val = part_min;
			best_pos = part_pos;
			has_value = true;
		}
	}

	return has_value;
}

double TrackExpressionVars::aggregate_min_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_min_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos;
}

double TrackExpressionVars::aggregate_min_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double best_pos = numeric_limits<double>::quiet_NaN();
	if (!find_best_min_pos_with_filter(track, parts, best_pos))
		return numeric_limits<double>::quiet_NaN();

	return best_pos - base_start;
}

double TrackExpressionVars::aggregate_stddev_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	long double total_weight = 0;
	long double total_weighted_sum = 0;
	long double M2 = 0;  // For variance calculation

	for (const auto& part : parts) {
		track.read_interval(part);
		double part_avg = track.last_avg();
		double part_stddev = track.last_stddev();
		double part_len = part.end - part.start;

		if (!std::isnan(part_avg) && !std::isnan(part_stddev)) {
			// Welford online algorithm for combining variances
			double delta = part_avg - (total_weight > 0 ? total_weighted_sum / total_weight : 0);
			total_weight += part_len;
			total_weighted_sum += part_avg * part_len;
			M2 += part_stddev * part_stddev * part_len + part_len * total_weight / (total_weight + part_len) * delta * delta;
		}
	}

	return (total_weight > 0) ? std::sqrt(M2 / total_weight) : numeric_limits<double>::quiet_NaN();
}

double TrackExpressionVars::aggregate_quantile_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, double percentile)
{
	// Create a combined percentiler and add all samples from all parts
	StreamPercentiler<float> combined_sp;
	combined_sp.init(m_iu.get_max_data_size(),
	                 m_iu.get_quantile_edge_data_size(),
	                 m_iu.get_quantile_edge_data_size());

	for (const auto& part : parts) {
		track.read_interval(part);

		const auto& sp = track.get_percentiler();
		const auto& samples = sp.samples();
		const auto& lowest = sp.lowest_vals();
		const auto& highest = sp.highest_vals();

		// Add all values to combined percentiler
		for (const auto& val : lowest) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
		for (const auto& val : samples) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
		for (const auto& val : highest) {
			if (!std::isnan(val)) {
				combined_sp.add(val, unif_rand);
			}
		}
	}

	if (combined_sp.stream_size() > 0) {
		bool is_estimated;
		return combined_sp.get_percentile(percentile, is_estimated);
	}

	return numeric_limits<double>::quiet_NaN();
}

double TrackExpressionVars::aggregate_exists_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Check if any value exists in any part
	for (const auto& part : parts) {
		track.read_interval(part);
		float exists = track.last_exists();
		if (exists == 1.0f) {
			return 1.0;
		}
	}
	return 0.0;
}

double TrackExpressionVars::aggregate_size_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Sum the sizes across all parts
	double total_size = 0.0;
	for (const auto& part : parts) {
		track.read_interval(part);
		total_size += track.last_size();
	}
	return total_size;
}

double TrackExpressionVars::aggregate_sample_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Collect all non-NaN values from all parts and sample one
	vector<float> all_values;
	for (const auto& part : parts) {
		track.read_interval(part);
		float val = track.last_sample();
		if (!std::isnan(val)) {
			all_values.push_back(val);
		}
	}

	if (all_values.empty()) {
		return numeric_limits<double>::quiet_NaN();
	}

	// Use R's random number generator for reproducibility
	GetRNGstate();
	int idx = (int)(unif_rand() * all_values.size());
	PutRNGstate();

	return all_values[idx];
}

double TrackExpressionVars::aggregate_sample_pos_abs_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts)
{
	// Collect all non-NaN positions from all parts and sample one
	vector<double> all_positions;
	for (const auto& part : parts) {
		track.read_interval(part);
		double pos = track.last_sample_pos();
		if (!std::isnan(pos)) {
			all_positions.push_back(pos);
		}
	}

	if (all_positions.empty()) {
		return numeric_limits<double>::quiet_NaN();
	}

	// Use R's random number generator for reproducibility
	GetRNGstate();
	int idx = (int)(unif_rand() * all_positions.size());
	PutRNGstate();

	return all_positions[idx];
}

double TrackExpressionVars::aggregate_sample_pos_rel_with_filter(GenomeTrack1D &track, const vector<GInterval> &parts, int64_t base_start)
{
	double abs_pos = aggregate_sample_pos_abs_with_filter(track, parts);
	if (std::isnan(abs_pos)) {
		return abs_pos;
	}
	return abs_pos - base_start;
}
