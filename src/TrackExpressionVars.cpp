#include <cstdint>
#include <cmath>
#include <unistd.h>
#include <unordered_map>

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
#include "GenomeSeqFetch.h"

const char *TrackExpressionVars::Track_var::FUNC_NAMES[TrackExpressionVars::Track_var::NUM_FUNCS] = {
	"avg", "min", "max", "nearest", "stddev", "sum", "quantile",
	"global.percentile", "global.percentile.min", "global.percentile.max",
	"weighted.sum", "area", "pwm", "pwm.max", "pwm.max.pos", "kmer.count", "kmer.frac"};

const char *TrackExpressionVars::Interv_var::FUNC_NAMES[TrackExpressionVars::Interv_var::NUM_FUNCS] = { "distance", "distance.center", "coverage" };

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
		if (ivar->val_func == Interv_var::DIST)
			ivar->eiinterv = ivar->eintervs.begin();
	}
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

	// First check if this is a PWM function
	SEXP rfunc = get_rvector_col(rvtrack, "func", vtrack.c_str(), false);
    if (!Rf_isNull(rfunc) && Rf_isString(rfunc)) {
        string func = CHAR(STRING_ELT(rfunc, 0));
        transform(func.begin(), func.end(), func.begin(), ::tolower);
        
        if (func == "pwm" || func == "pwm.max" || func == "pwm.max.pos") {
            // Create the Track_var without a Track_n_imdf
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            var.val_func = (func == "pwm" ? Track_var::PWM : func == "pwm.max" ? Track_var::PWM_MAX : Track_var::PWM_MAX_POS);
            var.track_n_imdf = nullptr;  // No track needed for PWM
            var.seq_imdf1d = nullptr;
            
            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);
			if (!Rf_isNewList(rparams))	{
				verror("Virtual track %s: PWM functions require a list parameter with pssm matrix", vtrack.c_str());
			}

			// Get PSSM matrix from params
			SEXP rpssm = VECTOR_ELT(rparams, findListElementIndex(rparams, "pssm"));
			if (!Rf_isMatrix(rpssm)){
				rdb::verror("Virtual track %s: PWM functions require a matrix parameter", vtrack.c_str());
			}

			// Get bidirect parameter
			SEXP rbidirect = VECTOR_ELT(rparams, findListElementIndex(rparams, "bidirect"));
			bool bidirect = true; // default value
			if (rbidirect != R_NilValue){
				if (!Rf_isLogical(rbidirect))
					rdb::verror("Virtual track %s: bidirect parameter must be logical", vtrack.c_str());
				bidirect = LOGICAL(rbidirect)[0];
			}

			// Get extend parameter
			SEXP rextend = VECTOR_ELT(rparams, findListElementIndex(rparams, "extend"));
			bool extend = false;
			if (rextend != R_NilValue){
				if (!Rf_isLogical(rextend))
					rdb::verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
				extend = LOGICAL(rextend)[0];
			}

			// Get strand parameter (numeric)
			SEXP rstrand = VECTOR_ELT(rparams, findListElementIndex(rparams, "strand"));
			char strand = 0;
			if (rstrand != R_NilValue){
				if (!Rf_isReal(rstrand) || Rf_length(rstrand) != 1)
					rdb::verror("Virtual track %s: strand parameter must be numeric", vtrack.c_str());
				strand = (char)REAL(rstrand)[0];
			}

			// Optional spatial parameters
			std::vector<float> spat_factor_vec;
			int spat_bin = 1;

			// Get spat_factor (numeric vector)
			int spat_idx = findListElementIndex(rparams, "spat_factor");
			SEXP rspat = R_NilValue;
			if (spat_idx >= 0) {
				rspat = VECTOR_ELT(rparams, spat_idx);
			}
			if (rspat != R_NilValue) {
				if (!Rf_isReal(rspat))
					rdb::verror("Virtual track %s: spat_factor must be a numeric vector", vtrack.c_str());
				int n = Rf_length(rspat);
				if (n <= 0)
					rdb::verror("Virtual track %s: spat_factor must have at least one element", vtrack.c_str());
				spat_factor_vec.resize(n);
				for (int i = 0; i < n; ++i) {
					spat_factor_vec[i] = REAL(rspat)[i];
					if (spat_factor_vec[i] <= 0)
						rdb::verror("Virtual track %s: all spat_factor values must be positive", vtrack.c_str());
				}

				// Get spat_bin (integer scalar)
				int bin_idx = findListElementIndex(rparams, "spat_bin");
				SEXP rbin = R_NilValue;
				if (bin_idx >= 0) {
					rbin = VECTOR_ELT(rparams, bin_idx);
				}
				if (rbin != R_NilValue) {
					if (!Rf_isInteger(rbin) && !Rf_isReal(rbin))
						rdb::verror("Virtual track %s: spat_bin must be numeric", vtrack.c_str());
					spat_bin = (int)(Rf_isReal(rbin) ? REAL(rbin)[0] : INTEGER(rbin)[0]);
					if (spat_bin <= 0)
						rdb::verror("Virtual track %s: spat_bin must be > 0", vtrack.c_str());
				}
			}

			// Optional spat_min/spat_max (integers)
			int spat_min = 0;
			int spat_max = 1000000;
			bool has_range = false;

			int smin_idx = findListElementIndex(rparams, "spat_min");
			SEXP rsmin = R_NilValue;
			if (smin_idx >= 0) {
				rsmin = VECTOR_ELT(rparams, smin_idx);
			}
			if (rsmin != R_NilValue) {
				if (!Rf_isInteger(rsmin) && !Rf_isReal(rsmin))
					rdb::verror("Virtual track %s: spat_min must be numeric", vtrack.c_str());
				spat_min = (int)(Rf_isReal(rsmin) ? REAL(rsmin)[0] : INTEGER(rsmin)[0]);
				has_range = true;
			}

			int smax_idx = findListElementIndex(rparams, "spat_max");
			SEXP rsmax = R_NilValue;
			if (smax_idx >= 0) {
				rsmax = VECTOR_ELT(rparams, smax_idx);
			}
			if (rsmax != R_NilValue) {
				if (!Rf_isInteger(rsmax) && !Rf_isReal(rsmax))
					rdb::verror("Virtual track %s: spat_max must be numeric", vtrack.c_str());
				spat_max = (int)(Rf_isReal(rsmax) ? REAL(rsmax)[0] : INTEGER(rsmax)[0]);
				has_range = true;
			}

			// Create PSSM and initialize PWM scorer
			DnaPSSM pssm = PWMScorer::create_pssm_from_matrix(rpssm);
			pssm.set_bidirect(bidirect);

			// Apply optional scan range if provided
			if (has_range) {
				pssm.set_range(spat_min, spat_max);
			}

			// Construct scorer with shared sequence fetcher for caching
			var.pwm_scorer = std::make_unique<PWMScorer>(
				pssm,
				&m_shared_seqfetch,
				extend,
				func == "pwm" ? PWMScorer::TOTAL_LIKELIHOOD :
				func == "pwm.max" ? PWMScorer::MAX_LIKELIHOOD :
				PWMScorer::MAX_LIKELIHOOD_POS,
				strand,
				spat_factor_vec,
				spat_bin
			);

            // Parse optional iterator modifier (sshift/eshift) for sequence-based vtracks
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);
            
            var.percentile = numeric_limits<double>::quiet_NaN();
            var.requires_pv = false;
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

			if (Rf_isNull(rparams))
				rdb::verror("Virtual track %s: function %s requires a parameter (kmer string)", vtrack.c_str(), func.c_str());

			// Get extension parameter if exists, default to true (similar to PWM behavior)
			bool extend = true;
			char strand = 0;
			if (Rf_isNewList(rparams))
			{
				// Handle as list params
				SEXP rextend = VECTOR_ELT(rparams, findListElementIndex(rparams, "extend"));
				if (rextend != R_NilValue)
				{
					if (!Rf_isLogical(rextend))
						rdb::verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
					extend = LOGICAL(rextend)[0];
				}

				// Extract kmer string from the list parameters
				SEXP rkmer = VECTOR_ELT(rparams, findListElementIndex(rparams, "kmer"));
				if (rkmer == R_NilValue || !Rf_isString(rkmer) || Rf_length(rkmer) != 1)
					rdb::verror("Virtual track %s: invalid parameter used for function %s (must be a kmer string)",
								vtrack.c_str(), func.c_str());

				const char *kmer = CHAR(STRING_ELT(rkmer, 0));

				SEXP rstrand = VECTOR_ELT(rparams, findListElementIndex(rparams, "strand"));
				if (rstrand != R_NilValue)
				{
					if (!Rf_isNumeric(rstrand) || Rf_length(rstrand) != 1)
						rdb::verror("Virtual track %s: strand parameter must be -1, 0, or 1", vtrack.c_str());
					strand = (char)REAL(rstrand)[0];
					if (strand != -1 && strand != 0 && strand != 1)
						rdb::verror("Virtual track %s: strand parameter must be -1, 0, or 1", vtrack.c_str());
				}

				KmerCounter::CountMode mode = func == "kmer.count" ? KmerCounter::SUM : KmerCounter::FRACTION;
				var.kmer_counter = std::make_unique<KmerCounter>(kmer, &m_shared_seqfetch, mode, extend, strand);

			}
			else if (Rf_isString(rparams) && Rf_length(rparams) == 1)
			{
				// Handle direct string parameter (backward compatibility)
				const char *kmer = CHAR(STRING_ELT(rparams, 0));
				KmerCounter::CountMode mode = func == "kmer.count" ? KmerCounter::SUM : KmerCounter::FRACTION;
				var.kmer_counter = std::make_unique<KmerCounter>(kmer, &m_shared_seqfetch, mode, extend, strand);
			}
			else
			{
				rdb::verror("Virtual track %s: invalid parameter used for function %s (must be a kmer string)",
							vtrack.c_str(), func.c_str());
			}

            // Parse optional iterator modifier (sshift/eshift) for sequence-based vtracks
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);

            var.requires_pv = false;
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
				
			} else	if (ifunc == Track_var::PWM || ifunc == Track_var::PWM_MAX || ifunc == Track_var::PWM_MAX_POS) {
				var.percentile = numeric_limits<double>::quiet_NaN();
				if (!Rf_isMatrix(rparams))
					verror("Virtual track %s: PWM functions require a matrix parameter", vtrack.c_str());

				// Read extend parameter from vtrack
				SEXP rextend = get_rvector_col(rvtrack, "extend", vtrack.c_str(), false);
				bool extend = true; // default value (matching constructor default)
				if (!Rf_isNull(rextend)) {
					if (!Rf_isLogical(rextend))
						verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
					extend = LOGICAL(rextend)[0];
				}

				DnaPSSM pssm = PWMScorer::create_pssm_from_matrix(rparams);
				var.pwm_scorer = std::make_unique<PWMScorer>(
					pssm,
					&m_shared_seqfetch,
					extend,
					ifunc == Track_var::PWM ? PWMScorer::TOTAL_LIKELIHOOD : PWMScorer::MAX_LIKELIHOOD);
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
	
	} else {
		verror("Virtual track %s: invalid function %s used with intervals", vtrack.c_str(), func.c_str());
	}

	return var;
}

void TrackExpressionVars::register_track_functions()
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		// Skip PWM variables since they don't have associated tracks
        if (ivar->val_func == Track_var::PWM || ivar->val_func == Track_var::PWM_MAX || ivar->val_func == Track_var::PWM_MAX_POS || ivar->val_func == Track_var::KMER_COUNT || ivar->val_func == Track_var::KMER_FRAC) {
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
		case Track_var::PWM:
		case Track_var::PWM_MAX:
		case Track_var::PWM_MAX_POS:
		case Track_var::KMER_COUNT:
		case Track_var::KMER_FRAC:
			// PWM functions work directly on sequences, no need to register track functions
			break;
		default:
			verror("Unrecognized virtual track function");
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
	    // Skip iterator validation for PWM variables since they don't have tracks or imdf
        if (itrack_var->val_func == Track_var::PWM || itrack_var->val_func == Track_var::PWM_MAX || itrack_var->val_func == Track_var::PWM_MAX_POS || itrack_var->val_func == Track_var::KMER_COUNT || itrack_var->val_func == Track_var::KMER_FRAC) {
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
		// Skip PWM tracks
		if ((ivar->val_func == Track_var::PWM || ivar->val_func == Track_var::PWM_MAX || ivar->val_func == Track_var::PWM_MAX_POS || ivar->val_func == Track_var::KMER_COUNT || ivar->val_func == Track_var::KMER_FRAC))
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
}

void TrackExpressionVars::start_chrom(const GInterval &interval)
{
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		// Skip track initialization for PWM-only tracks
		bool is_pwm_track = false;
		for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar)
		{
			if (ivar->track_n_imdf == &(*itrack_n_imdf) &&
				(ivar->val_func == Track_var::PWM || ivar->val_func == Track_var::PWM_MAX || ivar->val_func == Track_var::PWM_MAX_POS || ivar->val_func == Track_var::KMER_COUNT || ivar->val_func == Track_var::KMER_FRAC))
			{
				is_pwm_track = true;
				break;
			}
		}
		if (is_pwm_track)
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
		// Skip track setup for PWM-only tracks
		Track_vars::iterator pwm_var = m_track_vars.begin();
		for (; pwm_var != m_track_vars.end(); ++pwm_var)
		{
			if (pwm_var->track_n_imdf == &(*itrack_n_imdf) &&
				(pwm_var->val_func == Track_var::PWM || pwm_var->val_func == Track_var::PWM_MAX || pwm_var->val_func == Track_var::PWM_MAX_POS || pwm_var->val_func == Track_var::KMER_COUNT || pwm_var->val_func == Track_var::KMER_FRAC))
				break;
		}
		if (pwm_var != m_track_vars.end())
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
        if (ivar->val_func == Track_var::PWM || ivar->val_func == Track_var::PWM_MAX || ivar->val_func == Track_var::PWM_MAX_POS) {
            pwm_vtracks.push_back(&*ivar);
        } else if (ivar->val_func == Track_var::KMER_COUNT || ivar->val_func == Track_var::KMER_FRAC) {
            kmer_vtracks.push_back(&*ivar);
        }
    }

    // Batch process if we have multiple sequence vtracks (threshold: 4+)
    if (kmer_vtracks.size() + pwm_vtracks.size() >= 4) {
        batch_process_sequence_vtracks(kmer_vtracks, pwm_vtracks, m_interval1d, idx);
    } else {
        // Process individually if too few
        for (Track_var* ivar : pwm_vtracks) {
            const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : m_interval1d;
            ivar->var[idx] = ivar->pwm_scorer->score_interval(seq_interval, m_iu.get_chromkey());
        }
        for (Track_var* ivar : kmer_vtracks) {
            const GInterval &seq_interval = ivar->seq_imdf1d ? ivar->seq_imdf1d->interval : m_interval1d;
            if (ivar->kmer_counter) {
                ivar->var[idx] = ivar->kmer_counter->score_interval(seq_interval, m_iu.get_chromkey());
            } else {
                ivar->var[idx] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    // Process regular track vtracks
    for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
        // Skip sequence-based vtracks (already processed above)
        if (ivar->val_func == Track_var::PWM || ivar->val_func == Track_var::PWM_MAX ||
            ivar->val_func == Track_var::PWM_MAX_POS || ivar->val_func == Track_var::KMER_COUNT ||
            ivar->val_func == Track_var::KMER_FRAC) {
            continue;
        }

		if (GenomeTrack::is_1d(ivar->track_n_imdf->type)) {
			GenomeTrack1D &track = *(GenomeTrack1D *)ivar->track_n_imdf->track;

			if (ivar->track_n_imdf->imdf1d && ivar->track_n_imdf->imdf1d->out_of_range)
				ivar->var[idx] = numeric_limits<double>::quiet_NaN();
			else {
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
				case Track_var::PWM_MAX:
				case Track_var::PWM_MAX_POS:
				case Track_var::PWM:
				case Track_var::KMER_COUNT:
				case Track_var::KMER_FRAC:
					break;
				default:
					verror("Internal error: unsupported function %d", ivar->val_func);
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
		}
		else if (ivar->val_func == Interv_var::COVERAGE)
		{
			const GInterval &interval = ivar->imdf1d ? ivar->imdf1d->interval : m_interval1d;

			if (ivar->imdf1d && ivar->imdf1d->out_of_range)
			{
				ivar->var[idx] = 0;
				continue;
			}

			int64_t total_overlap = 0;
			GIntervals::const_iterator iinterv;

			// For non-sequential access or first access
			if (ivar->imdf1d || ivar->siinterv == ivar->sintervs.end())
			{
				iinterv = lower_bound(ivar->sintervs.begin(), ivar->sintervs.end(), interval,
									  GIntervals::compare_by_start_coord);

				// Check previous interval too
				if (iinterv != ivar->sintervs.begin())
				{
					auto prev = iinterv - 1;
					if (prev->chromid == interval.chromid && prev->end > interval.start)
					{
						int64_t overlap_start = max(interval.start, prev->start);
						int64_t overlap_end = min(interval.end, prev->end);
						total_overlap += overlap_end - overlap_start;
					}
				}
			}
			else
			{
				// For sequential access, start from last position
				iinterv = ivar->siinterv;

				// Skip past intervals from previous chromosomes
				while (iinterv != ivar->sintervs.end() && iinterv->chromid < interval.chromid)
					++iinterv;

				// But check if we need to back up
				while (iinterv != ivar->sintervs.begin() &&
					   (iinterv - 1)->chromid == interval.chromid &&
					   (iinterv - 1)->end > interval.start)
				{
					--iinterv;
				}
			}

			// Check forward intervals
			while (iinterv != ivar->sintervs.end() &&
				   iinterv->chromid == interval.chromid &&
				   iinterv->start < interval.end)
			{
				if (iinterv->end > interval.start)
				{
					int64_t overlap_start = max(interval.start, iinterv->start);
					int64_t overlap_end = min(interval.end, iinterv->end);
					total_overlap += overlap_end - overlap_start;
				}
				++iinterv;
			}

			if (!ivar->imdf1d)
			{
				ivar->siinterv = iinterv;
			}

			ivar->var[idx] = (double)total_overlap / (interval.end - interval.start);
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
