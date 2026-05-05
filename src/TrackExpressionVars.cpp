#include <cstdint>
#include <cmath>
#include <map>
#include <sys/stat.h>
#include <unistd.h>
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
#include "TrackVarProcessor.h"
#include "IntervVarProcessor.h"
#include "ValueVarProcessor.h"
#include "SequenceVarProcessor.h"
#include "GlmVarProcessor.h"

namespace {
// Parse direction parameter from R params list. Returns ABOVE or BELOW.
PWMEditDistanceScorer::Direction parse_direction_param(SEXP rparams, const std::string& vtrack) {
    using rdb::verror;
    PWMEditDistanceScorer::Direction direction = PWMEditDistanceScorer::Direction::ABOVE;
    if (Rf_isNewList(rparams)) {
        int dir_idx = TrackExprParams::findListElementIndex(rparams, "direction");
        if (dir_idx >= 0) {
            SEXP rdir = VECTOR_ELT(rparams, dir_idx);
            if (!Rf_isNull(rdir)) {
                if (Rf_isString(rdir) && Rf_length(rdir) == 1) {
                    std::string dir_str = CHAR(STRING_ELT(rdir, 0));
                    if (dir_str == "above") {
                        direction = PWMEditDistanceScorer::Direction::ABOVE;
                    } else if (dir_str == "below") {
                        direction = PWMEditDistanceScorer::Direction::BELOW;
                    } else {
                        verror("direction parameter must be \"above\" or \"below\" for vtrack %s", vtrack.c_str());
                    }
                } else {
                    verror("direction parameter must be a single string (\"above\" or \"below\") for vtrack %s", vtrack.c_str());
                }
            }
        }
    }
    return direction;
}
} // anonymous namespace

const char *TrackExpressionVars::Track_var::FUNC_NAMES[TrackExpressionVars::Track_var::NUM_FUNCS] = {
	"avg", "min", "max", "nearest", "stddev", "sum", "lse", "quantile",
	"global.percentile", "global.percentile.min", "global.percentile.max",
	"weighted.sum", "area", "pwm", "pwm.max", "pwm.max.pos", "pwm.count", "kmer.count", "kmer.frac",
    "masked.count", "masked.frac",
    "max.pos.abs", "max.pos.relative", "min.pos.abs", "min.pos.relative",
    "exists", "size", "sample", "sample.pos.abs", "sample.pos.relative",
    "first", "first.pos.abs", "first.pos.relative", "last", "last.pos.abs", "last.pos.relative",
    "pwm.edit_distance", "pwm.edit_distance.pos", "pwm.max.edit_distance",
    "pwm.edit_distance.lse", "pwm.edit_distance.lse.pos",
    "pwm.n_mutations",
    "glm.predict"};

const char *TrackExpressionVars::Interv_var::FUNC_NAMES[TrackExpressionVars::Interv_var::NUM_FUNCS] = { "distance", "distance.center", "distance.edge", "coverage", "neighbor.count" };

const char *TrackExpressionVars::Value_var::FUNC_NAMES[TrackExpressionVars::Value_var::NUM_FUNCS] = {
	"avg", "min", "max", "stddev", "sum", "lse", "quantile",
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

namespace {
shared_ptr<GenomeTrack> create_and_init_1d_track(const string &filename, int chromid, GenomeTrack::Type type)
{
	shared_ptr<GenomeTrack> track;

	if (type == GenomeTrack::FIXED_BIN) {
		auto t = make_shared<GenomeTrackFixedBin>();
		t->init_read(filename.c_str(), chromid);
		track = t;
	} else if (type == GenomeTrack::SPARSE) {
		auto t = make_shared<GenomeTrackSparse>();
		t->init_read(filename.c_str(), chromid);
		track = t;
	} else if (type == GenomeTrack::ARRAYS) {
		auto t = make_shared<GenomeTrackArrays>();
		t->init_read(filename.c_str(), chromid);
		track = t;
	}

	return track;
}

bool supports_shared_1d_backend(GenomeTrack::Type type)
{
	return type == GenomeTrack::FIXED_BIN || type == GenomeTrack::SPARSE;
}

bool is_shared_1d_backend_enabled()
{
	SEXP opt = Rf_GetOption1(Rf_install("gshared.track.backend"));
	if (Rf_isNull(opt))
		return false;
	return Rf_asLogical(opt) == 1;
}

shared_ptr<GenomeTrack> create_dependent_1d_track(GenomeTrack::Type type, const shared_ptr<GenomeTrack> &master_track)
{
	if (type == GenomeTrack::FIXED_BIN) {
		GenomeTrackFixedBin *master = dynamic_cast<GenomeTrackFixedBin *>(master_track.get());
		if (!master)
			return shared_ptr<GenomeTrack>();

		auto dependent = make_shared<GenomeTrackFixedBin>();
		dependent->set_master_obj(master);
		return dependent;
	}

	if (type == GenomeTrack::SPARSE) {
		GenomeTrackSparse *master = dynamic_cast<GenomeTrackSparse *>(master_track.get());
		if (!master)
			return shared_ptr<GenomeTrack>();

		auto dependent = make_shared<GenomeTrackSparse>();
		dependent->set_master_obj(master);
		return dependent;
	}

	return shared_ptr<GenomeTrack>();
}

}  // namespace


TrackExpressionVars::TrackExpressionVars(rdb::IntervUtils &iu) :
	m_iu(iu)
{
	m_imdfs1d.reserve(10000);
	m_imdfs2d.reserve(10000);
	m_track_n_imdfs.reserve(10000);
	m_groot = get_groot(m_iu.get_env());
	m_shared_seqfetch.set_seqdir(m_groot + "/seq");
	m_track_processor = std::make_unique<TrackVarProcessor>(iu);
	m_interv_processor = std::make_unique<IntervVarProcessor>(iu);
	m_value_processor = std::make_unique<ValueVarProcessor>(iu);
	m_sequence_processor = std::make_unique<SequenceVarProcessor>(iu, m_shared_seqfetch);
	m_glm_processor = std::make_unique<GlmVarProcessor>(iu);
}

TrackExpressionVars::~TrackExpressionVars()
{
	for (Track_vars::iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ivar++)
		runprotect(ivar->rvar);
	for (Interv_vars::iterator ivar = m_interv_vars.begin(); ivar != m_interv_vars.end(); ivar++)
		runprotect(ivar->rvar);
	for (Value_vars::iterator ivar = m_value_vars.begin(); ivar != m_value_vars.end(); ivar++)
		runprotect(ivar->rvar);
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
		if (ivar->val_func == Interv_var::DIST || ivar->val_func == Interv_var::DIST_EDGE || ivar->val_func == Interv_var::NEIGHBOR_COUNT)
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
	track_n_imdf.track.reset();
	track_n_imdf.type = track_type;
	track_n_imdf.slice = slice;
	track_n_imdf.slice_func = slice_func;
	track_n_imdf.slice_percentile = slice_percentile;
	track_n_imdf.imdf1d = pimdf1d;
	track_n_imdf.imdf2d = pimdf2d;
	return track_n_imdf;
}

// Template implementation for attaching filter to any variable type with a 'filter' member
// Works with Track_var, Interv_var, and Value_var which all have:
//   std::shared_ptr<Genome1DFilter> filter
template<typename VarType>
void TrackExpressionVars::attach_filter_to_var(SEXP rvtrack, const string &vtrack, VarType &var)
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

// Explicit template instantiations for the three variable types
template void TrackExpressionVars::attach_filter_to_var<TrackExpressionVars::Track_var>(SEXP, const string &, Track_var &);
template void TrackExpressionVars::attach_filter_to_var<TrackExpressionVars::Interv_var>(SEXP, const string &, Interv_var &);
template void TrackExpressionVars::attach_filter_to_var<TrackExpressionVars::Value_var>(SEXP, const string &, Value_var &);

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
		} else if (func == "masked.count" || func == "masked.frac") {
			// Create the Track_var without a Track_n_imdf
			m_track_vars.push_back(Track_var());
			Track_var &var = m_track_vars.back();
			var.var_name = vtrack;
			var.val_func = func == "masked.count" ? Track_var::MASKED_COUNT : Track_var::MASKED_FRAC;
			var.track_n_imdf = nullptr; // No track needed for masked counting
			var.seq_imdf1d = nullptr;
			var.percentile = numeric_limits<double>::quiet_NaN();

			// No parameters needed for masked counting (unlike PWM/KMER)
			// Just create the counter with the appropriate mode
			MaskedBpCounter::CountMode mode =
				func == "masked.count" ? MaskedBpCounter::COUNT : MaskedBpCounter::FRACTION;
			var.masked_counter = std::make_unique<MaskedBpCounter>(&m_shared_seqfetch, mode);

			// Parse optional iterator modifier (sshift/eshift) for sequence-based vtracks
			Iterator_modifier1D imdf1d;
			parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
			var.seq_imdf1d = add_imdf(imdf1d);

			var.requires_pv = false;

			// Attach filter if present
			attach_filter_to_var(rvtrack, vtrack, var);
			return;
		} else if (func == "pwm.edit_distance" ||
		           func == "pwm.edit_distance.pos" ||
		           func == "pwm.max.edit_distance") {
            // Create the Track_var without a Track_n_imdf
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            if (func == "pwm.edit_distance") {
                var.val_func = Track_var::PWM_EDIT_DISTANCE;
            } else if (func == "pwm.edit_distance.pos") {
                var.val_func = Track_var::PWM_EDIT_DISTANCE_POS;
            } else {
                var.val_func = Track_var::PWM_MAX_EDIT_DISTANCE;
            }
            var.track_n_imdf = nullptr;  // No track needed for PWM edit distance
            var.seq_imdf1d = nullptr;

            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);

            // Parse PWM parameters using helper struct
            TrackExprParams::PWMParams pwm_params = TrackExprParams::PWMParams::parse(rparams, vtrack);

            // Extract max_edits parameter (optional, default -1 for exact computation)
            int max_edits = -1;
            if (Rf_isNewList(rparams)) {
                int me_idx = findListElementIndex(rparams, "max_edits");
                if (me_idx >= 0) {
                    SEXP rmax_edits = VECTOR_ELT(rparams, me_idx);
                    if (!Rf_isNull(rmax_edits)) {
                        if (Rf_isInteger(rmax_edits)) {
                            max_edits = INTEGER(rmax_edits)[0];
                        } else if (Rf_isReal(rmax_edits)) {
                            max_edits = static_cast<int>(REAL(rmax_edits)[0]);
                        } else {
                            verror("max_edits parameter must be NULL or a positive integer for vtrack %s", vtrack.c_str());
                        }
                        if (max_edits < 1) {
                            verror("max_edits must be a positive integer or NULL for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract max_indels parameter (optional, default 0 = substitutions only)
            int max_indels = 0;
            if (Rf_isNewList(rparams)) {
                int mi_idx = findListElementIndex(rparams, "max_indels");
                if (mi_idx >= 0) {
                    SEXP rmax_indels = VECTOR_ELT(rparams, mi_idx);
                    if (!Rf_isNull(rmax_indels)) {
                        if (Rf_isInteger(rmax_indels)) {
                            max_indels = INTEGER(rmax_indels)[0];
                        } else if (Rf_isReal(rmax_indels)) {
                            max_indels = static_cast<int>(REAL(rmax_indels)[0]);
                        } else {
                            verror("max_indels parameter must be NULL or a non-negative integer for vtrack %s", vtrack.c_str());
                        }
                        if (max_indels < 0) {
                            verror("max_indels must be a non-negative integer or NULL for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract score.min parameter (optional, default NaN = no filter)
            float score_min = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int sm_idx = findListElementIndex(rparams, "score.min");
                if (sm_idx >= 0) {
                    SEXP rscore_min = VECTOR_ELT(rparams, sm_idx);
                    if (!Rf_isNull(rscore_min)) {
                        if (Rf_isReal(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(REAL(rscore_min)[0]);
                        } else if (Rf_isInteger(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(INTEGER(rscore_min)[0]);
                        } else {
                            verror("score.min parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract score.max parameter (optional, default NaN = no filter)
            float score_max = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int smx_idx = findListElementIndex(rparams, "score.max");
                if (smx_idx >= 0) {
                    SEXP rscore_max = VECTOR_ELT(rparams, smx_idx);
                    if (!Rf_isNull(rscore_max)) {
                        if (Rf_isReal(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(REAL(rscore_max)[0]);
                        } else if (Rf_isInteger(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(INTEGER(rscore_max)[0]);
                        } else {
                            verror("score.max parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            auto direction = parse_direction_param(rparams, vtrack);

            PWMEditDistanceScorer::Mode mode = PWMEditDistanceScorer::Mode::MIN_EDITS;
            if (var.val_func == Track_var::PWM_EDIT_DISTANCE_POS) {
                mode = PWMEditDistanceScorer::Mode::MIN_EDITS_POSITION;
            } else if (var.val_func == Track_var::PWM_MAX_EDIT_DISTANCE) {
                mode = PWMEditDistanceScorer::Mode::PWM_MAX_EDITS;
            }

            // Construct scorer with shared sequence fetcher for caching
            var.pwm_edit_distance_scorer = std::make_unique<PWMEditDistanceScorer>(
                pwm_params.core.pssm,
                &m_shared_seqfetch,
                static_cast<float>(pwm_params.core.score_thresh),
                max_edits,
                pwm_params.extend_flag,
                static_cast<char>(pwm_params.core.strand_mode),
                mode,
                score_min,
                max_indels,
                score_max,
                direction
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
        } else if (func == "pwm.edit_distance.lse" ||
                   func == "pwm.edit_distance.lse.pos") {
            // LSE-mode edit distance: minimum edits to raise LSE score above threshold
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            if (func == "pwm.edit_distance.lse") {
                var.val_func = Track_var::PWM_EDIT_DISTANCE_LSE;
            } else {
                var.val_func = Track_var::PWM_EDIT_DISTANCE_LSE_POS;
            }
            var.track_n_imdf = nullptr;
            var.seq_imdf1d = nullptr;

            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);

            // Parse PWM parameters using helper struct
            TrackExprParams::PWMParams pwm_params = TrackExprParams::PWMParams::parse(rparams, vtrack);

            // Extract max_edits parameter (optional, default -1 for unlimited)
            int max_edits = -1;
            if (Rf_isNewList(rparams)) {
                int me_idx = findListElementIndex(rparams, "max_edits");
                if (me_idx >= 0) {
                    SEXP rmax_edits = VECTOR_ELT(rparams, me_idx);
                    if (!Rf_isNull(rmax_edits)) {
                        if (Rf_isInteger(rmax_edits)) {
                            max_edits = INTEGER(rmax_edits)[0];
                        } else if (Rf_isReal(rmax_edits)) {
                            max_edits = static_cast<int>(REAL(rmax_edits)[0]);
                        } else {
                            verror("max_edits parameter must be NULL or a positive integer for vtrack %s", vtrack.c_str());
                        }
                        if (max_edits < 1) {
                            verror("max_edits must be a positive integer or NULL for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract score.min parameter (optional, default NaN = no filter)
            float score_min = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int sm_idx = findListElementIndex(rparams, "score.min");
                if (sm_idx >= 0) {
                    SEXP rscore_min = VECTOR_ELT(rparams, sm_idx);
                    if (!Rf_isNull(rscore_min)) {
                        if (Rf_isReal(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(REAL(rscore_min)[0]);
                        } else if (Rf_isInteger(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(INTEGER(rscore_min)[0]);
                        } else {
                            verror("score.min parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract score.max parameter (optional, default NaN = no filter)
            float score_max = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int smx_idx = findListElementIndex(rparams, "score.max");
                if (smx_idx >= 0) {
                    SEXP rscore_max = VECTOR_ELT(rparams, smx_idx);
                    if (!Rf_isNull(rscore_max)) {
                        if (Rf_isReal(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(REAL(rscore_max)[0]);
                        } else if (Rf_isInteger(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(INTEGER(rscore_max)[0]);
                        } else {
                            verror("score.max parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            PWMLseEditDistanceScorer::Mode lse_mode = PWMLseEditDistanceScorer::Mode::LSE_EDIT_DISTANCE;
            if (var.val_func == Track_var::PWM_EDIT_DISTANCE_LSE_POS) {
                lse_mode = PWMLseEditDistanceScorer::Mode::LSE_EDIT_DISTANCE_POS;
            }

            auto lse_direction = parse_direction_param(rparams, vtrack);

            // Construct LSE scorer with shared sequence fetcher
            var.pwm_lse_edit_distance_scorer = std::make_unique<PWMLseEditDistanceScorer>(
                pwm_params.core.pssm,
                &m_shared_seqfetch,
                static_cast<float>(pwm_params.core.score_thresh),
                max_edits,
                pwm_params.extend_flag,
                static_cast<char>(pwm_params.core.strand_mode),
                lse_mode,
                score_min,
                score_max,
                lse_direction
            );

            // Parse optional iterator modifier
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);

            var.percentile = numeric_limits<double>::quiet_NaN();
            var.requires_pv = false;

            // Attach filter if present
            attach_filter_to_var(rvtrack, vtrack, var);
            return;
        } else if (func == "pwm.n_mutations") {
            // N_MUTATIONS mode: count single-base substitutions that independently cross threshold
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            var.val_func = Track_var::PWM_N_MUTATIONS;
            var.track_n_imdf = nullptr;
            var.seq_imdf1d = nullptr;

            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);

            // Parse PWM parameters using helper struct
            TrackExprParams::PWMParams pwm_params = TrackExprParams::PWMParams::parse(rparams, vtrack);

            // Extract score.min parameter (optional, default NaN = no filter)
            float score_min = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int sm_idx = findListElementIndex(rparams, "score.min");
                if (sm_idx >= 0) {
                    SEXP rscore_min = VECTOR_ELT(rparams, sm_idx);
                    if (!Rf_isNull(rscore_min)) {
                        if (Rf_isReal(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(REAL(rscore_min)[0]);
                        } else if (Rf_isInteger(rscore_min) && Rf_length(rscore_min) == 1) {
                            score_min = static_cast<float>(INTEGER(rscore_min)[0]);
                        } else {
                            verror("score.min parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            // Extract score.max parameter (optional, default NaN = no filter)
            float score_max = std::numeric_limits<float>::quiet_NaN();
            if (Rf_isNewList(rparams)) {
                int smx_idx = findListElementIndex(rparams, "score.max");
                if (smx_idx >= 0) {
                    SEXP rscore_max = VECTOR_ELT(rparams, smx_idx);
                    if (!Rf_isNull(rscore_max)) {
                        if (Rf_isReal(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(REAL(rscore_max)[0]);
                        } else if (Rf_isInteger(rscore_max) && Rf_length(rscore_max) == 1) {
                            score_max = static_cast<float>(INTEGER(rscore_max)[0]);
                        } else {
                            verror("score.max parameter must be NULL or a single numeric value for vtrack %s", vtrack.c_str());
                        }
                    }
                }
            }

            auto direction = parse_direction_param(rparams, vtrack);

            // Construct scorer with N_MUTATIONS mode, no max_edits or max_indels
            var.pwm_edit_distance_scorer = std::make_unique<PWMEditDistanceScorer>(
                pwm_params.core.pssm,
                &m_shared_seqfetch,
                static_cast<float>(pwm_params.core.score_thresh),
                -1,  // max_edits: not applicable for n_mutations
                pwm_params.extend_flag,
                static_cast<char>(pwm_params.core.strand_mode),
                PWMEditDistanceScorer::Mode::N_MUTATIONS,
                score_min,
                0,   // max_indels: not applicable for n_mutations
                score_max,
                direction
            );

            // Parse optional iterator modifier
            Iterator_modifier1D imdf1d;
            parse_imdf(rvtrack, vtrack, &imdf1d, NULL);
            var.seq_imdf1d = add_imdf(imdf1d);

            var.percentile = numeric_limits<double>::quiet_NaN();
            var.requires_pv = false;

            // Attach filter if present
            attach_filter_to_var(rvtrack, vtrack, var);
            return;
        } else if (func == "glm.predict") {
            // Fused GLM predictor virtual track (sourceless)
            m_track_vars.push_back(Track_var());
            Track_var &var = m_track_vars.back();
            var.var_name = vtrack;
            var.val_func = Track_var::GLM_PREDICT;
            var.track_n_imdf = nullptr;

            SEXP rparams = get_rvector_col(rvtrack, "params", vtrack.c_str(), false);
            if (Rf_isNull(rparams) || !Rf_isNewList(rparams))
                verror("Virtual track %s: glm.predict requires params list", vtrack.c_str());

            // --- tracks (character vector, N entries) ---
            int t_idx = findListElementIndex(rparams, "tracks");
            if (t_idx < 0)
                verror("Virtual track %s: 'tracks' parameter is required", vtrack.c_str());
            SEXP rtracks = VECTOR_ELT(rparams, t_idx);
            if (!Rf_isString(rtracks) || Rf_length(rtracks) < 1)
                verror("Virtual track %s: 'tracks' must be a non-empty character vector", vtrack.c_str());
            int N = Rf_length(rtracks);

            // --- inner_func (character vector, N) ---
            int if_idx = findListElementIndex(rparams, "inner_func");
            if (if_idx < 0)
                verror("Virtual track %s: 'inner_func' parameter is required", vtrack.c_str());
            SEXP rif = VECTOR_ELT(rparams, if_idx);
            if (!Rf_isString(rif) || Rf_length(rif) != N)
                verror("Virtual track %s: 'inner_func' must be a character vector of length %d", vtrack.c_str(), N);

            // --- weights (numeric vector N, or matrix N×K flattened column-major) ---
            int w_idx = findListElementIndex(rparams, "weights");
            if (w_idx < 0)
                verror("Virtual track %s: 'weights' parameter is required", vtrack.c_str());
            SEXP rweights = VECTOR_ELT(rparams, w_idx);
            if (!Rf_isReal(rweights))
                verror("Virtual track %s: 'weights' must be numeric", vtrack.c_str());

            // Detect K from num_bins parameter (default 1)
            int K = 1;
            int nb_idx = findListElementIndex(rparams, "num_bins");
            if (nb_idx >= 0) {
                SEXP rnb = VECTOR_ELT(rparams, nb_idx);
                if (!Rf_isNull(rnb)) {
                    if (Rf_isInteger(rnb)) K = INTEGER(rnb)[0];
                    else K = (int)Rf_asReal(rnb);
                }
            }
            var.glm_num_bins = K;

            if (Rf_length(rweights) != N * K)
                verror("Virtual track %s: 'weights' length must be %d (N=%d x K=%d)",
                       vtrack.c_str(), N * K, N, K);

            // --- bias (numeric scalar or vector of length K) ---
            var.glm_bias.assign(K, 0.0);
            int b_idx = findListElementIndex(rparams, "bias");
            if (b_idx >= 0) {
                SEXP rbias = VECTOR_ELT(rparams, b_idx);
                if (!Rf_isNull(rbias) && Rf_isReal(rbias)) {
                    double *bp = REAL(rbias);
                    int blen = Rf_length(rbias);
                    for (int i = 0; i < K && i < blen; i++)
                        var.glm_bias[i] = bp[i];
                }
            }

            // --- scale_factor (numeric scalar) ---
            var.glm_scale_factor = 10.0;
            int sf_idx = findListElementIndex(rparams, "scale_factor");
            if (sf_idx >= 0) {
                SEXP rsf = VECTOR_ELT(rparams, sf_idx);
                if (!Rf_isNull(rsf))
                    var.glm_scale_factor = Rf_asReal(rsf);
            }

            // --- sshifts / eshifts (numeric vectors, N) ---
            SEXP rsshifts = R_NilValue, reshifts = R_NilValue;
            int ss_idx = findListElementIndex(rparams, "sshifts");
            int es_idx = findListElementIndex(rparams, "eshifts");
            if (ss_idx >= 0) rsshifts = VECTOR_ELT(rparams, ss_idx);
            if (es_idx >= 0) reshifts = VECTOR_ELT(rparams, es_idx);

            // --- max_cap / dis_from_cap (numeric vectors, N, may contain NA) ---
            SEXP rmaxcap = R_NilValue, rdiscap = R_NilValue;
            int mc_idx = findListElementIndex(rparams, "max_cap");
            int dc_idx = findListElementIndex(rparams, "dis_from_cap");
            if (mc_idx >= 0) rmaxcap = VECTOR_ELT(rparams, mc_idx);
            if (dc_idx >= 0) rdiscap = VECTOR_ELT(rparams, dc_idx);

            // --- trans_family (character vector, N, may contain NA) ---
            SEXP rtfam = R_NilValue;
            int tf_idx = findListElementIndex(rparams, "trans_family");
            if (tf_idx >= 0) rtfam = VECTOR_ELT(rparams, tf_idx);

            // --- trans param vectors (5 parallel numeric vectors, N each) ---
            SEXP rtL = R_NilValue, rtk = R_NilValue, rtx0 = R_NilValue;
            SEXP rtpre = R_NilValue, rtpost = R_NilValue;
            int tL_idx = findListElementIndex(rparams, "trans_L");
            int tk_idx = findListElementIndex(rparams, "trans_k");
            int tx0_idx = findListElementIndex(rparams, "trans_x0");
            int tpre_idx = findListElementIndex(rparams, "trans_pre");
            int tpost_idx = findListElementIndex(rparams, "trans_post");
            if (tL_idx >= 0) rtL = VECTOR_ELT(rparams, tL_idx);
            if (tk_idx >= 0) rtk = VECTOR_ELT(rparams, tk_idx);
            if (tx0_idx >= 0) rtx0 = VECTOR_ELT(rparams, tx0_idx);
            if (tpre_idx >= 0) rtpre = VECTOR_ELT(rparams, tpre_idx);
            if (tpost_idx >= 0) rtpost = VECTOR_ELT(rparams, tpost_idx);

            // --- kernel_bins (optional numeric vector, B) ---
            SEXP rkbins = R_NilValue;
            int kb_idx = findListElementIndex(rparams, "kernel_bins");
            if (kb_idx >= 0) rkbins = VECTOR_ELT(rparams, kb_idx);
            if (!Rf_isNull(rkbins) && Rf_isReal(rkbins)) {
                int B = Rf_length(rkbins);
                var.glm_kernel_bins.resize(B);
                double *kbp = REAL(rkbins);
                for (int b = 0; b < B; b++)
                    var.glm_kernel_bins[b] = kbp[b];
            }

            // --- kernels (optional list of numeric vectors, N or recycled from length 1) ---
            SEXP rkernels = R_NilValue;
            int kn_idx = findListElementIndex(rparams, "kernels");
            if (kn_idx >= 0) rkernels = VECTOR_ELT(rparams, kn_idx);

            // --- Build GlmEntry vector ---
            var.glm_entries.resize(N);
            // Copy weights directly — R's column-major layout is already bin-major: wp[b*N + i]
            var.glm_all_weights.assign(REAL(rweights), REAL(rweights) + N * K);
            double *ssp = (!Rf_isNull(rsshifts) && Rf_isReal(rsshifts)) ? REAL(rsshifts) : nullptr;
            double *esp = (!Rf_isNull(reshifts) && Rf_isReal(reshifts)) ? REAL(reshifts) : nullptr;
            double *mcp = (!Rf_isNull(rmaxcap) && Rf_isReal(rmaxcap)) ? REAL(rmaxcap) : nullptr;
            double *dcp = (!Rf_isNull(rdiscap) && Rf_isReal(rdiscap)) ? REAL(rdiscap) : nullptr;
            double *tLp = (!Rf_isNull(rtL) && Rf_isReal(rtL)) ? REAL(rtL) : nullptr;
            double *tkp = (!Rf_isNull(rtk) && Rf_isReal(rtk)) ? REAL(rtk) : nullptr;
            double *tx0p = (!Rf_isNull(rtx0) && Rf_isReal(rtx0)) ? REAL(rtx0) : nullptr;
            double *tprep = (!Rf_isNull(rtpre) && Rf_isReal(rtpre)) ? REAL(rtpre) : nullptr;
            double *tpostp = (!Rf_isNull(rtpost) && Rf_isReal(rtpost)) ? REAL(rtpost) : nullptr;

            for (int i = 0; i < N; i++) {
                auto &entry = var.glm_entries[i];
                entry.track_name = CHAR(STRING_ELT(rtracks, i));

                // inner_func
                string ifunc = CHAR(STRING_ELT(rif, i));
                entry.inner_is_lse = (ifunc == "lse");

                // weight(s) — bin-major layout matches R column-major: wp[b*N + i]
                entry.weight_offset = i;

                // shifts
                entry.sshift = ssp ? (int64_t)ssp[i] : 0;
                entry.eshift = esp ? (int64_t)esp[i] : 0;

                // scaling
                if (mcp && !ISNA(mcp[i]) && dcp && !ISNA(dcp[i])) {
                    entry.scaling.enabled = true;
                    entry.scaling.max_cap = mcp[i];
                    entry.scaling.dis_from_cap = dcp[i];
                }

                // simple_cap flag (optional logical vector)
                int sc_idx = findListElementIndex(rparams, "simple_cap");
                if (sc_idx >= 0) {
                    SEXP rsc = VECTOR_ELT(rparams, sc_idx);
                    if (!Rf_isNull(rsc) && Rf_isLogical(rsc) && i < Rf_length(rsc)) {
                        entry.scaling.simple_cap = (LOGICAL(rsc)[i] == TRUE);
                    }
                }

                // transform
                if (rtfam != R_NilValue && Rf_isString(rtfam) && Rf_length(rtfam) > i) {
                    SEXP tfam_elt = STRING_ELT(rtfam, i);
                    if (tfam_elt != NA_STRING) {
                        string tfam = CHAR(tfam_elt);
                        if (tfam == "logist") {
                            entry.transform.enabled = true;
                            if (tLp && !ISNA(tLp[i])) entry.transform.L = tLp[i];
                            if (tkp && !ISNA(tkp[i])) entry.transform.k = tkp[i];
                            if (tx0p && !ISNA(tx0p[i])) entry.transform.x_0 = tx0p[i];
                            if (tprep && !ISNA(tprep[i])) entry.transform.pre_shift = tprep[i];
                            if (tpostp && !ISNA(tpostp[i])) entry.transform.post_shift = tpostp[i];
                        }
                    }
                }

                // kernel weights per entry (optional)
                if (!Rf_isNull(rkernels) && Rf_isNewList(rkernels)) {
                    int kn_len = Rf_length(rkernels);
                    int ki = (kn_len == 1) ? 0 : i; // recycle if length 1
                    SEXP rk = VECTOR_ELT(rkernels, ki);
                    if (Rf_isReal(rk)) {
                        int B = Rf_length(rk);
                        entry.kernel_weights.resize(B);
                        double *kp = REAL(rk);
                        for (int b = 0; b < B; b++)
                            entry.kernel_weights[b] = kp[b];
                    }
                }
            }

            // --- Interactions (optional) ---
            int ii_idx = findListElementIndex(rparams, "inter_i");
            int ij_idx = findListElementIndex(rparams, "inter_j");
            int iw_idx = findListElementIndex(rparams, "inter_weights");
            if (ii_idx >= 0 && ij_idx >= 0 && iw_idx >= 0) {
                SEXP rii = VECTOR_ELT(rparams, ii_idx);
                SEXP rij = VECTOR_ELT(rparams, ij_idx);
                SEXP riw = VECTOR_ELT(rparams, iw_idx);
                if (!Rf_isNull(rii) && !Rf_isNull(rij) && !Rf_isNull(riw)) {
                    int M = Rf_length(rii);  // inter_i always has exactly M elements
                    int iw_len = Rf_length(riw);
                    if (iw_len != M * K)
                        verror("Virtual track %s: 'inter_weights' length must be %d (M=%d x K=%d)",
                               vtrack.c_str(), M * K, M, K);
                    // inter_i/inter_j can be integer or real
                    var.glm_interactions.resize(M);
                    // Copy interaction weights directly — R column-major = bin-major: iwp[b*M + m]
                    var.glm_all_inter_weights.assign(REAL(riw), REAL(riw) + M * K);

                    // Get interaction transform params
                    SEXP ritfam = R_NilValue;
                    SEXP ritL = R_NilValue, ritk = R_NilValue, ritx0 = R_NilValue;
                    SEXP ritpre = R_NilValue, ritpost = R_NilValue;
                    int itf_idx = findListElementIndex(rparams, "inter_trans_family");
                    if (itf_idx >= 0) ritfam = VECTOR_ELT(rparams, itf_idx);
                    int itL_idx = findListElementIndex(rparams, "inter_trans_L");
                    int itk_idx = findListElementIndex(rparams, "inter_trans_k");
                    int itx0_idx = findListElementIndex(rparams, "inter_trans_x0");
                    int itpre_idx = findListElementIndex(rparams, "inter_trans_pre");
                    int itpost_idx = findListElementIndex(rparams, "inter_trans_post");
                    if (itL_idx >= 0) ritL = VECTOR_ELT(rparams, itL_idx);
                    if (itk_idx >= 0) ritk = VECTOR_ELT(rparams, itk_idx);
                    if (itx0_idx >= 0) ritx0 = VECTOR_ELT(rparams, itx0_idx);
                    if (itpre_idx >= 0) ritpre = VECTOR_ELT(rparams, itpre_idx);
                    if (itpost_idx >= 0) ritpost = VECTOR_ELT(rparams, itpost_idx);

                    double *iip = Rf_isInteger(rii) ? nullptr : REAL(rii);
                    int *iipi = Rf_isInteger(rii) ? INTEGER(rii) : nullptr;
                    double *ijp = Rf_isInteger(rij) ? nullptr : REAL(rij);
                    int *ijpi = Rf_isInteger(rij) ? INTEGER(rij) : nullptr;
                    double *itLp = (!Rf_isNull(ritL) && Rf_isReal(ritL)) ? REAL(ritL) : nullptr;
                    double *itkp = (!Rf_isNull(ritk) && Rf_isReal(ritk)) ? REAL(ritk) : nullptr;
                    double *itx0p = (!Rf_isNull(ritx0) && Rf_isReal(ritx0)) ? REAL(ritx0) : nullptr;
                    double *itprep = (!Rf_isNull(ritpre) && Rf_isReal(ritpre)) ? REAL(ritpre) : nullptr;
                    double *itpostp = (!Rf_isNull(ritpost) && Rf_isReal(ritpost)) ? REAL(ritpost) : nullptr;

                    for (int m = 0; m < M; m++) {
                        auto &inter = var.glm_interactions[m];
                        // 1-based R index to 0-based C++ index
                        inter.entry_i = (iipi ? iipi[m] : (int)iip[m]) - 1;
                        inter.entry_j = (ijpi ? ijpi[m] : (int)ijp[m]) - 1;
                        inter.weight_offset = m;  // bin-major: access as [b*M + m]

                        if (inter.entry_i < 0 || inter.entry_i >= N ||
                            inter.entry_j < 0 || inter.entry_j >= N)
                            verror("Virtual track %s: interaction index out of range [1, %d]",
                                   vtrack.c_str(), N);

                        // Interaction transform
                        if (ritfam != R_NilValue && Rf_isString(ritfam) && Rf_length(ritfam) > m) {
                            SEXP itfam_elt = STRING_ELT(ritfam, m);
                            if (itfam_elt != NA_STRING) {
                                string itfam = CHAR(itfam_elt);
                                if (itfam == "logist") {
                                    inter.transform.enabled = true;
                                    if (itLp && !ISNA(itLp[m])) inter.transform.L = itLp[m];
                                    if (itkp && !ISNA(itkp[m])) inter.transform.k = itkp[m];
                                    if (itx0p && !ISNA(itx0p[m])) inter.transform.x_0 = itx0p[m];
                                    if (itprep && !ISNA(itprep[m])) inter.transform.pre_shift = itprep[m];
                                    if (itpostp && !ISNA(itpostp[m])) inter.transform.post_shift = itpostp[m];
                                }
                            }
                        }
                    }
                }
            }

            // --- selector_tracks (CHARACTER vector len M) and selector_breaks (LIST of REAL vectors len M) ---
            int st_idx = findListElementIndex(rparams, "selector_tracks");
            int sb_idx = findListElementIndex(rparams, "selector_breaks");
            if (st_idx >= 0) {
                SEXP rst = VECTOR_ELT(rparams, st_idx);
                if (!Rf_isNull(rst)) {
                    if (!Rf_isString(rst))
                        verror("'selector_tracks' must be a character vector");
                    int M = Rf_length(rst);
                    if (M < 1)
                        verror("'selector_tracks' must have length >= 1 when present");
                    if (sb_idx < 0)
                        verror("'selector_breaks' is required when 'selector_tracks' is present");
                    SEXP rsb = VECTOR_ELT(rparams, sb_idx);
                    if (!Rf_isVectorList(rsb) || Rf_length(rsb) != M)
                        verror("'selector_breaks' must be a list of length %d (one per selector)", M);

                    var.glm_selector_track_names.resize(M);
                    var.glm_selector_binfinders.resize(M);
                    var.glm_selector_strides.assign(M, 0);

                    int64_t K_total = 1;
                    for (int m = 0; m < M; ++m) {
                        var.glm_selector_track_names[m] = CHAR(STRING_ELT(rst, m));
                        SEXP rbm = VECTOR_ELT(rsb, m);
                        if (!Rf_isReal(rbm) || Rf_length(rbm) < 2)
                            verror("'selector_breaks[[%d]]' must be a numeric vector of length >= 2", m + 1);
                        int nbr = Rf_length(rbm);
                        std::vector<double> breaks(REAL(rbm), REAL(rbm) + nbr);
                        var.glm_selector_binfinders[m].init(breaks, true, true);
                        int K_m = nbr - 1;
                        var.glm_selector_strides[m] = (size_t)K_total;
                        K_total *= K_m;
                    }
                    if (K_total != var.glm_num_bins)
                        verror("Compound bin count from selectors (%lld) does not match weights/bias (%d)",
                               (long long)K_total, var.glm_num_bins);
                }
            }

            // Allocate scaled cache if interactions exist
            if (!var.glm_interactions.empty()) {
                var.glm_scaled_cache.resize(N, 0.0);
            }

            // Build scaling groups: entries with identical
            // (track_name, sshift, eshift, inner_is_lse, scaling) share one
            // aggregate+scale computation. Only the transform+weight differ.
            if (var.glm_kernel_bins.empty()) {
                // Group key: (track_name, sshift, eshift, inner_is_lse, scaling_enabled, max_cap, dis_from_cap)
                // Use a map with string key for simplicity at parse time (runs once).
                map<string, vector<int>> group_map;
                for (int i = 0; i < N; i++) {
                    auto &e = var.glm_entries[i];
                    // Build a string key from the grouping fields
                    char keybuf[256];
                    snprintf(keybuf, sizeof(keybuf),
                             "%s|%lld|%lld|%d|%d|%.17g|%.17g|%d",
                             e.track_name.c_str(),
                             (long long)e.sshift, (long long)e.eshift,
                             (int)e.inner_is_lse,
                             (int)e.scaling.enabled,
                             e.scaling.enabled ? e.scaling.max_cap : 0.0,
                             e.scaling.enabled ? e.scaling.dis_from_cap : 0.0,
                             (int)e.scaling.simple_cap);
                    group_map[string(keybuf)].push_back(i);
                }
                var.glm_scaling_groups.reserve(group_map.size());
                for (auto &kv : group_map) {
                    Track_var::GlmScalingGroup g;
                    g.representative_idx = kv.second[0];
                    g.entry_indices = std::move(kv.second);
                    var.glm_scaling_groups.push_back(std::move(g));
                }

                // Build track-level groups: cluster scaling groups by track name
                // for cache-friendly traversal (read one track's data, process all groups).
                map<string, vector<int>> tg_map;
                for (int sg = 0; sg < (int)var.glm_scaling_groups.size(); sg++) {
                    auto &rep = var.glm_entries[var.glm_scaling_groups[sg].representative_idx];
                    tg_map[rep.track_name].push_back(sg);
                }
                var.glm_track_groups.reserve(tg_map.size());
                for (auto &kv : tg_map) {
                    Track_var::GlmTrackGroup tg;
                    tg.track_entry_idx = var.glm_scaling_groups[kv.second[0]].representative_idx;
                    tg.min_sshift = numeric_limits<int64_t>::max();
                    tg.max_eshift = numeric_limits<int64_t>::min();
                    for (int sg_idx : kv.second) {
                        auto &rep = var.glm_entries[var.glm_scaling_groups[sg_idx].representative_idx];
                        tg.min_sshift = min(tg.min_sshift, rep.sshift);
                        tg.max_eshift = max(tg.max_eshift, rep.eshift);
                    }
                    tg.scaling_group_indices = std::move(kv.second);
                    var.glm_track_groups.push_back(std::move(tg));
                }
            }

            var.percentile = numeric_limits<double>::quiet_NaN();
            var.requires_pv = false;
            return;
        }
	}

	// Only check source if not a sequence-based function
	SEXP rsrc = get_rvector_col(rvtrack, "src", vtrack.c_str(), false);
	if (Rf_isNull(rsrc))
	{
		verror("Source must be specified for non-sequence-based virtual tracks");
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
					(GenomeTrack::is_2d(track_type) && (ifunc == Track_var::REG_NEAREST || ifunc == Track_var::STDDEV || ifunc == Track_var::SUM || ifunc == Track_var::LSE || ifunc == Track_var::QUANTILE)) ||
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
	} else if (!strcmp(func.c_str(), Interv_var::FUNC_NAMES[Interv_var::DIST_EDGE])) {
		var.val_func = Interv_var::DIST_EDGE;

		if (!Rf_isNull(rparams))
			verror("Virtual track %s: function %s does not accept any parameters", vtrack.c_str(), func.c_str());

		var.dist_margin = 0.;
		var.sintervs.swap(intervs1d);
		var.sintervs.sort();
		var.eintervs = var.sintervs;
		var.eintervs.sort(GIntervals::compare_by_end_coord);
		// Overlapping intervals are allowed for edge-to-edge distance
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
		case Value_var::LSE:
			var.track->register_function(GenomeTrack1D::LSE);
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
		// Skip GLM predict variables - they manage their own track reads
		if (TrackExpressionVars::is_glm_function(ivar->val_func)) {
			continue;
		}
		GenomeTrack1D *track1d = GenomeTrack::is_1d(ivar->track_n_imdf->type) ? (GenomeTrack1D *)ivar->track_n_imdf->track.get() : NULL;
		GenomeTrack2D *track2d = GenomeTrack::is_2d(ivar->track_n_imdf->type) ? (GenomeTrack2D *)ivar->track_n_imdf->track.get() : NULL;

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
		case Track_var::LSE:
			if (!track1d)
				verror("vtrack function 'lse' can only be used on 1D tracks");
			track1d->register_function(GenomeTrack1D::LSE);
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
		// Sequence-based functions bypass normal track function registration
		default:
			if (!TrackExpressionVars::is_sequence_based_function((Track_var::Val_func)ivar->val_func))
				verror("Unrecognized virtual track function");
			break;
		}

			if (ivar->track_n_imdf->type == GenomeTrack::ARRAYS) {
				GenomeTrackArrays *track = (GenomeTrackArrays *)ivar->track_n_imdf->track.get();
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
	    // Skip iterator validation for sourceless variables (sequence-based, glm.predict)
        if (TrackExpressionVars::is_sequence_based_function(itrack_var->val_func)) {
            continue;
        }
        // GLM predict is 1D-only: it reads source tracks directly without dimension projection
        if (TrackExpressionVars::is_glm_function(itrack_var->val_func)) {
            if (expr_itr.is_2d())
                verror("Virtual track %s: glm.predict does not support 2D iterators", itrack_var->var_name.c_str());
            continue;
        }
        if (!itrack_var->track_n_imdf) {
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
		// Skip sourceless tracks (sequence-based, glm.predict)
		if (TrackExpressionVars::is_sequence_based_function(ivar->val_func))
		{
			continue;
		}
		if (!ivar->track_n_imdf) {
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
					GenomeTrackFixedBin gtrack_fbin;

					for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename)
					{
						int chromid = -1;

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
							// Metadata-only: this loop only validates bin_size
							// consistency, so skip the mmap setup that init_read
							// otherwise does. With many tracks × many
							// chromosomes the per-call open + mmap + madvise +
							// close was dominating gextract setup time.
							gtrack_fbin.init_read_metadata((trackpath + "/" + *ifilename).c_str(), chromid);
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

bool TrackExpressionVars::is_sequence_track_n_imdf(const Track_n_imdf &track_n_imdf) const
{
	for (Track_vars::const_iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ++ivar) {
		if (ivar->track_n_imdf == &track_n_imdf &&
			TrackExpressionVars::is_sequence_based_function(ivar->val_func))
			return true;
	}

	return false;
}

void TrackExpressionVars::reset_shared_1d_track_masters()
{
	m_shared_1d_track_masters.clear();
}

shared_ptr<GenomeTrack> TrackExpressionVars::init_1d_track_with_shared_backend(const string &filename, int chromid, GenomeTrack::Type type)
{
	if (!is_shared_1d_backend_enabled() || !supports_shared_1d_backend(type))
		return create_and_init_1d_track(filename, chromid, type);

	const BackendKey backend_key{type, chromid, filename};

	shared_ptr<GenomeTrack> master_track;
	auto imaster = m_shared_1d_track_masters.find(backend_key);
	if (imaster == m_shared_1d_track_masters.end()) {
		master_track = create_and_init_1d_track(filename, chromid, type);
		if (!master_track)
			return shared_ptr<GenomeTrack>();
		m_shared_1d_track_masters.emplace(backend_key, master_track);
	} else
		master_track = imaster->second;

	shared_ptr<GenomeTrack> dependent_track = create_dependent_1d_track(type, master_track);
	return dependent_track ? dependent_track : master_track;
}

void TrackExpressionVars::start_chrom(const GInterval &interval)
{
	reset_shared_1d_track_masters();

	// Precompute sequence-only flag (avoids O(N*M) scan in set_vars hot loop)
	for (Track_n_imdfs::iterator itn = m_track_n_imdfs.begin(); itn != m_track_n_imdfs.end(); ++itn) {
		itn->has_only_sequence_functions = is_sequence_track_n_imdf(*itn);
	}

	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		if (itrack_n_imdf->has_only_sequence_functions)
			continue;

		try
		{
			string track_dir = track2path(m_iu.get_env(), itrack_n_imdf->name);
			string resolved = GenomeTrack::find_existing_1d_filename(m_iu.get_chromkey(), track_dir, interval.chromid);
			string filename(track_dir + "/" + resolved);

			// For indexed tracks, reuse the persistent backend to avoid
			// re-mmapping the entire track.dat on every chromosome transition.
			const IndexedBackendKey idx_key{itrack_n_imdf->type, track_dir, 0};
			auto idx_it = m_indexed_track_backends.find(idx_key);

			if (idx_it != m_indexed_track_backends.end()) {
				// Reuse existing indexed track object — init_read reuses the mmap
				shared_ptr<GenomeTrack> master = idx_it->second;
				if (itrack_n_imdf->type == GenomeTrack::FIXED_BIN) {
					static_cast<GenomeTrackFixedBin *>(master.get())->init_read(filename.c_str(), interval.chromid);
				} else if (itrack_n_imdf->type == GenomeTrack::SPARSE) {
					static_cast<GenomeTrackSparse *>(master.get())->init_read(filename.c_str(), interval.chromid);
				}
				const BackendKey bk{itrack_n_imdf->type, interval.chromid, filename};
				m_shared_1d_track_masters.emplace(bk, master);
				shared_ptr<GenomeTrack> dependent = create_dependent_1d_track(itrack_n_imdf->type, master);
				itrack_n_imdf->track = dependent ? dependent : master;
			} else {
				shared_ptr<GenomeTrack> new_track = init_1d_track_with_shared_backend(filename, interval.chromid, itrack_n_imdf->type);
				if (!new_track) {
					verror("Internal error: track %s of type %s is not supported by 1D iterators",
						   itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
				}
				itrack_n_imdf->track = new_track;

				// If this is an indexed track, add to persistent cache
				struct stat idx_st;
				if (stat((track_dir + "/track.idx").c_str(), &idx_st) == 0 && supports_shared_1d_backend(itrack_n_imdf->type)) {
					const BackendKey bk{itrack_n_imdf->type, interval.chromid, filename};
					auto master_it = m_shared_1d_track_masters.find(bk);
					if (master_it != m_shared_1d_track_masters.end()) {
						m_indexed_track_backends.emplace(idx_key, master_it->second);
					} else {
						// Shared backends disabled — cache the track directly
						m_indexed_track_backends.emplace(idx_key, new_track);
					}
				}
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

	// Open track handles for GLM predict multi-track vtracks.
	// Many entries reference the same underlying track (e.g. 2430 entries
	// from 190 unique motifs). Cache handles by track_name so we open each
	// file only once, avoiding FD exhaustion under heavy parallelism.
	{
		struct GlmTrackCache {
			shared_ptr<GenomeTrack> handle;
			GenomeTrackFixedBin *fixedbin;
			GenomeTrackSparse *sparse;
			unsigned bin_size;
		};
		unordered_map<string, GlmTrackCache> glm_track_cache;

		for (auto &var : m_track_vars) {
			if (var.val_func != Track_var::GLM_PREDICT) continue;
			for (auto &entry : var.glm_entries) {
				try {
					auto it = glm_track_cache.find(entry.track_name);
					if (it != glm_track_cache.end()) {
						entry.track_handle   = it->second.handle;
						entry.cached_fixedbin = it->second.fixedbin;
						entry.cached_sparse   = it->second.sparse;
						entry.cached_bin_size = it->second.bin_size;
						continue;
					}

					string track_dir = track2path(m_iu.get_env(), entry.track_name);
					GenomeTrack::Type ttype = GenomeTrack::get_type(track_dir.c_str(), m_iu.get_chromkey(), false);
					string resolved = GenomeTrack::find_existing_1d_filename(m_iu.get_chromkey(), track_dir, interval.chromid);
					string filename(track_dir + "/" + resolved);
					// GLM reads bin_size/mmap directly without read_interval(),
					// so shared-backend dependents (which defer init) would be uninitialized.
					entry.track_handle = create_and_init_1d_track(filename, interval.chromid, ttype);
					if (entry.track_handle) {
						GenomeTrack *raw = entry.track_handle.get();
						entry.cached_fixedbin = dynamic_cast<GenomeTrackFixedBin*>(raw);
						entry.cached_sparse = entry.cached_fixedbin ? nullptr
						                      : dynamic_cast<GenomeTrackSparse*>(raw);
						if (entry.cached_fixedbin) {
							entry.cached_bin_size = entry.cached_fixedbin->get_bin_size();
						} else {
							entry.cached_bin_size = 0;
						}
					} else {
						entry.cached_fixedbin = nullptr;
						entry.cached_sparse = nullptr;
						entry.cached_bin_size = 0;
					}
					glm_track_cache.emplace(entry.track_name,
						GlmTrackCache{entry.track_handle, entry.cached_fixedbin,
						              entry.cached_sparse, entry.cached_bin_size});
				} catch (TGLException &) {
					entry.track_handle.reset();
					entry.cached_fixedbin = nullptr;
					entry.cached_sparse = nullptr;
					entry.cached_bin_size = 0;
				}
			}
		}
	}

	// Open selector track for GLM predict vars
	for (auto &var : m_track_vars) {
		if (var.val_func != Track_var::GLM_PREDICT) continue;
		if (var.glm_selector_track_name.empty()) continue;
		try {
			string track_dir = track2path(m_iu.get_env(), var.glm_selector_track_name);
			GenomeTrack::Type ttype = GenomeTrack::get_type(track_dir.c_str(), m_iu.get_chromkey(), false);
			string resolved = GenomeTrack::find_existing_1d_filename(m_iu.get_chromkey(), track_dir, interval.chromid);
			string filename(track_dir + "/" + resolved);
			var.glm_selector_handle = init_1d_track_with_shared_backend(filename, interval.chromid, ttype);
			if (!var.glm_selector_handle)
				verror("GLM selector track '%s' returned null handle for chrom %s",
				       var.glm_selector_track_name.c_str(),
				       m_iu.get_chromkey().id2chrom(interval.chromid).c_str());
			GenomeTrack *raw = var.glm_selector_handle.get();
			var.glm_selector_fixedbin = dynamic_cast<GenomeTrackFixedBin*>(raw);
			if (!var.glm_selector_fixedbin)
				verror("GLM selector track '%s' is not a fixed-bin track for chrom %s",
				       var.glm_selector_track_name.c_str(),
				       m_iu.get_chromkey().id2chrom(interval.chromid).c_str());
			var.glm_selector_bin_size = var.glm_selector_fixedbin->get_bin_size();
		} catch (TGLException &e) {
			verror("GLM selector track '%s' failed to open for chrom %s: %s",
			       var.glm_selector_track_name.c_str(),
			       m_iu.get_chromkey().id2chrom(interval.chromid).c_str(),
			       e.msg());
		}
	}

	// Pre-build GLM var batch list for the processor
	m_glm_processor->prepare_batch(m_track_vars);
}

void TrackExpressionVars::start_chrom(const GInterval2D &interval)
{
	unordered_map<BackendKey, shared_ptr<GenomeTrack>, BackendKeyHash> prev_shared_1d_track_masters;
	prev_shared_1d_track_masters.swap(m_shared_1d_track_masters);
	const bool shared_1d_backend_enabled = is_shared_1d_backend_enabled();

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

				const bool track_supports_shared_backend = supports_shared_1d_backend(itrack_n_imdf->type);
				const bool should_manage_shared_masters = track_supports_shared_backend &&
					(shared_1d_backend_enabled || !prev_shared_1d_track_masters.empty() || !m_shared_1d_track_masters.empty());

				string track_dir;
				string resolved;
				string filename;
				if (chromid != itrack_n_imdf->imdf1d->interval.chromid || should_manage_shared_masters) {
					track_dir = track2path(m_iu.get_env(), itrack_n_imdf->name);
					resolved = GenomeTrack::find_existing_1d_filename(m_iu.get_chromkey(), track_dir, chromid);
					filename = track_dir + "/" + resolved;
				}

				if (chromid != itrack_n_imdf->imdf1d->interval.chromid) {
					// For indexed tracks, reuse the persistent backend to avoid re-mmapping.
					// Key includes dim so DIM1/DIM2 projections get independent masters.
					const int dim = static_cast<int>(itrack_n_imdf->imdf1d->dim);
					const IndexedBackendKey idx_key{itrack_n_imdf->type, track_dir, dim};
					auto idx_it = m_indexed_track_backends.find(idx_key);

					if (idx_it != m_indexed_track_backends.end()) {
						shared_ptr<GenomeTrack> master = idx_it->second;
						if (itrack_n_imdf->type == GenomeTrack::FIXED_BIN)
							static_cast<GenomeTrackFixedBin *>(master.get())->init_read(filename.c_str(), chromid);
						else if (itrack_n_imdf->type == GenomeTrack::SPARSE)
							static_cast<GenomeTrackSparse *>(master.get())->init_read(filename.c_str(), chromid);
						const BackendKey bk{itrack_n_imdf->type, chromid, filename};
						m_shared_1d_track_masters.emplace(bk, master);
						shared_ptr<GenomeTrack> dependent = create_dependent_1d_track(itrack_n_imdf->type, master);
						itrack_n_imdf->track = dependent ? dependent : master;
					} else {
						shared_ptr<GenomeTrack> new_track = init_1d_track_with_shared_backend(filename, chromid, itrack_n_imdf->type);
						if (!new_track) {
							verror("Internal error: track %s of type %s is not supported by 1D iterators (projected from 2D)",
								   itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
						}
						itrack_n_imdf->track = new_track;

						// If indexed, add to persistent cache
						struct stat idx_st;
						if (stat((track_dir + "/track.idx").c_str(), &idx_st) == 0 && supports_shared_1d_backend(itrack_n_imdf->type)) {
							const BackendKey bk{itrack_n_imdf->type, chromid, filename};
							auto master_it = m_shared_1d_track_masters.find(bk);
							if (master_it != m_shared_1d_track_masters.end()) {
								m_indexed_track_backends.emplace(idx_key, master_it->second);
							} else {
								// Shared backends disabled — cache the track directly
								m_indexed_track_backends.emplace(idx_key, new_track);
							}
						}
					}
				} else if (should_manage_shared_masters) {
					const BackendKey backend_key{itrack_n_imdf->type, chromid, filename};
					auto imaster = m_shared_1d_track_masters.find(backend_key);

					if (imaster != m_shared_1d_track_masters.end()) {
						shared_ptr<GenomeTrack> new_track = create_dependent_1d_track(itrack_n_imdf->type, imaster->second);
						itrack_n_imdf->track = new_track ? new_track : imaster->second;
					} else {
						auto iprev_master = prev_shared_1d_track_masters.find(backend_key);
						if (iprev_master != prev_shared_1d_track_masters.end()) {
							m_shared_1d_track_masters.emplace(backend_key, iprev_master->second);
						} else {
							shared_ptr<GenomeTrack> new_track = init_1d_track_with_shared_backend(filename, chromid, itrack_n_imdf->type);
							if (!new_track) {
								verror("Internal error: track %s of type %s is not supported by 1D iterators (projected from 2D)",
									   itrack_n_imdf->name.c_str(), GenomeTrack::TYPE_NAMES[itrack_n_imdf->type]);
							}
							itrack_n_imdf->track = new_track;
						}
					}
				}
			} else if (!m_interval2d.is_same_chrom(interval)) {
				string filename(track2path(m_iu.get_env(), itrack_n_imdf->name) + "/" + GenomeTrack::get_2d_filename(m_iu.get_chromkey(), interval.chromid1(), interval.chromid2()));
				if (itrack_n_imdf->type == GenomeTrack::RECTS) {
					auto t = make_shared<GenomeTrackRectsRects>(m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					t->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
					itrack_n_imdf->track = t;
				} else if (itrack_n_imdf->type == GenomeTrack::POINTS) {
					auto t = make_shared<GenomeTrackRectsPoints>(m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					t->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
					itrack_n_imdf->track = t;
				} else if (itrack_n_imdf->type == GenomeTrack::COMPUTED) {
					auto t = make_shared<GenomeTrackComputed>(m_groot, m_iu.get_track_chunk_size(), m_iu.get_track_num_chunks());
					t->init_read(filename.c_str(), interval.chromid1(), interval.chromid2());
					itrack_n_imdf->track = t;
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
	// Setup tracks (read intervals for non-sequence tracks)
	for (Track_n_imdfs::iterator itrack_n_imdf = m_track_n_imdfs.begin(); itrack_n_imdf != m_track_n_imdfs.end(); ++itrack_n_imdf)
	{
		if (itrack_n_imdf->has_only_sequence_functions)
			continue;
		try {
				if (GenomeTrack::is_2d(itrack_n_imdf->type)) {
					if (itrack_n_imdf->imdf2d) {
						if (!itrack_n_imdf->imdf2d->out_of_range)
							((GenomeTrack2D *)itrack_n_imdf->track.get())->read_interval(itrack_n_imdf->imdf2d->interval, m_band);
					} else
						((GenomeTrack2D *)itrack_n_imdf->track.get())->read_interval(m_interval2d, m_band);
				} else {
					if (itrack_n_imdf->imdf1d) {
						if (!itrack_n_imdf->imdf1d->out_of_range)
							((GenomeTrack1D *)itrack_n_imdf->track.get())->read_interval(itrack_n_imdf->imdf1d->interval);
					} else
						((GenomeTrack1D *)itrack_n_imdf->track.get())->read_interval(m_interval1d);
				}
		} catch (TGLException &e) {
			verror("%s", e.msg());
		}
		}

	// Process sequence-based variables first (PWM, kmer, masked)
	m_sequence_processor->process_sequence_vars(m_track_vars, m_interval1d, idx);

	// Process regular track variables
	m_track_processor->process_track_vars(m_track_vars, m_interval1d, m_interval2d, m_band, idx);

	// Process interval variables
	m_interv_processor->process_interv_vars(m_interv_vars, m_interval1d, idx);

	// Process value variables
	m_value_processor->process_value_vars(m_value_vars, m_interval1d, idx);

	// Process GLM predict variables
	m_glm_processor->process_glm_vars(m_track_vars, m_interval1d, idx);
}
