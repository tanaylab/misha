// gintervals.pwm_hits: the windows a pwm.count virtual track counts, as intervals.

#include <cstdint>
#include <memory>
#include <vector>

#include "port.h"

#include "GenomeSeqFetch.h"
#include "PWMScorer.h"
#include "TrackExpressionParams.h"
#include "rdbinterval.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

extern "C" {

// _params: a pwm.count parameter list as .vtrack_params_pwm() builds it.
// Returns intervals (chrom, start, end) plus strand, score and intervalID (the
// 1-based row of the input interval the hit's anchor lies in).
SEXP C_gintervals_pwm_hits(SEXP _intervals, SEXP _params, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		IntervUtils iu(_envir);
		GIntervalsFetcher1D *intervals = NULL;
		iu.convert_rintervs(_intervals, &intervals, NULL);
		unique_ptr<GIntervalsFetcher1D> intervals_guard(intervals);
		intervals->sort();

		TrackExprParams::PWMParams params = TrackExprParams::PWMParams::parse(_params, "gintervals.pwm_hits");

		GenomeSeqFetch seqfetch;
		seqfetch.set_seqdir(string(get_groot(_envir)) + "/seq");

		PWMScorer scorer(params.core.pssm, &seqfetch, params.extend_flag, PWMScorer::MOTIF_COUNT,
		                 static_cast<char>(params.core.strand_mode), params.core.spat_factor,
		                 params.core.spat_bin_size, static_cast<float>(params.core.score_thresh));

		GIntervals out;
		vector<int> strands;
		vector<double> scores;
		vector<int> ids;
		vector<PWMScorer::Hit> hits;

		for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
			check_interrupt();
			const GInterval &interval = intervals->cur_interval();
			hits.clear();
			scorer.collect_hits(interval, iu.get_chromkey(), hits);
			const int id = (int)iu.get_orig_interv_idx(interval) + 1;
			for (const PWMScorer::Hit &h : hits) {
				out.push_back(GInterval(interval.chromid, h.start, h.end, (char)h.strand));
				strands.push_back(h.strand);
				scores.push_back(h.score);
				ids.push_back(id);
			}
			iu.verify_max_data_size(out.size(), "Result");
		}

		if (out.empty())
			return R_NilValue;

		const unsigned ncols = GInterval::NUM_COLS + 3;
		SEXP answer = iu.convert_intervs(&out, ncols, false);
		SEXP rstrand = rprotect_ptr(RSaneAllocVector(INTSXP, strands.size()));
		SEXP rscore = rprotect_ptr(RSaneAllocVector(REALSXP, scores.size()));
		SEXP rid = rprotect_ptr(RSaneAllocVector(INTSXP, ids.size()));
		for (size_t i = 0; i < strands.size(); ++i) {
			INTEGER(rstrand)[i] = strands[i];
			REAL(rscore)[i] = scores[i];
			INTEGER(rid)[i] = ids[i];
		}
		SET_VECTOR_ELT(answer, GInterval::NUM_COLS, rstrand);
		SET_VECTOR_ELT(answer, GInterval::NUM_COLS + 1, rscore);
		SET_VECTOR_ELT(answer, GInterval::NUM_COLS + 2, rid);
		SEXP colnames = Rf_getAttrib(answer, R_NamesSymbol);
		SET_STRING_ELT(colnames, GInterval::NUM_COLS, RSaneMkChar("strand"));
		SET_STRING_ELT(colnames, GInterval::NUM_COLS + 1, RSaneMkChar("score"));
		SET_STRING_ELT(colnames, GInterval::NUM_COLS + 2, RSaneMkChar("intervalID"));
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

}
