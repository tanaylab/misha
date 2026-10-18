// gextract.variants: which bases each sequence vtrack in an expression reads,
// and the expression's value with one base substituted per row.

#include <cstdint>
#include <string>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionVars.h"

using namespace std;
using namespace rdb;

extern "C" {

// c(lo..., hi...) over the sequence vtracks of `_expr`; R reshapes it.
SEXP C_gvariants_reach(SEXP _expr, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_expr) || Rf_length(_expr) != 1)
			verror("The track expression must be a single string");

		IntervUtils iu(_envir);
		TrackExpressionVars vars(iu);
		vars.parse_exprs(vector<string>(1, CHAR(STRING_ELT(_expr, 0))));
		const vector<pair<int64_t, int64_t>> offs = vars.seq_read_offsets();

		SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, offs.size() * 2));
		for (size_t i = 0; i < offs.size(); ++i) {
			REAL(answer)[i] = (double)offs[i].first;
			REAL(answer)[i + offs.size()] = (double)offs[i].second;
		}
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

}
