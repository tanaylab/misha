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

// The value of `_expr` on each row [start, end) of chromosome `chrom`, with
// base `_sub_base` (0..3 = A, C, G, T; -1 = none) at `_sub_pos` while the row
// is computed. Rows must be sorted by chromosome and start.
SEXP C_gextract_variants(SEXP _expr, SEXP _chroms, SEXP _starts, SEXP _ends, SEXP _sub_pos, SEXP _sub_base, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_expr) || Rf_length(_expr) != 1)
			verror("The track expression must be a single string");
		const R_xlen_t n = Rf_xlength(_chroms);
		if (!Rf_isString(_chroms) || !Rf_isReal(_starts) || !Rf_isReal(_ends) || !Rf_isReal(_sub_pos) ||
		    !Rf_isInteger(_sub_base) || Rf_xlength(_starts) != n || Rf_xlength(_ends) != n ||
		    Rf_xlength(_sub_pos) != n || Rf_xlength(_sub_base) != n)
			verror("Invalid rows: chroms, starts, ends, positions and bases must have the same length");

		static const char BASES[4] = {'A', 'C', 'G', 'T'};
		IntervUtils iu(_envir);
		GIntervals rows;
		rows.reserve(n);
		vector<SeqSubstitution> subs(n);
		for (R_xlen_t i = 0; i < n; ++i) {
			const int chromid = iu.chrom2id(CHAR(STRING_ELT(_chroms, i)));
			rows.push_back(GInterval(chromid, (int64_t)REAL(_starts)[i], (int64_t)REAL(_ends)[i], 0));
			if (i && (chromid < rows[i - 1].chromid ||
			          (chromid == rows[i - 1].chromid && rows[i].start < rows[i - 1].start)))
				verror("Rows must be sorted by chromosome and start");
			const int b = INTEGER(_sub_base)[i];
			if (b == NA_INTEGER || b < -1 || b > 3)
				verror("Invalid substituted base code %d", b);
			if (b >= 0) {
				subs[i].chromid = chromid;
				subs[i].pos = (int64_t)REAL(_sub_pos)[i];
				subs[i].base = BASES[b];
			}
		}

		GIntervals all_genome;
		iu.get_all_genome_intervs(all_genome);

		SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, n));
		double *out = REAL(answer);
		R_xlen_t i = 0;
		if (n) {
			TrackExprScanner scanner(iu);
			scanner.report_progress(false);
			scanner.set_rows(&rows, &subs);
			for (scanner.begin(_expr, &all_genome, NULL, R_NilValue, R_NilValue); !scanner.isend(); scanner.next())
				out[i++] = scanner.last_real(0);
		}
		if (i != n)
			verror("Internal error: %ld rows evaluated, %ld expected", (long)i, (long)n);
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

}
