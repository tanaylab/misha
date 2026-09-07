#ifndef POTTS_PARAMS_H_
#define POTTS_PARAMS_H_

#include <cstring> // strcmp, in find_elt
#include <string>
#include <vector>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif
#ifdef warning
#undef warning
#endif

#include "PottsModel.h"
#include "rdbutils.h"

// The one SEXP -> PottsModel conversion, shared by C_gseq_potts() and
// TrackExpressionVars' potts branch. Two copies of the index arithmetic below
// would be two chances to transpose a coupling block, and the symptom of that
// is a plausible wrong number rather than a crash.
//
// `e`, `J`, `pairs` and `intercept` arrive pre-validated and pre-shaped from
// .coerce_potts_model() in R. This still checks types and shapes, because a
// .Call is a trust boundary, but it does not re-derive anything.
struct PottsParams {
    PottsModel model;
    bool bidirect = true;
    bool extend_flag = true;
    char strand_mode = 1;
    double score_thresh = 0.0;

    static int find_elt(SEXP list, const char *name) {
        SEXP names = Rf_getAttrib(list, R_NamesSymbol);
        if (names == R_NilValue)
            return -1;
        const int len = Rf_length(names);
        for (int i = 0; i < len; i++)
            if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
                return i;
        return -1;
    }

    static PottsParams parse(SEXP rparams, const std::string &who) {
        PottsParams params;

        if (!Rf_isNewList(rparams))
            rdb::verror("%s: potts parameters must be a named list", who.c_str());

        const int ie = find_elt(rparams, "e");
        if (ie < 0)
            rdb::verror("%s: potts parameters require an 'e' matrix", who.c_str());
        SEXP re = VECTOR_ELT(rparams, ie);
        if (!Rf_isMatrix(re) || !Rf_isReal(re))
            rdb::verror("%s: 'e' must be a numeric matrix", who.c_str());
        SEXP dim_e = Rf_getAttrib(re, R_DimSymbol);
        const int W = INTEGER(dim_e)[0];
        if (INTEGER(dim_e)[1] != 4)
            rdb::verror("%s: 'e' must have 4 columns (A, C, G, T)", who.c_str());
        if (W < 1)
            rdb::verror("%s: 'e' has no rows", who.c_str());

        const int ip = find_elt(rparams, "pairs");
        int npair = 0;
        SEXP rp = R_NilValue;
        if (ip >= 0) {
            rp = VECTOR_ELT(rparams, ip);
            if (rp != R_NilValue) {
                if (!Rf_isMatrix(rp) || !Rf_isInteger(rp))
                    rdb::verror("%s: 'pairs' must be an integer matrix", who.c_str());
                if (INTEGER(Rf_getAttrib(rp, R_DimSymbol))[1] != 2)
                    rdb::verror("%s: 'pairs' must have 2 columns", who.c_str());
                npair = INTEGER(Rf_getAttrib(rp, R_DimSymbol))[0];
            }
        }

        const int ij = find_elt(rparams, "J");
        SEXP rj = (ij >= 0) ? VECTOR_ELT(rparams, ij) : R_NilValue;
        if (npair) {
            if (rj == R_NilValue || !Rf_isMatrix(rj) || !Rf_isReal(rj) ||
                INTEGER(Rf_getAttrib(rj, R_DimSymbol))[0] != npair ||
                INTEGER(Rf_getAttrib(rj, R_DimSymbol))[1] != 16)
                rdb::verror("%s: 'J' must be a %d x 16 numeric matrix", who.c_str(), npair);
        }

        // R matrices are column-major. PottsModel wants `e` row-major
        // ([i*4 + base]) and each J block contiguous ([k*16 + b*4 + a]).
        std::vector<double> e_row((size_t)W * 4);
        const double *ep = REAL(re);
        for (int i = 0; i < W; ++i)
            for (int b = 0; b < 4; ++b)
                e_row[(size_t)i * 4 + b] = ep[i + (size_t)W * b];

        std::vector<int> p1(npair), p2(npair);
        if (npair) {
            const int *pp = INTEGER(rp);
            for (int k = 0; k < npair; ++k) {
                p1[k] = pp[k] - 1;                 // column 1, 1-based -> 0-based
                p2[k] = pp[k + (size_t)npair] - 1; // column 2
                if (p1[k] < 0 || p2[k] >= W || p1[k] >= p2[k])
                    rdb::verror("%s: 'pairs' row %d is out of 1..%d or not ascending",
                                who.c_str(), k + 1, W);
            }
        }

        std::vector<double> j_row((size_t)npair * 16);
        if (npair) {
            const double *jp = REAL(rj);
            for (int k = 0; k < npair; ++k)
                for (int f = 0; f < 16; ++f)
                    j_row[(size_t)k * 16 + f] = jp[k + (size_t)npair * f];
        }

        double intercept = 0.0;
        const int ii = find_elt(rparams, "intercept");
        if (ii >= 0 && VECTOR_ELT(rparams, ii) != R_NilValue) {
            SEXP ri = VECTOR_ELT(rparams, ii);
            if (!Rf_isReal(ri) || Rf_length(ri) != 1)
                rdb::verror("%s: 'intercept' must be a single number", who.c_str());
            intercept = REAL(ri)[0];
        }

        params.model = PottsModel(W, intercept, e_row.data(),
                                  npair ? j_row.data() : NULL, p1, p2);

        const int ib = find_elt(rparams, "bidirect");
        if (ib >= 0 && VECTOR_ELT(rparams, ib) != R_NilValue)
            params.bidirect = LOGICAL(VECTOR_ELT(rparams, ib))[0] == 1;

        const int ix = find_elt(rparams, "extend");
        if (ix >= 0 && VECTOR_ELT(rparams, ix) != R_NilValue)
            params.extend_flag = LOGICAL(VECTOR_ELT(rparams, ix))[0] == 1;

        const int is = find_elt(rparams, "strand");
        if (is >= 0 && VECTOR_ELT(rparams, is) != R_NilValue)
            params.strand_mode = (char)Rf_asInteger(VECTOR_ELT(rparams, is));

        const int it = find_elt(rparams, "score.thresh");
        if (it >= 0 && VECTOR_ELT(rparams, it) != R_NilValue)
            params.score_thresh = Rf_asReal(VECTOR_ELT(rparams, it));

        return params;
    }
};

#endif // POTTS_PARAMS_H_
