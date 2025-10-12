/*
 * TrackExpressionParams.h
 *
 * Parameter parsing structures for virtual track expressions
 * Consolidates parameter extraction logic for PWM, KMER, and other track types
 */

#ifndef TRACKEXPRESSIONPARAMS_H_
#define TRACKEXPRESSIONPARAMS_H_

#include <vector>
#include <string>
#include <limits>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "DnaPSSM.h"
#include "PWMScorer.h"
#include "rdbutils.h"

using namespace std;

namespace TrackExprParams {

// Helper to find list element by name
inline int findListElementIndex(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue)
        rdb::verror("List must have named elements");

    int len = Rf_length(list);
    for (int i = 0; i < len; i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
            return i;
    }
    return -1;  // Element not found
}

// PWM parameters structure
struct PWMParams {
    DnaPSSM pssm;
    bool bidirect;
    bool extend;
    char strand;
    std::vector<float> spat_factor;
    int spat_bin;
    int spat_min;
    int spat_max;
    bool has_range;
    float score_thresh;

    PWMParams()
        : bidirect(true), extend(true), strand(1),
          spat_bin(1), spat_min(0), spat_max(1000000),
          has_range(false), score_thresh(0.0f) {}

    // Parse PWM parameters from R list
    static PWMParams parse(SEXP rparams, const string& vtrack) {
        PWMParams params;

        if (!Rf_isNewList(rparams)) {
            rdb::verror("Virtual track %s: PWM functions require a list parameter with pssm matrix", vtrack.c_str());
        }

        // Get PSSM matrix (required)
        int pssm_idx = findListElementIndex(rparams, "pssm");
        if (pssm_idx < 0) {
            rdb::verror("Virtual track %s: PWM functions require a 'pssm' parameter", vtrack.c_str());
        }
        SEXP rpssm = VECTOR_ELT(rparams, pssm_idx);
        if (!Rf_isMatrix(rpssm)) {
            rdb::verror("Virtual track %s: PWM functions require a matrix parameter", vtrack.c_str());
        }
        params.pssm = PWMScorer::create_pssm_from_matrix(rpssm);

        // Get bidirect parameter (optional, default: true)
        int bidirect_idx = findListElementIndex(rparams, "bidirect");
        if (bidirect_idx >= 0) {
            SEXP rbidirect = VECTOR_ELT(rparams, bidirect_idx);
            if (rbidirect != R_NilValue) {
                if (!Rf_isLogical(rbidirect))
                    rdb::verror("Virtual track %s: bidirect parameter must be logical", vtrack.c_str());
                params.bidirect = LOGICAL(rbidirect)[0];
            }
        }
        params.pssm.set_bidirect(params.bidirect);

        // Get extend parameter (optional, default: true)
        int extend_idx = findListElementIndex(rparams, "extend");
        if (extend_idx >= 0) {
            SEXP rextend = VECTOR_ELT(rparams, extend_idx);
            if (rextend != R_NilValue) {
                if (!Rf_isLogical(rextend))
                    rdb::verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
                params.extend = LOGICAL(rextend)[0];
            }
        }

        // Get strand parameter (optional, default: 1)
        int strand_idx = findListElementIndex(rparams, "strand");
        if (strand_idx >= 0) {
            SEXP rstrand = VECTOR_ELT(rparams, strand_idx);
            if (rstrand != R_NilValue) {
                if (!Rf_isReal(rstrand) || Rf_length(rstrand) != 1)
                    rdb::verror("Virtual track %s: strand parameter must be numeric", vtrack.c_str());
                params.strand = (char)REAL(rstrand)[0];
            }
        }

        // Get spat_factor (optional)
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
            params.spat_factor.resize(n);
            for (int i = 0; i < n; ++i) {
                params.spat_factor[i] = REAL(rspat)[i];
                if (params.spat_factor[i] <= 0)
                    rdb::verror("Virtual track %s: all spat_factor values must be positive", vtrack.c_str());
            }

            // Get spat_bin (optional, default: 1)
            int bin_idx = findListElementIndex(rparams, "spat_bin");
            SEXP rbin = R_NilValue;
            if (bin_idx >= 0) {
                rbin = VECTOR_ELT(rparams, bin_idx);
            }
            if (rbin != R_NilValue) {
                if (!Rf_isInteger(rbin) && !Rf_isReal(rbin))
                    rdb::verror("Virtual track %s: spat_bin must be numeric", vtrack.c_str());
                params.spat_bin = (int)(Rf_isReal(rbin) ? REAL(rbin)[0] : INTEGER(rbin)[0]);
                if (params.spat_bin <= 0)
                    rdb::verror("Virtual track %s: spat_bin must be > 0", vtrack.c_str());
            }
        }

        // Get spat_min/spat_max (optional)
        int smin_idx = findListElementIndex(rparams, "spat_min");
        SEXP rsmin = R_NilValue;
        if (smin_idx >= 0) {
            rsmin = VECTOR_ELT(rparams, smin_idx);
        }
        if (rsmin != R_NilValue) {
            if (!Rf_isInteger(rsmin) && !Rf_isReal(rsmin))
                rdb::verror("Virtual track %s: spat_min must be numeric", vtrack.c_str());
            params.spat_min = (int)(Rf_isReal(rsmin) ? REAL(rsmin)[0] : INTEGER(rsmin)[0]);
            params.has_range = true;
        }

        int smax_idx = findListElementIndex(rparams, "spat_max");
        SEXP rsmax = R_NilValue;
        if (smax_idx >= 0) {
            rsmax = VECTOR_ELT(rparams, smax_idx);
        }
        if (rsmax != R_NilValue) {
            if (!Rf_isInteger(rsmax) && !Rf_isReal(rsmax))
                rdb::verror("Virtual track %s: spat_max must be numeric", vtrack.c_str());
            params.spat_max = (int)(Rf_isReal(rsmax) ? REAL(rsmax)[0] : INTEGER(rsmax)[0]);
            params.has_range = true;
        }

        // Apply range if provided
        if (params.has_range) {
            params.pssm.set_range(params.spat_min, params.spat_max);
        }

        // Get score_thresh parameter (for pwm.count)
        int thresh_idx = findListElementIndex(rparams, "score.thresh");
        SEXP rthresh = R_NilValue;
        if (thresh_idx >= 0) {
            rthresh = VECTOR_ELT(rparams, thresh_idx);
        }
        if (rthresh != R_NilValue) {
            if (!Rf_isReal(rthresh) || Rf_length(rthresh) != 1)
                rdb::verror("Virtual track %s: score.thresh parameter must be numeric", vtrack.c_str());
            params.score_thresh = REAL(rthresh)[0];
        }

        return params;
    }
};

// KMER parameters structure
struct KmerParams {
    std::string kmer;
    bool extend;
    char strand;

    KmerParams() : extend(true), strand(0) {}

    // Parse KMER parameters from R list or string
    static KmerParams parse(SEXP rparams, const string& vtrack) {
        KmerParams params;

        if (Rf_isNull(rparams)) {
            rdb::verror("Virtual track %s: kmer functions require a parameter (kmer string)", vtrack.c_str());
        }

        if (Rf_isNewList(rparams)) {
            // Handle as list params
            SEXP rextend = VECTOR_ELT(rparams, findListElementIndex(rparams, "extend"));
            if (rextend != R_NilValue) {
                if (!Rf_isLogical(rextend))
                    rdb::verror("Virtual track %s: extend parameter must be logical", vtrack.c_str());
                params.extend = LOGICAL(rextend)[0];
            }

            // Extract kmer string from the list parameters
            SEXP rkmer = VECTOR_ELT(rparams, findListElementIndex(rparams, "kmer"));
            if (rkmer == R_NilValue || !Rf_isString(rkmer) || Rf_length(rkmer) != 1)
                rdb::verror("Virtual track %s: invalid parameter used for kmer functions (must be a kmer string)", vtrack.c_str());
            params.kmer = CHAR(STRING_ELT(rkmer, 0));

            SEXP rstrand = VECTOR_ELT(rparams, findListElementIndex(rparams, "strand"));
            if (rstrand != R_NilValue) {
                if (!Rf_isNumeric(rstrand) || Rf_length(rstrand) != 1)
                    rdb::verror("Virtual track %s: strand parameter must be -1, 0, or 1", vtrack.c_str());
                params.strand = (char)REAL(rstrand)[0];
                if (params.strand != -1 && params.strand != 0 && params.strand != 1)
                    rdb::verror("Virtual track %s: strand parameter must be -1, 0, or 1", vtrack.c_str());
            }
        } else if (Rf_isString(rparams) && Rf_length(rparams) == 1) {
            // Handle direct string parameter (backward compatibility)
            params.kmer = CHAR(STRING_ELT(rparams, 0));
        } else {
            rdb::verror("Virtual track %s: invalid parameter used for kmer functions (must be a kmer string)", vtrack.c_str());
        }

        return params;
    }
};

} // namespace TrackExprParams

#endif /* TRACKEXPRESSIONPARAMS_H_ */
