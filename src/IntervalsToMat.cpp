#include <R.h>
#include <Rinternals.h>
#include <cstdio>

// Build "chrom:start-end" strings for n intervals.
//   chrom: STRSXP or factor (INTSXP with levels)
//   start, end: INTSXP, length n
// Returns: STRSXP of length n.
extern "C" SEXP C_intervals_coord_strings(SEXP chrom_sxp, SEXP start_sxp, SEXP end_sxp) {
    if (TYPEOF(start_sxp) != INTSXP || TYPEOF(end_sxp) != INTSXP) {
        Rf_error("`start` and `end` must be integer vectors");
    }
    R_xlen_t n = XLENGTH(start_sxp);
    if (XLENGTH(end_sxp) != n) {
        Rf_error("`start` and `end` length mismatch");
    }

    int chrom_is_factor = Rf_isFactor(chrom_sxp);
    int chrom_is_char = (TYPEOF(chrom_sxp) == STRSXP);
    if (!chrom_is_factor && !chrom_is_char) {
        Rf_error("`chrom` must be character or factor");
    }
    if (!chrom_is_factor && XLENGTH(chrom_sxp) != n) {
        Rf_error("`chrom` length mismatch");
    }
    if (chrom_is_factor && XLENGTH(chrom_sxp) != n) {
        Rf_error("`chrom` length mismatch");
    }

    SEXP chrom_levels = R_NilValue;
    const int *chrom_ix = nullptr;
    if (chrom_is_factor) {
        chrom_levels = Rf_getAttrib(chrom_sxp, R_LevelsSymbol);
        chrom_ix = INTEGER(chrom_sxp);
    }

    const int *start = INTEGER(start_sxp);
    const int *end = INTEGER(end_sxp);

    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));

    // Chromosome names are short; coordinates fit in int32 (max ~2.1e9).
    // 64 bytes covers chrom (max ~40 chars in genomics) + ":" + 2 ints + "-".
    char buf[128];

    for (R_xlen_t i = 0; i < n; ++i) {
        const char *chrom_str;
        if (chrom_is_factor) {
            int lvl = chrom_ix[i];
            if (lvl == NA_INTEGER) {
                SET_STRING_ELT(out, i, NA_STRING);
                continue;
            }
            chrom_str = CHAR(STRING_ELT(chrom_levels, lvl - 1));
        } else {
            SEXP s = STRING_ELT(chrom_sxp, i);
            if (s == NA_STRING) {
                SET_STRING_ELT(out, i, NA_STRING);
                continue;
            }
            chrom_str = CHAR(s);
        }
        if (start[i] == NA_INTEGER || end[i] == NA_INTEGER) {
            SET_STRING_ELT(out, i, NA_STRING);
            continue;
        }
        std::snprintf(buf, sizeof(buf), "%s:%d-%d", chrom_str, start[i], end[i]);
        SET_STRING_ELT(out, i, Rf_mkChar(buf));
    }

    UNPROTECT(1);
    return out;
}
