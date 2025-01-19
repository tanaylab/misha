#include <R.h>
#include <Rinternals.h>
#include "GenomeUtils.h"

extern "C" {

// R wrapper for seq2reverse_complementary
SEXP C_revcomp(SEXP _seq) {
    try {
        if (!Rf_isString(_seq))
            Rf_error("Sequence must be a character string");

        SEXP ans;
        int n = Rf_length(_seq);
        PROTECT(ans = Rf_allocVector(STRSXP, n));

        for (int i = 0; i < n; i++) {
            std::string seq = CHAR(STRING_ELT(_seq, i));
            std::string revcomp = seq2reverse_complementary(seq);
            SET_STRING_ELT(ans, i, Rf_mkChar(revcomp.c_str()));
        }

        UNPROTECT(1);
        return ans;
    } catch(std::exception &e) {
        Rf_error("Error in revcomp: %s", e.what());
    } catch(...) {
        Rf_error("Unknown error in revcomp");
    }
    return R_NilValue; // Never reached
}

} 