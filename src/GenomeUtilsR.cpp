#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include "GenomeUtils.h"
#include <vector>
#include <string>
#include <map>

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

// R wrapper for seq2reverse
SEXP C_rev(SEXP _seq) {
    try {
        if (!Rf_isString(_seq))
            Rf_error("Sequence must be a character string");

        SEXP ans;
        int n = Rf_length(_seq);
        PROTECT(ans = Rf_allocVector(STRSXP, n));

        for (int i = 0; i < n; i++) {
            std::string seq = CHAR(STRING_ELT(_seq, i));
            std::string rev = seq2reverse(seq);
            SET_STRING_ELT(ans, i, Rf_mkChar(rev.c_str()));
        }

        UNPROTECT(1);
        return ans;
    } catch(std::exception &e) {
        Rf_error("Error in rev: %s", e.what());
    } catch(...) {
        Rf_error("Unknown error in rev");
    }
    return R_NilValue; // Never reached
}

// R wrapper for seq2complementary
SEXP C_comp(SEXP _seq) {
    try {
        if (!Rf_isString(_seq))
            Rf_error("Sequence must be a character string");

        SEXP ans;
        int n = Rf_length(_seq);
        PROTECT(ans = Rf_allocVector(STRSXP, n));

        for (int i = 0; i < n; i++) {
            std::string seq = CHAR(STRING_ELT(_seq, i));
            std::string comp = seq2complementary(seq);
            SET_STRING_ELT(ans, i, Rf_mkChar(comp.c_str()));
        }

        UNPROTECT(1);
        return ans;
    } catch(std::exception &e) {
        Rf_error("Error in comp: %s", e.what());
    } catch(...) {
        Rf_error("Unknown error in comp");
    }
    return R_NilValue; // Never reached
}

// Generate random genome intervals
SEXP C_grandom_genome(SEXP _size, SEXP _n, SEXP _dist_from_edge, SEXP _chrom_df) {
    try {
        // Validate inputs
        if (!Rf_isInteger(_size) && !Rf_isReal(_size))
            Rf_error("size must be numeric");
        if (!Rf_isInteger(_n) && !Rf_isReal(_n))
            Rf_error("n must be numeric");
        if (!Rf_isReal(_dist_from_edge) && !Rf_isInteger(_dist_from_edge))
            Rf_error("dist_from_edge must be numeric");
        if (!Rf_isFrame(_chrom_df))
            Rf_error("chrom_df must be a data frame");

        // Parse parameters
        int64_t size = Rf_isInteger(_size) ? INTEGER(_size)[0] : (int64_t)REAL(_size)[0];
        int64_t n = Rf_isInteger(_n) ? INTEGER(_n)[0] : (int64_t)REAL(_n)[0];
        double dist_from_edge = Rf_isReal(_dist_from_edge) ? REAL(_dist_from_edge)[0] : (double)INTEGER(_dist_from_edge)[0];

        if (size <= 0)
            Rf_error("size must be positive");
        if (n <= 0)
            Rf_error("n must be positive");
        if (dist_from_edge < 0)
            Rf_error("dist_from_edge must be non-negative");

        // Parse chromosome data frame
        SEXP chroms = VECTOR_ELT(_chrom_df, 0); // chrom column
        SEXP starts = VECTOR_ELT(_chrom_df, 1); // start column
        SEXP ends = VECTOR_ELT(_chrom_df, 2);   // end column
        SEXP chrom_levels = Rf_getAttrib(chroms, R_LevelsSymbol);

        int num_chroms = Rf_length(chroms);
        if (num_chroms == 0)
            Rf_error("No chromosomes provided");

        // Build chromosome info: name, length, and sampling probabilities
        struct ChromInfo {
            std::string name;
            int64_t length;
            double prob;
            double cum_prob;
        };

        std::vector<ChromInfo> chrom_info;
        chrom_info.reserve(num_chroms);

        double total_length = 0.0;

        for (int i = 0; i < num_chroms; i++) {
            ChromInfo info;

            // Get chromosome name
            if (Rf_isString(chroms)) {
                info.name = CHAR(STRING_ELT(chroms, i));
            } else if (Rf_isFactor(chroms)) {
                int level_idx = INTEGER(chroms)[i] - 1; // R factors are 1-indexed
                if (level_idx < 0 || level_idx >= Rf_length(chrom_levels))
                    Rf_error("Invalid factor level in chrom column");
                info.name = CHAR(STRING_ELT(chrom_levels, level_idx));
            } else {
                Rf_error("chrom column must be character or factor");
            }

            // Get chromosome length (end - start)
            int64_t start = Rf_isReal(starts) ? (int64_t)REAL(starts)[i] : (int64_t)INTEGER(starts)[i];
            int64_t end = Rf_isReal(ends) ? (int64_t)REAL(ends)[i] : (int64_t)INTEGER(ends)[i];
            info.length = end - start;

            // Validate chromosome length
            if (info.length < size + 2 * dist_from_edge) {
                Rf_error("Chromosome %s is too short (length: %lld, required: %lld)",
                         info.name.c_str(), (long long)info.length,
                         (long long)(size + 2 * dist_from_edge));
            }

            total_length += info.length;
            chrom_info.push_back(info);
        }

        // Calculate sampling probabilities and cumulative probabilities
        double cum_prob = 0.0;
        for (size_t i = 0; i < chrom_info.size(); i++) {
            chrom_info[i].prob = chrom_info[i].length / total_length;
            cum_prob += chrom_info[i].prob;
            chrom_info[i].cum_prob = cum_prob;
        }

        // Allocate output vectors
        SEXP out_chroms_idx, out_chroms_levels, out_starts, out_ends;
        SEXP out_df, col_names, row_names;

        PROTECT(out_df = Rf_allocVector(VECSXP, 3));
        PROTECT(out_chroms_idx = Rf_allocVector(INTSXP, n));
        PROTECT(out_starts = Rf_allocVector(REALSXP, n));
        PROTECT(out_ends = Rf_allocVector(REALSXP, n));
        PROTECT(out_chroms_levels = Rf_allocVector(STRSXP, chrom_info.size()));
        PROTECT(col_names = Rf_allocVector(STRSXP, 3));
        PROTECT(row_names = Rf_allocVector(INTSXP, n));

        // Set chromosome levels
        for (size_t i = 0; i < chrom_info.size(); i++) {
            SET_STRING_ELT(out_chroms_levels, i, Rf_mkChar(chrom_info[i].name.c_str()));
        }

        // Generate random intervals
        GetRNGstate(); // Get R's random number generator state

        for (int64_t i = 0; i < n; i++) {
            // Sample chromosome using inverse transform sampling
            double u = unif_rand();
            size_t chrom_idx = 0;
            for (size_t j = 0; j < chrom_info.size(); j++) {
                if (u <= chrom_info[j].cum_prob) {
                    chrom_idx = j;
                    break;
                }
            }

            // Generate random start position
            int64_t chrom_len = chrom_info[chrom_idx].length;
            int64_t valid_range = chrom_len - size - 2 * (int64_t)dist_from_edge;
            int64_t start = (int64_t)dist_from_edge + (int64_t)(unif_rand() * valid_range);
            int64_t end = start + size;

            // Store in output vectors (R is 1-indexed for factors)
            INTEGER(out_chroms_idx)[i] = chrom_idx + 1;
            REAL(out_starts)[i] = (double)start;
            REAL(out_ends)[i] = (double)end;
            INTEGER(row_names)[i] = i + 1;
        }

        PutRNGstate(); // Update R's random number generator state

        // Set up chromosome factor
        Rf_setAttrib(out_chroms_idx, R_LevelsSymbol, out_chroms_levels);
        Rf_setAttrib(out_chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

        // Build data frame
        SET_VECTOR_ELT(out_df, 0, out_chroms_idx);
        SET_VECTOR_ELT(out_df, 1, out_starts);
        SET_VECTOR_ELT(out_df, 2, out_ends);

        // Set column names
        SET_STRING_ELT(col_names, 0, Rf_mkChar("chrom"));
        SET_STRING_ELT(col_names, 1, Rf_mkChar("start"));
        SET_STRING_ELT(col_names, 2, Rf_mkChar("end"));

        // Set data frame attributes
        Rf_setAttrib(out_df, R_NamesSymbol, col_names);
        Rf_setAttrib(out_df, R_ClassSymbol, Rf_mkString("data.frame"));
        Rf_setAttrib(out_df, R_RowNamesSymbol, row_names);

        UNPROTECT(7);
        return out_df;

    } catch(std::exception &e) {
        Rf_error("Error in grandom_genome: %s", e.what());
    } catch(...) {
        Rf_error("Unknown error in grandom_genome");
    }
    return R_NilValue; // Never reached
}

} 