#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

// Undefine R macros that conflict with C++ standard library
#undef length

#include "GenomeUtils.h"
#include "GIntervals.h"
#include "GInterval.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>

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

// Helper structure for valid genome segments
struct ValidSegment {
    int chromid;
    std::string chrom_name;
    int64_t start;
    int64_t end;
    int64_t length;
    double cum_prob;
};

// Generate random genome intervals
SEXP C_grandom_genome(SEXP _size, SEXP _n, SEXP _dist_from_edge, SEXP _chrom_df, SEXP _filter) {
    try {
        // Validate inputs
        if (!Rf_isInteger(_size) && !Rf_isReal(_size))
            Rf_error("size must be numeric");
        if (!Rf_isInteger(_n) && !Rf_isReal(_n))
            Rf_error("n must be numeric");
        if (!Rf_isReal(_dist_from_edge) && !Rf_isInteger(_dist_from_edge))
            Rf_error("dist_from_edge must be numeric");
        if (!Rf_inherits(_chrom_df, "data.frame"))        
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

        // Build chromosome name to id mapping
        std::map<std::string, int> chrom_name_to_id;
        std::vector<std::string> chrom_id_to_name(num_chroms);
        std::vector<int64_t> chrom_starts(num_chroms);
        std::vector<int64_t> chrom_ends(num_chroms);

        for (int i = 0; i < num_chroms; i++) {
            std::string chrom_name;
            if (Rf_isString(chroms)) {
                chrom_name = CHAR(STRING_ELT(chroms, i));
            } else if (Rf_isFactor(chroms)) {
                int level_idx = INTEGER(chroms)[i] - 1;
                if (level_idx < 0 || level_idx >= Rf_length(chrom_levels))
                    Rf_error("Invalid factor level in chrom column");
                chrom_name = CHAR(STRING_ELT(chrom_levels, level_idx));
            } else {
                Rf_error("chrom column must be character or factor");
            }

            chrom_name_to_id[chrom_name] = i;
            chrom_id_to_name[i] = chrom_name;
            chrom_starts[i] = Rf_isReal(starts) ? (int64_t)REAL(starts)[i] : (int64_t)INTEGER(starts)[i];
            chrom_ends[i] = Rf_isReal(ends) ? (int64_t)REAL(ends)[i] : (int64_t)INTEGER(ends)[i];
        }

        // Parse filter intervals if provided
        std::vector<ValidSegment> valid_segments;
        bool use_filter = !Rf_isNull(_filter) && Rf_isFrame(_filter) && Rf_length(VECTOR_ELT(_filter, 0)) > 0;

        if (use_filter) {
            // Parse filter data frame
            SEXP filter_chroms = VECTOR_ELT(_filter, 0);
            SEXP filter_starts = VECTOR_ELT(_filter, 1);
            SEXP filter_ends = VECTOR_ELT(_filter, 2);
            SEXP filter_chrom_levels = Rf_getAttrib(filter_chroms, R_LevelsSymbol);
            int num_filter = Rf_length(filter_chroms);

            // Build filter intervals per chromosome
            std::map<int, GIntervals> filter_by_chrom;

            for (int i = 0; i < num_filter; i++) {
                std::string filt_chrom_name;
                if (Rf_isString(filter_chroms)) {
                    filt_chrom_name = CHAR(STRING_ELT(filter_chroms, i));
                } else if (Rf_isFactor(filter_chroms)) {
                    int level_idx = INTEGER(filter_chroms)[i] - 1;
                    if (level_idx >= 0 && level_idx < Rf_length(filter_chrom_levels))
                        filt_chrom_name = CHAR(STRING_ELT(filter_chrom_levels, level_idx));
                }

                // Skip if chromosome not in our list
                if (chrom_name_to_id.find(filt_chrom_name) == chrom_name_to_id.end())
                    continue;

                int chromid = chrom_name_to_id[filt_chrom_name];
                int64_t filt_start = Rf_isReal(filter_starts) ? (int64_t)REAL(filter_starts)[i] : (int64_t)INTEGER(filter_starts)[i];
                int64_t filt_end = Rf_isReal(filter_ends) ? (int64_t)REAL(filter_ends)[i] : (int64_t)INTEGER(filter_ends)[i];

                // Expand filter by size to ensure no overlap
                filt_start = std::max((int64_t)0, filt_start - size + 1);

                filter_by_chrom[chromid].push_back(GInterval(chromid, filt_start, filt_end, 0));
            }

            // Compute valid segments for each chromosome
            for (int chromid = 0; chromid < num_chroms; chromid++) {
                int64_t chrom_start = chrom_starts[chromid];
                int64_t chrom_end = chrom_ends[chromid];

                // Apply dist_from_edge constraints
                int64_t valid_start = chrom_start + (int64_t)dist_from_edge;
                int64_t valid_end = chrom_end - (int64_t)dist_from_edge - size;

                if (valid_start >= valid_end)
                    continue; // Chromosome too short

                // Create full chromosome interval
                GIntervals chrom_interval;
                chrom_interval.push_back(GInterval(chromid, valid_start, valid_end, 0));

                GIntervals result;
                if (filter_by_chrom.find(chromid) != filter_by_chrom.end()) {
                    // Sort and unify filter intervals
                    filter_by_chrom[chromid].sort(GIntervals::compare_by_start_coord);
                    filter_by_chrom[chromid].unify_overlaps();

                    // Compute difference: chromosome - filter
                    GIntervals::diff(chrom_interval, filter_by_chrom[chromid], result);
                } else {
                    result = chrom_interval;
                }

                // Add valid segments
                for (size_t i = 0; i < result.size(); i++) {
                    int64_t seg_len = result[i].end - result[i].start;
                    if (seg_len >= size) {
                        ValidSegment seg;
                        seg.chromid = chromid;
                        seg.chrom_name = chrom_id_to_name[chromid];
                        seg.start = result[i].start;
                        seg.end = result[i].end;
                        seg.length = seg_len;
                        valid_segments.push_back(seg);
                    }
                }
            }

            if (valid_segments.empty()) {
                Rf_error("No valid regions available after applying filter");
            }

            // Calculate cumulative probabilities
            double total_length = 0.0;
            for (size_t i = 0; i < valid_segments.size(); i++) {
                total_length += valid_segments[i].length;
            }

            double cum_prob = 0.0;
            for (size_t i = 0; i < valid_segments.size(); i++) {
                cum_prob += valid_segments[i].length / total_length;
                valid_segments[i].cum_prob = cum_prob;
            }

        } else {
            // No filter: use simple chromosome-based sampling
            for (int chromid = 0; chromid < num_chroms; chromid++) {
                int64_t chrom_start = chrom_starts[chromid];
                int64_t chrom_end = chrom_ends[chromid];
                int64_t chrom_len = chrom_end - chrom_start;

                int64_t valid_start = chrom_start + (int64_t)dist_from_edge;
                int64_t valid_end = chrom_end - (int64_t)dist_from_edge - size;

                if (valid_start >= valid_end) {
                    Rf_error("Chromosome %s is too short (length: %lld, required: %lld)",
                             chrom_id_to_name[chromid].c_str(), (long long)chrom_len,
                             (long long)(size + 2 * dist_from_edge));
                }

                ValidSegment seg;
                seg.chromid = chromid;
                seg.chrom_name = chrom_id_to_name[chromid];
                seg.start = valid_start;
                seg.end = valid_end;
                seg.length = valid_end - valid_start;
                valid_segments.push_back(seg);
            }

            // Calculate cumulative probabilities
            double total_length = 0.0;
            for (size_t i = 0; i < valid_segments.size(); i++) {
                total_length += valid_segments[i].length;
            }

            double cum_prob = 0.0;
            for (size_t i = 0; i < valid_segments.size(); i++) {
                cum_prob += valid_segments[i].length / total_length;
                valid_segments[i].cum_prob = cum_prob;
            }
        }

        // Allocate output vectors
        SEXP out_chroms_idx, out_chroms_levels, out_starts, out_ends;
        SEXP out_df, col_names, row_names;

        PROTECT(out_df = Rf_allocVector(VECSXP, 3));
        PROTECT(out_chroms_idx = Rf_allocVector(INTSXP, n));
        PROTECT(out_starts = Rf_allocVector(REALSXP, n));
        PROTECT(out_ends = Rf_allocVector(REALSXP, n));
        PROTECT(out_chroms_levels = Rf_allocVector(STRSXP, num_chroms));
        PROTECT(col_names = Rf_allocVector(STRSXP, 3));
        PROTECT(row_names = Rf_allocVector(INTSXP, n));

        // Set chromosome levels
        for (int i = 0; i < num_chroms; i++) {
            SET_STRING_ELT(out_chroms_levels, i, Rf_mkChar(chrom_id_to_name[i].c_str()));
        }

        // Generate random intervals
        GetRNGstate();

        for (int64_t i = 0; i < n; i++) {
            // Sample segment using inverse transform sampling
            double u = unif_rand();
            size_t seg_idx = 0;
            for (size_t j = 0; j < valid_segments.size(); j++) {
                if (u <= valid_segments[j].cum_prob) {
                    seg_idx = j;
                    break;
                }
            }

            const ValidSegment &seg = valid_segments[seg_idx];

            // Generate random position within segment
            int64_t valid_range = seg.length - size;
            int64_t start = seg.start + (int64_t)(unif_rand() * valid_range);
            int64_t end = start + size;

            // Store in output vectors
            INTEGER(out_chroms_idx)[i] = seg.chromid + 1; // R is 1-indexed
            REAL(out_starts)[i] = (double)start;
            REAL(out_ends)[i] = (double)end;
            INTEGER(row_names)[i] = i + 1;
        }

        PutRNGstate();

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
        Rf_error("Error in gintervals.random: %s", e.what());
    } catch(...) {
        Rf_error("Unknown error in gintervals.random");
    }
    return R_NilValue; // Never reached
}

} 