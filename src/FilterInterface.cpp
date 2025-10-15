#include <cstdint>
#include <algorithm>
#include <memory>
#include <string>

#include "Filter.h"
#include "FilterRegistry.h"
#include "rdbinterval.h"
#include "rdbutils.h"
#include "TGLException.h"

using namespace std;
using namespace rdb;

extern "C" {

/**
 * Register a compiled filter from an R data frame.
 *
 * @param r_df Data frame with columns: chrom (character), start (numeric), end (numeric)
 *             Expected to be already canonicalized (sorted, merged overlaps)
 * @param r_key Unique string key for caching this filter
 * @param envir R environment
 * @return R_NilValue on success, error on failure
 */
SEXP C_register_filter(SEXP r_df, SEXP r_key, SEXP envir)
{
    try {
        RdbInitializer rdb_init;

        // Validate arguments
        if (!Rf_isString(r_key) || Rf_length(r_key) != 1)
            verror("Key argument must be a single string");

        if (!Rf_inherits(r_df, "data.frame"))
            verror("Filter argument must be a data frame");

        const char *key = CHAR(STRING_ELT(r_key, 0));

        // Check if filter already registered
        auto existing_filter = FilterRegistry::instance().get(key);
        if (existing_filter != nullptr) {
            return R_NilValue; // Already registered
        }

        IntervUtils iu(envir);

        // Convert R intervals to GIntervals
        GIntervals intervals;
        iu.convert_rintervs(r_df, &intervals, nullptr, false, &iu.get_chromkey(), "Filter: ");

        // Build Genome1DFilter
        auto filter = std::make_shared<Genome1DFilter>();

        if (intervals.empty()) {
            filter->empty = true;
            filter->per_chrom.resize(iu.get_chromkey().get_num_chroms());
        } else {
            filter->empty = false;

            // Size per_chrom vector based on chromkey
            int max_chromid = -1;
            for (const auto& interv : intervals) {
                if (interv.chromid > max_chromid) {
                    max_chromid = interv.chromid;
                }
            }
            filter->per_chrom.resize(max_chromid + 1);

            // Group intervals by chromosome
            // Intervals are already sorted by chromid, then start
            for (const auto& interv : intervals) {
                filter->per_chrom[interv.chromid].push_back(interv);
            }

            // Verify and merge overlaps per chromosome
            for (size_t chromid = 0; chromid < filter->per_chrom.size(); ++chromid) {
                auto& chrom_intervs = filter->per_chrom[chromid];
                if (chrom_intervs.empty()) continue;

                // Sort by start (should already be sorted, but ensure)
                std::sort(chrom_intervs.begin(), chrom_intervs.end(),
                         [](const GInterval& a, const GInterval& b) {
                             return a.start < b.start;
                         });

                // Merge overlapping/adjacent intervals
                size_t write_idx = 0;
                for (size_t read_idx = 1; read_idx < chrom_intervs.size(); ++read_idx) {
                    // If current interval overlaps or touches previous, extend previous
                    if (chrom_intervs[read_idx].start <= chrom_intervs[write_idx].end) {
                        chrom_intervs[write_idx].end = std::max(
                            chrom_intervs[write_idx].end,
                            chrom_intervs[read_idx].end
                        );
                    } else {
                        // Non-overlapping, advance write pointer
                        ++write_idx;
                        if (write_idx != read_idx) {
                            chrom_intervs[write_idx] = chrom_intervs[read_idx];
                        }
                    }
                }
                chrom_intervs.resize(write_idx + 1);
            }
        }

        // Register filter
        FilterRegistry::instance().put(key, filter);

        return R_NilValue;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

/**
 * Get filter statistics for display.
 *
 * @param r_key Filter key
 * @return Named list with: num_chroms, total_bases, empty
 */
SEXP C_get_filter_info(SEXP r_key)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(r_key) || Rf_length(r_key) != 1)
            verror("Key argument must be a single string");

        const char *key = CHAR(STRING_ELT(r_key, 0));
        auto filter = FilterRegistry::instance().get(key);

        if (filter == nullptr) {
            return R_NilValue;
        }

        SEXP result;
        SEXP names;

        rprotect(result = RSaneAllocVector(VECSXP, 3));
        rprotect(names = RSaneAllocVector(STRSXP, 3));

        SET_STRING_ELT(names, 0, Rf_mkChar("num_chroms"));
        SET_STRING_ELT(names, 1, Rf_mkChar("total_bases"));
        SET_STRING_ELT(names, 2, Rf_mkChar("empty"));

        SEXP num_chroms_sexp, total_bases_sexp, empty_sexp;
        rprotect(num_chroms_sexp = RSaneAllocVector(INTSXP, 1));
        rprotect(total_bases_sexp = RSaneAllocVector(REALSXP, 1));
        rprotect(empty_sexp = RSaneAllocVector(LGLSXP, 1));

        INTEGER(num_chroms_sexp)[0] = filter->num_masked_chroms();
        REAL(total_bases_sexp)[0] = filter->total_masked_bases();
        LOGICAL(empty_sexp)[0] = filter->empty;

        SET_VECTOR_ELT(result, 0, num_chroms_sexp);
        SET_VECTOR_ELT(result, 1, total_bases_sexp);
        SET_VECTOR_ELT(result, 2, empty_sexp);

        Rf_setAttrib(result, R_NamesSymbol, names);

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

/**
 * Remove a filter from the registry.
 *
 * @param r_key Filter key
 * @return TRUE if removed, FALSE if not found
 */
SEXP C_remove_filter(SEXP r_key)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(r_key) || Rf_length(r_key) != 1)
            verror("Key argument must be a single string");

        const char *key = CHAR(STRING_ELT(r_key, 0));
        bool removed = FilterRegistry::instance().remove(key);

        SEXP result;
        rprotect(result = RSaneAllocVector(LGLSXP, 1));
        LOGICAL(result)[0] = removed;

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

/**
 * Clear all filters from the registry.
 */
SEXP C_clear_filters()
{
    try {
        RdbInitializer rdb_init;
        FilterRegistry::instance().clear();
        return R_NilValue;
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
