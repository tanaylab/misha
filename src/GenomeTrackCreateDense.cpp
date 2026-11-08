#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "AggregationHelpers.h"
#include "GenomeTrackFixedBin.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrack_create_dense(SEXP _track, SEXP _data_frame, SEXP _binsize, SEXP _defvalue, SEXP _func, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track) || Rf_length(_track) != 1)
            verror("Track argument is not a string");

        if (!Rf_inherits(_data_frame, "data.frame")){
            verror("Data frame argument is not a data frame");
        }

        if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
            verror("Binsize argument is not a number");

        if ((!Rf_isReal(_defvalue) && !Rf_isInteger(_defvalue)) || Rf_length(_defvalue) != 1)
            verror("Defvalue argument is not a number");

        if (!Rf_isString(_func) || Rf_length(_func) != 1)
            verror("func argument is not a string");

        const char *track = CHAR(STRING_ELT(_track, 0));
        double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
        unsigned binsize = (unsigned)dbinsize;
        double defvalue = Rf_isReal(_defvalue) ? REAL(_defvalue)[0] : INTEGER(_defvalue)[0];
        const char *func_str = CHAR(STRING_ELT(_func, 0));

        if (dbinsize <= 0 || binsize != dbinsize)
            verror("Invalid value of binsize argument: %g\n", dbinsize);

        AggregationType agg_type = parse_aggregation_type(func_str);
        AggregationConfig agg_cfg{ agg_type, /*na_rm=*/true, /*min_n=*/-1, /*nth_index=*/-1 };

        // Extract columns from the data frame
        SEXP chrom_col = VECTOR_ELT(_data_frame, 0);
        SEXP start_col = VECTOR_ELT(_data_frame, 1);
        SEXP end_col = VECTOR_ELT(_data_frame, 2);
        SEXP value_col = VECTOR_ELT(_data_frame, 3);

        if (!Rf_isString(chrom_col))
            verror("Chromosome column is not a character vector");

        if ((!Rf_isReal(start_col) && !Rf_isInteger(start_col)))
            verror("Start column is not numeric");

        if ((!Rf_isReal(end_col) && !Rf_isInteger(end_col)))
            verror("End column is not numeric");

        if ((!Rf_isReal(value_col) && !Rf_isInteger(value_col)))
            verror("Value column is not numeric");

        int num_rows = Rf_length(chrom_col);
        if (Rf_length(start_col) != num_rows || Rf_length(end_col) != num_rows || Rf_length(value_col) != num_rows)
            verror("Data frame columns have inconsistent lengths");

        string dirname = create_track_dir(_envir, track);
        IntervUtils iu(_envir);

        // Convert the data frame to intervals and values. convert_rintervs validates the intervals
        // (known chromosome, no NA, 0 <= start < end <= chromosome size) as gtrack.create_sparse does.
        GIntervals data;
        vector<float> values;

        iu.convert_rintervs(_data_frame, &data, NULL);
        for (const GInterval &interval : data) {
            uint64_t i = iu.get_orig_interv_idx(interval);
            values.push_back(Rf_isReal(value_col) ? REAL(value_col)[i] : INTEGER(value_col)[i]);
        }

        // Sort intervals
        auto less_interv = [](const GInterval &a, const GInterval &b) {
            // Compare chromid first
            if (a.chromid != b.chromid)
                return a.chromid < b.chromid;
            // Then start coordinate
            if (a.start != b.start)
                return a.start < b.start;
            // Then end coordinate
            return a.end < b.end;
        };

        GIntervals sorted_data;
        vector<float> sorted_values;
        // Intervals from gextract or giterator.intervals come sorted and distinct: no sort then.
        // Equal intervals still go through the sort, whose order of them is the order their values
        // are summed in.
        if (adjacent_find(data.begin(), data.end(),
                [&less_interv](const GInterval &a, const GInterval &b) { return !less_interv(a, b); }) == data.end()) {
            sorted_data.swap(data);
            sorted_values.swap(values);
        } else {
            vector<int> sorted_idx(data.size());
            for (size_t i = 0; i < sorted_idx.size(); ++i)
                sorted_idx[i] = i;

            sort(sorted_idx.begin(), sorted_idx.end(),
                [&data, &less_interv](int i1, int i2) { return less_interv(data[i1], data[i2]); });

            sorted_data.reserve(data.size());
            sorted_values.reserve(data.size());
            for (int idx : sorted_idx) {
                sorted_data.push_back(data[idx]);
                sorted_values.push_back(values[idx]);
            }
        }

        // Get all genome intervals
        GIntervals all_genome_intervs;
        iu.get_all_genome_intervs(all_genome_intervs);

        Progress_reporter progress;
        progress.init(iu.get_chromkey().get_num_chroms(), 1);

        // The first interval that has not reached a bin yet
        size_t data_idx = 0;
        // The intervals that overlap the current bin, in sorted order. An interval stays here for every bin it
        // overlaps, so a long interval does not make the bins after it rescan the short intervals it spans.
        vector<size_t> active;

        // Create track files for each chromosome
        for (int chromid = 0; chromid < (int)iu.get_chromkey().get_num_chroms(); ++chromid) {
            uint64_t chromsize = iu.get_chromkey().get_chrom_size(chromid);
            char filename[FILENAME_MAX];

            snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(),
                     GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());

            GenomeTrackFixedBin gtrack;
            gtrack.init_write(filename, binsize, chromid);

            // Batch writing: accumulate bins and write in batches
            vector<float> bin_batch;
            bin_batch.reserve(10000);
            uint64_t bins_since_progress = 0;

            // Skip intervals from previous chromosomes
            while (data_idx < sorted_data.size() && sorted_data[data_idx].chromid < chromid) {
                data_idx++;
            }
            active.clear();

            // Fast path: the chromosome's intervals lie on the bin grid (one per bin at most, starting
            // at a multiple of binsize and ending at the bin's end), so a bin holds one contribution:
            // its interval's value, or defval where there is none or it is NaN. aggregate_single_value
            // computes what aggregate_values does for one contribution.
            size_t chrom_end_idx = data_idx;
            while (chrom_end_idx < sorted_data.size() && sorted_data[chrom_end_idx].chromid == chromid)
                ++chrom_end_idx;
            bool on_grid = has_single_value_form(agg_type);
            for (size_t i = data_idx; on_grid && i < chrom_end_idx; ++i) {
                const GInterval &interv = sorted_data[i];
                if (interv.start < 0 || (uint64_t)interv.start % binsize ||
                    (uint64_t)interv.end != min((uint64_t)interv.start + binsize, chromsize) ||
                    (i > data_idx && sorted_data[i - 1].start == interv.start))
                    on_grid = false;
            }
            if (on_grid) {
                size_t i = data_idx;
                for (uint64_t start_coord = 0; start_coord < chromsize; ) {
                    for (; start_coord < chromsize && bin_batch.size() < 10000; start_coord += binsize) {
                        double len = (double)(min(start_coord + (uint64_t)binsize, chromsize) - start_coord);
                        double v = numeric_limits<double>::quiet_NaN();
                        if (i < chrom_end_idx && (uint64_t)sorted_data[i].start == start_coord)
                            v = (double)sorted_values[i++];
                        if (std::isnan(v))
                            v = defvalue;
                        double agg_value = std::isnan(v) ? v : aggregate_single_value(agg_type, v, len);
                        bin_batch.push_back(std::isnan(agg_value) ? numeric_limits<float>::quiet_NaN() : (float)agg_value);
                    }
                    gtrack.write_next_bins(bin_batch.data(), bin_batch.size());
                    bin_batch.clear();
                    progress.report(0);
                    check_interrupt();
                }
                data_idx = chrom_end_idx;
                progress.report(1);
                continue;
            }

            AggregationState state;

            // Process each bin in the chromosome
            for (uint64_t start_coord = 0; start_coord < chromsize; start_coord += binsize) {
                uint64_t end_coord = min(start_coord + (uint64_t)binsize, chromsize);
                uint64_t bin_size_actual = end_coord - start_coord;

                // Drop the intervals that ended before this bin, and add those that start in it
                size_t num_active = 0;
                for (size_t idx : active) {
                    if ((uint64_t)sorted_data[idx].end > start_coord)
                        active[num_active++] = idx;
                }
                active.resize(num_active);
                while (data_idx < sorted_data.size() &&
                       sorted_data[data_idx].chromid == chromid &&
                       (uint64_t)sorted_data[data_idx].start < end_coord) {
                    active.push_back(data_idx++);
                }

                state.reset();
                // Bases of the bin under at least one interval with a value. Intervals may overlap, so this is
                // the length of their union, not the sum of their overlaps.
                uint64_t covered_bases = 0;
                uint64_t covered_end = start_coord;
                int64_t chain_id = 0;

                for (size_t idx : active) {
                    uint64_t overlap_start = max(start_coord, (uint64_t)sorted_data[idx].start);
                    uint64_t overlap_end = min(end_coord, (uint64_t)sorted_data[idx].end);

                    if (!std::isnan((double)sorted_values[idx])) {
                        uint64_t overlap_size = overlap_end - overlap_start;
                        aggregation_state_add(state,
                                              (double)sorted_values[idx],
                                              (double)overlap_size,
                                              (double)bin_size_actual,
                                              (int64_t)overlap_start,
                                              (int64_t)overlap_end,
                                              chain_id++);
                        // overlap_start does not decrease along the sorted intervals
                        if (overlap_end > covered_end) {
                            covered_bases += overlap_end - max(overlap_start, covered_end);
                            covered_end = overlap_end;
                        }
                    }
                }

                // Synthetic "uncovered" contribution if defvalue is set.
                // Skip for COUNT so empty/partial bins report only real interval count.
                if (covered_bases < bin_size_actual &&
                    !std::isnan(defvalue) &&
                    agg_type != AggregationType::COUNT)
                {
                    uint64_t uncovered = bin_size_actual - covered_bases;
                    aggregation_state_add(state,
                                          defvalue,
                                          (double)uncovered,
                                          (double)bin_size_actual,
                                          (int64_t)start_coord,
                                          (int64_t)end_coord,
                                          chain_id++);
                }

                double agg_value = aggregate_values(agg_cfg, state);
                float bin_value;
                if (std::isnan(agg_value)) {
                    bin_value = numeric_limits<float>::quiet_NaN();
                } else {
                    bin_value = (float)agg_value;
                }
                bin_batch.push_back(bin_value);

                // Write in batches of 10000 bins
                if (bin_batch.size() >= 10000) {
                    gtrack.write_next_bins(bin_batch.data(), bin_batch.size());
                    bin_batch.clear();
                    progress.report(0);
                    check_interrupt();
                    bins_since_progress = 0;
                } else {
                    ++bins_since_progress;
                    // Report progress every 10000 bins even if batch not full
                    if (bins_since_progress >= 10000) {
                        progress.report(0);
                        check_interrupt();
                        bins_since_progress = 0;
                    }
                }
            }

            // Write remaining bins
            if (!bin_batch.empty()) {
                gtrack.write_next_bins(bin_batch.data(), bin_batch.size());
            }

            progress.report(1);
        }

        progress.report_last();
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    return R_NilValue;
}

}
