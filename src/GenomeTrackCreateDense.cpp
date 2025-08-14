#include <cstdint>
#include <cmath>
#include <cstring>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "GenomeTrackFixedBin.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrack_create_dense(SEXP _track, SEXP _data_frame, SEXP _binsize, SEXP _defvalue, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(_track) || Rf_length(_track) != 1)
            verror("Track argument is not a string");

        if (!Rf_isDataFrame(_data_frame))
            verror("Data frame argument is not a data frame");
        }

        if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
            verror("Binsize argument is not a number");

        if ((!Rf_isReal(_defvalue) && !Rf_isInteger(_defvalue)) || Rf_length(_defvalue) != 1)
            verror("Defvalue argument is not a number");

        const char *track = CHAR(STRING_ELT(_track, 0));
        double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
        unsigned binsize = (unsigned)dbinsize;
        double defvalue = Rf_isReal(_defvalue) ? REAL(_defvalue)[0] : INTEGER(_defvalue)[0];

        if (dbinsize <= 0 || binsize != dbinsize)
            verror("Invalid value of binsize argument: %g\n", dbinsize);

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

        // Process data and create intervals
        GIntervals data;
        vector<float> values;
        
        // Convert data frame to GIntervals and values
        for (int i = 0; i < num_rows; i++) {
            const char *chrom = CHAR(STRING_ELT(chrom_col, i));
            double start = Rf_isReal(start_col) ? REAL(start_col)[i] : INTEGER(start_col)[i];
            double end = Rf_isReal(end_col) ? REAL(end_col)[i] : INTEGER(end_col)[i];
            float value = Rf_isReal(value_col) ? REAL(value_col)[i] : INTEGER(value_col)[i];
            
            int chromid = iu.get_chromkey().chrom2id(chrom);
            if (chromid < 0)
                continue; // Skip unknown chromosomes
                
            GInterval interval;
            interval.chromid = chromid;
            interval.start = (int64_t)start;
            interval.end = (int64_t)end;
            data.push_back(interval);
            values.push_back(value);
        }
        
        // Sort intervals
        vector<int> sorted_idx(data.size());
        for (size_t i = 0; i < sorted_idx.size(); ++i)
            sorted_idx[i] = i;
            
        sort(sorted_idx.begin(), sorted_idx.end(), 
            [&data](int i1, int i2) { 
                // Compare chromid first
                if (data[i1].chromid != data[i2].chromid)
                    return data[i1].chromid < data[i2].chromid;
                // Then start coordinate
                if (data[i1].start != data[i2].start)
                    return data[i1].start < data[i2].start;
                // Then end coordinate
                return data[i1].end < data[i2].end;
            });
            
        GIntervals sorted_data;
        vector<float> sorted_values;
        for (int idx : sorted_idx) {
            sorted_data.push_back(data[idx]);
            sorted_values.push_back(values[idx]);
        }
        
        // Get all genome intervals
        GIntervals all_genome_intervs;
        iu.get_all_genome_intervs(all_genome_intervs);

        Progress_reporter progress;
        progress.init(iu.get_chromkey().get_num_chroms(), 1);

        // Initialize an index for data intervals
        size_t data_idx = 0;

        // Create track files for each chromosome
        for (int chromid = 0; chromid < (int)iu.get_chromkey().get_num_chroms(); ++chromid) {
            uint64_t chromsize = iu.get_chromkey().get_chrom_size(chromid);
            char filename[FILENAME_MAX];

            snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), 
                     GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());

            GenomeTrackFixedBin gtrack;
            gtrack.init_write(filename, binsize, chromid);

            // Process each bin in the chromosome
            for (uint64_t start_coord = 0; start_coord < chromsize; start_coord += binsize) {
                double sum = 0;
                uint64_t num_non_nans = 0;
                uint64_t end_coord = min(start_coord + binsize, chromsize);

                // Advance to the first interval that might overlap with the current bin
                while (data_idx < sorted_data.size() && 
                       (sorted_data[data_idx].chromid < chromid || 
                        (sorted_data[data_idx].chromid == chromid && (uint64_t)sorted_data[data_idx].end <= start_coord))) {
                    data_idx++;
                }

                // Process all intervals that overlap with the current bin
                size_t cur_idx = data_idx;
                while (cur_idx < sorted_data.size() && 
                       sorted_data[cur_idx].chromid == chromid && 
                       (uint64_t)sorted_data[cur_idx].start < end_coord) {
                    
                    // Calculate overlap between interval and current bin
                    uint64_t overlap_start = max(start_coord, (uint64_t)sorted_data[cur_idx].start);
                    uint64_t overlap_end = min(end_coord, (uint64_t)sorted_data[cur_idx].end);
                    uint64_t overlap_size = overlap_end - overlap_start;
                    
                    if (overlap_size > 0 && !std::isnan(sorted_values[cur_idx])) {
                        sum += sorted_values[cur_idx] * overlap_size;
                        num_non_nans += overlap_size;
                    }
                    
                    cur_idx++;
                }

                // If no overlapping intervals with values, use default value
                float bin_value;
                if (num_non_nans > 0) {
                    bin_value = sum / num_non_nans;
                } else {
                    bin_value = defvalue;
                }

                gtrack.write_next_bin(bin_value);
                progress.report(0);
                check_interrupt();
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