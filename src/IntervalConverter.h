#ifndef INTERVALCONVERTER_H_
#define INTERVALCONVERTER_H_

#include <R.h>
#include <Rinternals.h>
#include <cstdint>

// Forward declarations
class GIntervals;
class GIntervals2D;
class GIntervalsFetcher1D;
class GIntervalsFetcher2D;
class GenomeChromKey;

namespace rdb {
	class IntervUtils;
}

// Converter class for interval conversions between R and C++
class IntervalConverter {
public:
	// Constructor takes a reference to IntervUtils for accessing chromosome key and environment
	explicit IntervalConverter(rdb::IntervUtils &iu);

	// Returns intervals type mask (INTERVS1D, INTERVS2D, or both)
	unsigned get_rintervs_type_mask(SEXP rintervals, const char *error_msg_prefix = "") const;

	// Converts R intervals to C++ intervals (returns type mask)
	// Returns: type mask (INTERVS1D, INTERVS2D, or both)
	unsigned convert_rintervs(SEXP rintervals, GIntervalsFetcher1D **intervals1d, GIntervalsFetcher2D **intervals2d, 
							  bool null_if_interv_nonexist = false, const GenomeChromKey *chromkey = NULL, 
							  const char *error_msg_prefix = "", bool verify = true) const;

	// Converts R intervals to C++ intervals (returns R data frame)
	// Returns: R data frame (or R_NilValue if null_if_interv_nonexist and interval doesn't exist)
	SEXP convert_rintervs(SEXP rintervals, GIntervals *intervals, GIntervals2D *intervals2d, 
						  bool null_if_interv_nonexist = false, const GenomeChromKey *chromkey = NULL, 
						  const char *error_msg_prefix = "", unsigned *pintervs_type_mask = NULL, 
						  bool verify = true, bool skip_missing_chroms = false) const;

	// Converts C++ 1D intervals to R data frame
	SEXP convert_intervs(GIntervalsFetcher1D *intervals, unsigned num_cols = 3, 
						 bool null_if_empty = true, bool use_original_index = false) const;

	// Converts C++ 2D intervals to R data frame
	SEXP convert_intervs(GIntervalsFetcher2D *intervals, unsigned num_cols = 6, 
						 bool null_if_empty = true, bool use_original_index = false) const;

private:
	rdb::IntervUtils &m_iu;
};

#endif /* INTERVALCONVERTER_H_ */

