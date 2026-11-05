#ifndef INTERVALCONVERTER_H_
#define INTERVALCONVERTER_H_

#include <R.h>
#include <Rinternals.h>

// Undefine R macros that conflict with C++ standard library
#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif
#ifdef warning
#undef warning
#endif

#include <cstddef>
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

// True for a chrom column the converters can read: character, a factor with character levels,
// or numbers (1 is the chromosome named "1").
bool is_chrom_col(SEXP col);

// The name of the chromosome in row i of a column accepted by is_chrom_col(). A number is printed
// into buf as a whole number; NA, a fraction or a bad factor code is a verror.
const char *chrom_name(SEXP chroms, SEXP chrom_levels, unsigned i, char *buf, size_t bufsize, const char *error_msg_prefix);

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

