#ifndef DATAFRAMEUTILS_H_
#define DATAFRAMEUTILS_H_

#include <R.h>
#include <Rinternals.h>
#include <vector>

// Static utility class for R data frame operations
class DataFrameUtils {
public:
	// Creates a data frame with given number or rows and columns. The data frame returned is still half baked.
	// Column names and the columns themselves must be defined later manually or via define_data_frame_cols.
	// If attrs_src is not R_NilValue, the attributes of the new data frame are copied from it.
	static SEXP create_data_frame(int numrows, int numcols, SEXP attrs_src = R_NilValue);

	// Copies columns definitions from src (must be a data frame) to tgt starting from column 'tgt_col_offset'.
	// tgt must be created by create_data_frame().
	// No column values are copied though. This function creates only the vectors of columns and copies column names.
	static void define_data_frame_cols(SEXP src, std::vector<SEXP> &src_cols, SEXP tgt, std::vector<SEXP> &tgt_cols, int tgt_col_offset);

	// Copies a row (values, not the definition) from src data frame to tgt data frame.
	// Before calling this function src columns must be defined in tgt by calling define_data_frame_cols().
	static void copy_data_frame_row(const std::vector<SEXP> &src_cols, int src_row, const std::vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset);

	// Copies multiple rows (values, not the definition) from src data frame to tgt data frame.
	// Before calling this function src columns must be defined in tgt by calling define_data_frame_cols().
	static void copy_data_frame_rows(const std::vector<SEXP> &src_cols, int src_row, int num_rows, const std::vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset);

	// Sets NAN at the given row and column of a data frame
	static void set_data_frame_val_nan(const std::vector<SEXP> &tgt_cols, int tgt_row, int tgt_col);
};

#endif /* DATAFRAMEUTILS_H_ */

