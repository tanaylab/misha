#include "DataFrameUtils.h"
#include "rdbutils.h" // For rprotect_ptr, verror, RSaneAllocVector
#include <cstring> // For strcmp
#include <vector>

using namespace std;
using namespace rdb;

SEXP DataFrameUtils::create_data_frame(int numrows, int numcols, SEXP attrs_src)
{
	SEXP answer, row_names, col_names;

    answer = rprotect_ptr(RSaneAllocVector(VECSXP, numcols));
    col_names = rprotect_ptr(RSaneAllocVector(STRSXP, numcols));
    row_names = rprotect_ptr(RSaneAllocVector(INTSXP, numrows));

	for (int i = 0; i < numrows; ++i)
		INTEGER(row_names)[i] = i + 1;

    if (attrs_src != R_NilValue) 
        Rf_copyMostAttrib(attrs_src, answer);

    Rf_setAttrib(answer, R_NamesSymbol, col_names);
    Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
    Rf_setAttrib(answer, R_RowNamesSymbol, row_names);

	return answer;
}

void DataFrameUtils::define_data_frame_cols(SEXP src, vector<SEXP> &src_cols, SEXP tgt, vector<SEXP> &tgt_cols, int tgt_col_offset)
{
	SEXP src_class = Rf_getAttrib(src, R_ClassSymbol);

	if (Rf_isNull(src_class) || !Rf_isString(src_class) ||
		(!(Rf_length(src_class) == 1 && !strcmp(CHAR(STRING_ELT(src_class, 0)), "data.frame")) &&
		!(Rf_length(src_class) == 3 && !strcmp(CHAR(STRING_ELT(src_class, 0)), "tbl_df")  &&
          !strcmp(CHAR(STRING_ELT(src_class, 1)), "tbl") && !strcmp(CHAR(STRING_ELT(src_class, 2)), "data.frame"))))
		verror("Copied object is not a data frame or tibble");

	if (Rf_length(tgt) < Rf_length(src) + tgt_col_offset)
		verror("Attempt to copy data frame columns beyond the valid size");

	int numrows = Rf_length(Rf_getAttrib(tgt, R_RowNamesSymbol));
    SEXP src_colnames = rprotect_ptr(Rf_getAttrib(src, R_NamesSymbol));
    SEXP tgt_colnames = rprotect_ptr(Rf_getAttrib(tgt, R_NamesSymbol));

	if (Rf_isNull(src_colnames) || !Rf_isString(src_colnames))
		verror("Invalid source data frame for a copy");

	src_cols.resize(Rf_length(src));
	if (tgt_cols.size() < (uint64_t)(Rf_length(tgt) + tgt_col_offset)){ 
		tgt_cols.resize(Rf_length(tgt) + tgt_col_offset);
	}

	for (int col = 0; col < Rf_length(src); ++col) {
		SEXP src_col = VECTOR_ELT(src, col);
		SEXP tgt_col;

        tgt_col = rprotect_ptr(RSaneAllocVector(TYPEOF(src_col), numrows));

		if (!Rf_isInteger(src_col) && !Rf_isReal(src_col) && !Rf_isLogical(src_col) && !Rf_isString(src_col) && !Rf_isFactor(src_col))
			verror("Unsupported type found in a data frame: %s", Rf_type2char(TYPEOF(src_col)));

		Rf_copyMostAttrib(src_col, tgt_col);
		SET_STRING_ELT(tgt_colnames, col + tgt_col_offset, STRING_ELT(src_colnames, col));
		src_cols[col] = src_col;
		tgt_cols[col + tgt_col_offset] = tgt_col;

        SET_VECTOR_ELT(tgt, col + tgt_col_offset, tgt_col);
    }
}

void DataFrameUtils::copy_data_frame_row(const vector<SEXP> &src_cols, int src_row, const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset)
{
	for (uint64_t col = 0; col < src_cols.size(); ++col) {
		SEXP src_col = src_cols[col];
		SEXP tgt_col = tgt_cols[col + tgt_col_offset];

		if (Rf_isInteger(src_col) || Rf_isFactor(src_col))
			INTEGER(tgt_col)[tgt_row] = INTEGER(src_col)[src_row];
		else if (Rf_isReal(src_col))
			REAL(tgt_col)[tgt_row] = REAL(src_col)[src_row];
		else if (Rf_isLogical(src_col))
			LOGICAL(tgt_col)[tgt_row] = LOGICAL(src_col)[src_row];
		else if (Rf_isString(src_col))
			SET_STRING_ELT(tgt_col, tgt_row, Rf_mkChar(CHAR(STRING_ELT(src_col, src_row))));
	}
}

void DataFrameUtils::copy_data_frame_rows(const vector<SEXP> &src_cols, int src_row, int num_rows, const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col_offset)
{
	for (uint64_t col = 0; col < src_cols.size(); ++col) {
		SEXP src_col = src_cols[col];
		SEXP tgt_col = tgt_cols[col + tgt_col_offset];

		if (Rf_isInteger(src_col) || Rf_isFactor(src_col)) {
			int *src_vals = INTEGER(src_col);
			int *tgt_vals = INTEGER(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isReal(src_col)) {
			double *src_vals = REAL(src_col);
			double *tgt_vals = REAL(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isLogical(src_col)) {
			int *src_vals = LOGICAL(src_col);
			int *tgt_vals = LOGICAL(tgt_col);
			for (int i = 0; i < num_rows; ++i) 
				tgt_vals[tgt_row + i] = src_vals[src_row + i];
		} else if (Rf_isString(src_col)) {
			for (int i = 0; i < num_rows; ++i) 
				SET_STRING_ELT(tgt_col, tgt_row + i, Rf_mkChar(CHAR(STRING_ELT(src_col, src_row + i))));
		}
	}
}

void DataFrameUtils::set_data_frame_val_nan(const vector<SEXP> &tgt_cols, int tgt_row, int tgt_col)
{
	SEXP rtgt_col = tgt_cols[tgt_col];

	if (Rf_isInteger(rtgt_col) || Rf_isFactor(rtgt_col))
		INTEGER(rtgt_col)[tgt_row] = NA_INTEGER;
	else if (Rf_isReal(rtgt_col))
		REAL(rtgt_col)[tgt_row] = NA_REAL;
	else if (Rf_isLogical(rtgt_col))
		LOGICAL(rtgt_col)[tgt_row] = NA_LOGICAL;
	else if (Rf_isString(rtgt_col))
		SET_STRING_ELT(rtgt_col, tgt_row, NA_STRING);
}

