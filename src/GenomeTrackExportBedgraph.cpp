/*
 * GenomeTrackExportBedgraph.cpp
 *
 * gtrack.export_bedgraph: the data lines of a bedGraph file. They are formatted as the R code
 * formatted them, paste(chrom, sprintf("%.0f", start), sprintf("%.0f", end), value, sep = "\t"),
 * so the files are the same byte for byte; building 85M strings in R took minutes.
 */

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "rdbutils.h"

using namespace std;
using namespace rdb;

namespace {

// as.character() of a double, for a finite value: R's formatReal() for one value with
// R_print.digits = DBL_DIG, then EncodeRealDrop0() (src/main/format.c and printutils.c, the same
// in R 4.3 and 4.4). This is the code path of R built with a long double wider than double
// (capabilities("long.double")); the R caller does not use this function otherwise.

const int R_DIGITS = 15;  // DBL_DIG, as coerceToString() sets R_print.digits
const int KP_MAX = 27;
const int NB = 1000;

// double literals, as in R: 1e23 and above are not exact powers of ten
const long double tbl[] = {
	1e00, 1e01, 1e02, 1e03, 1e04, 1e05, 1e06, 1e07, 1e08, 1e09,
	1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
	1e20, 1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27
};

void scientific(double x, int *neg, int *kpower, int *nsig, bool *roundingwidens)
{
	double alpha, r;
	int kp;

	if (x == 0.0) {
		*kpower = 0;
		*nsig = 1;
		*neg = 0;
		*roundingwidens = false;
		return;
	}

	if (x < 0.0) {
		*neg = 1;
		r = -x;
	} else {
		*neg = 0;
		r = x;
	}

	kp = (int)floor(log10(r)) - R_DIGITS + 1;
	long double r_prec = r;
	if (abs(kp) <= KP_MAX) {
		if (kp > 0)
			r_prec /= tbl[kp];
		else if (kp < 0)
			r_prec *= tbl[-kp];
	} else
		r_prec /= powl(10.0, (long double)kp);
	if (r_prec < tbl[R_DIGITS - 1]) {
		r_prec *= 10.0;
		kp--;
	}
	alpha = (double)nearbyintl(r_prec);

	*nsig = R_DIGITS;
	for (int j = 1; j <= R_DIGITS; j++) {
		alpha /= 10.0;
		if (alpha == floor(alpha))
			(*nsig)--;
		else
			break;
	}
	if (*nsig == 0) {
		*nsig = 1;
		kp += 1;
	}
	*kpower = kp + R_DIGITS - 1;

	int rgt = R_DIGITS - *kpower;
	rgt = rgt < 0 ? 0 : rgt > KP_MAX ? KP_MAX : rgt;
	double fuzz = 0.5 / (double)tbl[rgt];
	*roundingwidens = *kpower > 0 && *kpower <= KP_MAX && r < tbl[*kpower] - fuzz;
}

// Appends x (finite) as as.character(x) writes it; dec is getOption("OutDec")
void append_double(string &out, double x, int scipen, const string &dec)
{
	int neg, kpower, nsig;
	bool roundingwidens;

	scientific(x, &neg, &kpower, &nsig, &roundingwidens);

	int left = kpower + 1;
	if (roundingwidens)
		left--;
	int sleft = neg + ((left <= 0) ? 1 : left);
	int rgt = nsig - left;
	int mxsl = sleft;

	if (left < 0)
		mxsl = 1 + neg;
	if (rgt < 0)
		rgt = 0;
	int wF = mxsl + rgt + (rgt != 0);

	int e = (left > 100 || left <= -99) ? 2 : 1;
	int d = nsig - 1;
	int w = neg + (d > 0) + d + 4 + e;
	if (wF <= w + scipen) {
		e = 0;
		d = rgt;
		w = wF;
	}
	w = min(w, NB - 1);

	char buf[NB];
	if (x == 0.0)
		x = 0.0;  // no "-0"
	if (e)
		snprintf(buf, NB, d ? "%#*.*e" : "%*.*e", w, d, x);
	else
		snprintf(buf, NB, "%*.*f", w, d, x);
	buf[NB - 1] = '\0';

	// drop the trailing zeros of the fraction, as EncodeRealDrop0() does
	for (char *p = buf; *p; p++) {
		if (*p == '.') {
			char *replace = p++;
			while ('0' <= *p && *p <= '9') {
				if (*(p++) != '0')
					replace = p;
			}
			if (replace != p) {
				while ((*(replace++) = *(p++)))
					;
			}
			break;
		}
	}

	if (dec == ".")
		out += buf;
	else {
		for (const char *p = buf; *p; ++p) {
			if (*p == '.')
				out += dec;
			else
				out += *p;
		}
	}
}

// as.character() of a double, NaN aside (the caller drops NaN rows)
void append_value(string &out, double x, int scipen, const string &dec)
{
	if (R_IsNA(x))
		out += "NA";
	else if (x == R_PosInf)
		out += "Inf";
	else if (x == R_NegInf)
		out += "-Inf";
	else
		append_double(out, x, scipen, dec);
}

// sprintf("%.0f", x) of R
void append_coord(string &out, double x)
{
	char buf[NB];

	if (x == floor(x) && fabs(x) < 9e15 && !(x == 0 && signbit(x))) {
		int64_t v = (int64_t)x;
		char *p = buf + sizeof(buf);
		bool negative = v < 0;
		uint64_t u = negative ? -(uint64_t)v : (uint64_t)v;

		do {
			*--p = '0' + u % 10;
			u /= 10;
		} while (u);
		if (negative)
			*--p = '-';
		out.append(p, buf + sizeof(buf) - p);
	} else if (R_IsNA(x))
		out += "NA";
	else {
		snprintf(buf, sizeof(buf), "%.0f", x);
		out += buf;
	}
}

}

extern "C" {

// The bedGraph lines of the rows order[from - 1], order[from], ... (1-based row numbers; NULL
// order: the rows in order) whose value is not NaN (NA is written as "NA", as paste() writes it),
// until the text exceeds max_bytes. chrom is a factor or a character vector, start and end
// numeric. Returns list(text, next): next is the 1-based position in order to continue from,
// length + 1 when done.
SEXP gtrack_bedgraph_lines(SEXP _chrom, SEXP _start, SEXP _end, SEXP _value, SEXP _order, SEXP _from, SEXP _max_bytes,
                           SEXP _scipen, SEXP _dec, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		R_xlen_t n = Rf_xlength(_value);
		bool chrom_is_factor = Rf_isFactor(_chrom);

		if (!Rf_isReal(_value))
			verror("value must be a double vector");
		if (Rf_xlength(_chrom) != n || Rf_xlength(_start) != n || Rf_xlength(_end) != n)
			verror("chrom, start, end and value must have the same length");
		if (!chrom_is_factor && !Rf_isString(_chrom))
			verror("chrom must be a factor or a character vector");
		if ((!Rf_isReal(_start) && !Rf_isInteger(_start)) || (!Rf_isReal(_end) && !Rf_isInteger(_end)))
			verror("start and end must be numeric");
		if (!Rf_isNull(_order) && !Rf_isInteger(_order))
			verror("order must be an integer vector or NULL");
		if (!Rf_isReal(_from) || Rf_length(_from) != 1 || !Rf_isReal(_max_bytes) || Rf_length(_max_bytes) != 1)
			verror("from and max_bytes must be numbers");
		if (!Rf_isInteger(_scipen) || Rf_length(_scipen) != 1 || !Rf_isString(_dec) || Rf_length(_dec) != 1)
			verror("scipen must be an integer and dec a string");

		R_xlen_t num_rows = Rf_isNull(_order) ? n : Rf_xlength(_order);
		const int *order = Rf_isNull(_order) ? NULL : INTEGER(_order);
		R_xlen_t from = (R_xlen_t)REAL(_from)[0];
		uint64_t max_bytes = (uint64_t)REAL(_max_bytes)[0];
		int scipen = INTEGER(_scipen)[0];
		string dec = CHAR(STRING_ELT(_dec, 0));
		const double *value = REAL(_value);

		if (from < 1)
			verror("from must be >= 1");
		if (scipen == NA_INTEGER)
			scipen = 0;

		vector<string> levels;
		if (chrom_is_factor) {
			SEXP lv = Rf_getAttrib(_chrom, R_LevelsSymbol);
			for (R_xlen_t i = 0; i < Rf_xlength(lv); ++i)
				levels.push_back(CHAR(STRING_ELT(lv, i)));
		}

		auto coord = [](SEXP v, R_xlen_t i) -> double {
			if (Rf_isReal(v))
				return REAL(v)[i];
			int x = INTEGER(v)[i];
			return x == NA_INTEGER ? NA_REAL : (double)x;
		};

		string out;
		out.reserve(max_bytes + 1024);

		R_xlen_t i;
		for (i = from - 1; i < num_rows && out.size() < max_bytes; ++i) {
			R_xlen_t row = order ? (R_xlen_t)order[i] - 1 : i;

			if (row < 0 || row >= n)
				verror("order holds a row out of range");

			double v = value[row];
			if (std::isnan(v) && !R_IsNA(v))
				continue;

			if (chrom_is_factor) {
				int code = INTEGER(_chrom)[row];
				if (code == NA_INTEGER || code < 1 || code > (int)levels.size())
					out += "NA";
				else
					out += levels[code - 1];
			} else {
				SEXP s = STRING_ELT(_chrom, row);
				out += s == NA_STRING ? "NA" : CHAR(s);
			}
			out += '\t';
			append_coord(out, coord(_start, row));
			out += '\t';
			append_coord(out, coord(_end, row));
			out += '\t';
			append_value(out, v, scipen, dec);
			out += '\n';

			if (!(i % 1000000))
				check_interrupt();
		}

		if (out.size() > INT_MAX)
			verror("A bedGraph chunk of %llu bytes is too long for an R string", (unsigned long long)out.size());

		SEXP answer = rprotect_ptr(RSaneAllocVector(VECSXP, 2));
		SEXP text = rprotect_ptr(RSaneAllocVector(STRSXP, 1));
		SET_STRING_ELT(text, 0, RSaneMkCharLenCE(out.data(), (int)out.size(), CE_NATIVE));
		SEXP next = rprotect_ptr(RSaneAllocVector(REALSXP, 1));
		REAL(next)[0] = (double)i + 1;
		SET_VECTOR_ELT(answer, 0, text);
		SET_VECTOR_ELT(answer, 1, next);
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

}
