/*
 * GenomeTrackExportBedgraph.cpp
 *
 * gtrack.export_bedgraph: the data lines of a bedGraph file. They are formatted as the R code
 * formatted them, paste(chrom, sprintf("%.0f", start), sprintf("%.0f", end), value, sep = "\t"),
 * so the files are the same byte for byte; building 85M strings in R took minutes. Rows are formatted
 * in threads.
 */

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include "rdbinterval.h"
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

// The low 32 bits of NA_real_, which R_IsNA() compares on a NaN. Set by the calling thread before
// any worker starts: the workers use no R API.
uint32_t s_na_low_word = 1954;

// R_IsNA(), on the bits
inline bool is_na_real(double x)
{
	if (!std::isnan(x))
		return false;
	uint64_t bits;
	memcpy(&bits, &x, sizeof(bits));
	return (uint32_t)(bits & 0xffffffff) == s_na_low_word;
}

// as.character() of a double, NaN aside (the caller drops NaN rows)
void append_value(string &out, double x, int scipen, const string &dec)
{
	if (is_na_real(x))
		out += "NA";
	else if (std::isinf(x))
		out += x > 0 ? "Inf" : "-Inf";
	else
		append_double(out, x, scipen, dec);
}

// sprintf("%.0f", x) of R; "NA" if na
void append_coord(string &out, double x, bool na = false)
{
	if (na) {
		out += "NA";
		return;
	}

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
	} else if (is_na_real(x))
		out += "NA";
	else {
		snprintf(buf, sizeof(buf), "%.0f", x);
		out += buf;
	}
}

// The lines of rows [begin, end) of order (or of the data, without order) whose value is not
// NaN, appended to out. Pure C++ on the vectors' data: safe in a thread. Sets bad_row on a row
// number out of range.
struct BedgraphRows {
	R_xlen_t                n;
	const int              *order;
	const int              *chrom_codes;   // a factor's codes, or NULL
	SEXP                    chrom_strings; // a character vector, read only by the calling thread
	const vector<string>   *levels;
	const double           *start_real;
	const int              *start_int;
	const double           *end_real;
	const int              *end_int;
	const double           *value;
	int                     scipen;
	const string           *dec;
	int                     na_integer;    // NA_INTEGER, read by the calling thread

	void append(string &out, R_xlen_t begin, R_xlen_t end, bool &bad_row) const
	{
		for (R_xlen_t i = begin; i < end; ++i) {
			R_xlen_t row = order ? (R_xlen_t)order[i] - 1 : i;

			if (row < 0 || row >= n) {
				bad_row = true;
				return;
			}

			double v = value[row];
			if (std::isnan(v) && !is_na_real(v))
				continue;

			if (chrom_codes) {
				int code = chrom_codes[row];
				if (code == na_integer || code < 1 || code > (int)levels->size())
					out += "NA";
				else
					out += (*levels)[code - 1];
			} else {
				SEXP s = STRING_ELT(chrom_strings, row);
				out += s == NA_STRING ? "NA" : CHAR(s);
			}
			out += '\t';
			if (start_real)
				append_coord(out, start_real[row]);
			else
				append_coord(out, start_int[row], start_int[row] == na_integer);
			out += '\t';
			if (end_real)
				append_coord(out, end_real[row]);
			else
				append_coord(out, end_int[row], end_int[row] == na_integer);
			out += '\t';
			append_value(out, v, scipen, *dec);
			out += '\n';
		}
	}
};

}

extern "C" {

// The bedGraph lines of the rows order[from - 1], ..., order[from + max_rows - 2] (1-based row
// numbers; NULL order: the rows in order) whose value is not NaN (NA is written as "NA", as
// paste() writes it). chrom is a factor or a character vector, start and end numeric. Returns
// list(text, next): next is the 1-based position in order to continue from, length + 1 when
// done. The rows of a factor chrom are formatted in threads (at most gmax.processes and 8).
SEXP gtrack_bedgraph_lines(SEXP _chrom, SEXP _start, SEXP _end, SEXP _value, SEXP _order, SEXP _from, SEXP _max_rows,
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
		if (!Rf_isReal(_from) || Rf_length(_from) != 1 || !Rf_isReal(_max_rows) || Rf_length(_max_rows) != 1)
			verror("from and max_rows must be numbers");
		if (!Rf_isInteger(_scipen) || Rf_length(_scipen) != 1 || !Rf_isString(_dec) || Rf_length(_dec) != 1)
			verror("scipen must be an integer and dec a string");

		R_xlen_t num_rows = Rf_isNull(_order) ? n : Rf_xlength(_order);
		R_xlen_t from = (R_xlen_t)REAL(_from)[0];
		double max_rows = REAL(_max_rows)[0];
		string dec = CHAR(STRING_ELT(_dec, 0));

		if (from < 1 || !(max_rows >= 1))
			verror("from and max_rows must be >= 1");

		R_xlen_t begin = from - 1;
		R_xlen_t end = begin + (R_xlen_t)min(max_rows, (double)(num_rows - begin));
		if (end < begin)
			end = begin;

		vector<string> levels;
		if (chrom_is_factor) {
			SEXP lv = Rf_getAttrib(_chrom, R_LevelsSymbol);
			for (R_xlen_t i = 0; i < Rf_xlength(lv); ++i)
				levels.push_back(CHAR(STRING_ELT(lv, i)));
		}

		BedgraphRows rows;
		rows.n = n;
		rows.order = Rf_isNull(_order) ? NULL : INTEGER(_order);
		rows.chrom_codes = chrom_is_factor ? INTEGER(_chrom) : NULL;
		rows.chrom_strings = _chrom;
		rows.levels = &levels;
		rows.start_real = Rf_isReal(_start) ? REAL(_start) : NULL;
		rows.start_int = Rf_isReal(_start) ? NULL : INTEGER(_start);
		rows.end_real = Rf_isReal(_end) ? REAL(_end) : NULL;
		rows.end_int = Rf_isReal(_end) ? NULL : INTEGER(_end);
		rows.value = REAL(_value);
		rows.scipen = INTEGER(_scipen)[0] == NA_INTEGER ? 0 : INTEGER(_scipen)[0];
		rows.dec = &dec;
		rows.na_integer = NA_INTEGER;
		{
			double na = NA_REAL;
			uint64_t bits;
			memcpy(&bits, &na, sizeof(bits));
			s_na_low_word = (uint32_t)(bits & 0xffffffff);
		}

		// Threads only when the rows hold no R strings: STRING_ELT is not for threads
		unsigned num_threads = 1;
		if (chrom_is_factor && end - begin >= 200000) {
			IntervUtils iu(_envir);
			unsigned hw = std::thread::hardware_concurrency();
			num_threads = min<unsigned>({ hw ? hw : 1, (unsigned)max<uint64_t>(iu.get_max_processes(), 1), 8 });
		}

		vector<string> parts(num_threads);
		vector<char> bad(num_threads, 0);
		auto work = [&](unsigned t) {
			R_xlen_t b = begin + (end - begin) * t / num_threads;
			R_xlen_t e = begin + (end - begin) * (t + 1) / num_threads;
			bool bad_row = false;
			parts[t].reserve((e - b) * 40);
			rows.append(parts[t], b, e, bad_row);
			bad[t] = bad_row;
		};

		if (num_threads == 1)
			work(0);
		else {
			vector<thread> threads;
			for (unsigned t = 1; t < num_threads; ++t)
				threads.emplace_back(work, t);
			work(0);
			for (thread &th : threads)
				th.join();
		}

		uint64_t size = 0;
		for (unsigned t = 0; t < num_threads; ++t) {
			if (bad[t])
				verror("order holds a row out of range");
			size += parts[t].size();
		}
		if (size > INT_MAX)
			verror("A bedGraph chunk of %llu bytes is too long for an R string", (unsigned long long)size);

		string out;
		out.reserve(size);
		for (string &part : parts) {
			out += part;
			string().swap(part);
		}

		SEXP answer = rprotect_ptr(RSaneAllocVector(VECSXP, 2));
		SEXP text = rprotect_ptr(RSaneAllocVector(STRSXP, 1));
		SET_STRING_ELT(text, 0, RSaneMkCharLenCE(out.data(), (int)out.size(), CE_NATIVE));
		SEXP next = rprotect_ptr(RSaneAllocVector(REALSXP, 1));
		REAL(next)[0] = (double)end + 1;
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
