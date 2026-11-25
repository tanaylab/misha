// gseq.onehot / gextract.windows: one-hot arrays of genomic windows, read straight from the
// sequence store (no character strings in between), and of character strings. The encoding
// itself is in SeqOneHot.h, shared verbatim with pymisha.

#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "port.h"

#include "GenomeSeqFetch.h"
#include "SeqOneHot.h"
#include "rdbinterval.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// dtype codes, as R/sequence-onehot.R passes them
enum OneHotType { ONEHOT_RAW = 0, ONEHOT_INTEGER = 1, ONEHOT_DOUBLE = 2 };

static SEXPTYPE onehot_sexptype(SEXP _dtype)
{
	if (!Rf_isInteger(_dtype) || Rf_length(_dtype) != 1)
		verror("dtype must be an integer code");
	switch (INTEGER(_dtype)[0]) {
	case ONEHOT_RAW: return RAWSXP;
	case ONEHOT_INTEGER: return INTSXP;
	case ONEHOT_DOUBLE: return REALSXP;
	default: verror("Unknown dtype code %d", INTEGER(_dtype)[0]);
	}
	return RAWSXP;
}

// An array of dim c(num, 4, width) of the given type, protected.
static SEXP alloc_onehot(SEXPTYPE type, int64_t num, int64_t width)
{
	R_xlen_t len = (R_xlen_t)num * 4 * (R_xlen_t)width;
	SEXP answer;
	SEXP dim;

	// not initialized: the encoders write every element, zeros included
	rprotect(answer = RSaneAllocVector(type, len));

	rprotect(dim = RSaneAllocVector(INTSXP, 3));
	INTEGER(dim)[0] = (int)num;
	INTEGER(dim)[1] = 4;
	INTEGER(dim)[2] = (int)width;
	Rf_setAttrib(answer, R_DimSymbol, dim);
	runprotect(dim);
	return answer;
}

// _chroms: character (normalized chromosome names), _starts: numeric window starts, _width: the
// window width, _revcomp: logical (one per window), _pad: TRUE lets a window run past the ends of
// its chromosome (the outside positions are all zero), _dtype: see OneHotType.
extern "C" {

SEXP C_gseq_onehot(SEXP _chroms, SEXP _starts, SEXP _width, SEXP _revcomp, SEXP _pad, SEXP _dtype, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;
		IntervUtils iu(_envir);

		if (!Rf_isString(_chroms))
			verror("chroms must be a character vector");
		if (!Rf_isReal(_starts) && !Rf_isInteger(_starts))
			verror("starts must be numeric");
		if (!Rf_isLogical(_revcomp) || !Rf_isLogical(_pad) || Rf_length(_pad) != 1)
			verror("revcomp and pad must be logical");
		if ((!Rf_isReal(_width) && !Rf_isInteger(_width)) || Rf_length(_width) != 1)
			verror("width must be a number");

		int64_t num = Rf_xlength(_chroms);
		if (Rf_xlength(_starts) != num || Rf_xlength(_revcomp) != num)
			verror("chroms, starts and revcomp must have the same length");

		int64_t width = (int64_t)(Rf_isReal(_width) ? REAL(_width)[0] : INTEGER(_width)[0]);
		if (width <= 0)
			verror("width must be positive");
		if (num > INT_MAX || width > INT_MAX)
			verror("Too many windows or too wide windows for an R array");

		bool pad = LOGICAL(_pad)[0] == TRUE;
		SEXPTYPE type = onehot_sexptype(_dtype);
		const GenomeChromKey &chromkey = iu.get_chromkey();

		vector<int> chromids(num);
		vector<int64_t> starts(num);
		vector<unsigned char> revcomp(num);

		for (int64_t i = 0; i < num; ++i) {
			chromids[i] = chromkey.chrom2id(CHAR(STRING_ELT(_chroms, i)));
			double start = Rf_isReal(_starts) ? REAL(_starts)[i] : (INTEGER(_starts)[i] == NA_INTEGER ? NA_REAL : INTEGER(_starts)[i]);
			if (!R_FINITE(start) || start != (double)(int64_t)start)
				verror("Window %lld: start must be a finite whole number", (long long)(i + 1));
			starts[i] = (int64_t)start;
			if (!pad && (starts[i] < 0 || starts[i] + width > (int64_t)chromkey.get_chrom_size(chromids[i])))
				verror("Window (%s, %lld, %lld) runs past the end of its chromosome (use pad = TRUE to fill with zeros)",
					   chromkey.id2chrom(chromids[i]).c_str(), (long long)starts[i], (long long)(starts[i] + width));
			revcomp[i] = LOGICAL(_revcomp)[i] == TRUE;
		}

		iu.verify_max_data_size((uint64_t)num * width, "Result one-hot array");

		SEXP answer = alloc_onehot(type, num, width);
		seqonehot::Strides s = { 1, num, 4 * num };

		GenomeSeqFetch seqfetch;
		seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

		auto interrupt = []() { check_interrupt(); };
		if (type == RAWSXP)
			seqonehot::encode_windows(seqfetch, chromkey, chromids.data(), starts.data(), revcomp.data(), num, width, 0, RAW(answer), s, interrupt);
		else if (type == INTSXP)
			seqonehot::encode_windows(seqfetch, chromkey, chromids.data(), starts.data(), revcomp.data(), num, width, 0, INTEGER(answer), s, interrupt);
		else
			seqonehot::encode_windows(seqfetch, chromkey, chromids.data(), starts.data(), revcomp.data(), num, width, 0, REAL(answer), s, interrupt);

		// still protected: ~RdbInitializer unprotects it on the way out
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

// One-hot of character strings of equal length (NA is not allowed). _revcomp: logical, recycled.
SEXP C_gseq_onehot_strings(SEXP _seqs, SEXP _revcomp, SEXP _dtype)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_seqs))
			verror("seqs must be a character vector");
		if (!Rf_isLogical(_revcomp) || Rf_length(_revcomp) < 1)
			verror("revcomp must be logical");

		int64_t num = Rf_xlength(_seqs);
		SEXPTYPE type = onehot_sexptype(_dtype);
		int64_t width = 0;

		for (int64_t i = 0; i < num; ++i) {
			SEXP s = STRING_ELT(_seqs, i);
			if (s == NA_STRING)
				verror("Sequence %lld is NA", (long long)(i + 1));
			int64_t len = Rf_xlength(s);
			if (!i)
				width = len;
			else if (len != width)
				verror("All sequences must have the same length (sequence 1 has %lld bases, sequence %lld has %lld)",
					   (long long)width, (long long)(i + 1), (long long)len);
		}
		if (num > INT_MAX || width > INT_MAX)
			verror("Too many or too long sequences for an R array");

		SEXP answer = alloc_onehot(type, num, width);
		seqonehot::Strides st = { 1, num, 4 * num };
		int64_t nrc = Rf_xlength(_revcomp);

		for (int64_t i = 0; i < num; ++i) {
			const char *seq = CHAR(STRING_ELT(_seqs, i));
			bool rc = LOGICAL(_revcomp)[i % nrc] == TRUE;
			if (type == RAWSXP)
				seqonehot::encode(seq, width, rc, i, RAW(answer), st);
			else if (type == INTSXP)
				seqonehot::encode(seq, width, rc, i, INTEGER(answer), st);
			else
				seqonehot::encode(seq, width, rc, i, REAL(answer), st);
			if (!(i % 4096))
				check_interrupt();
		}

		// still protected: ~RdbInitializer unprotects it on the way out
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}
	return R_NilValue;
}

}
