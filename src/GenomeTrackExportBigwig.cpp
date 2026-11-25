/*
 * GenomeTrackExportBigwig.cpp
 *
 * gtrack.export_bigwig: the R entry points of BigWigWriter (BigWigWriter.h). R creates the file,
 * appends the records of one batch of chromosomes at a time, as gextract returns them, and finishes
 * it; the writer's state is a raw vector that R hands from one call to the next.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <new>
#include <string>
#include <thread>
#include <vector>

#include "BigWigWriter.h"
#include "rdbinterval.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

namespace {

string file_arg(SEXP _file)
{
	if (!Rf_isString(_file) || Rf_length(_file) != 1 || STRING_ELT(_file, 0) == NA_STRING)
		verror("file must be a string");
	return CHAR(STRING_ELT(_file, 0));
}

void state_arg(SEXP _state, vector<char> &state)
{
	if (TYPEOF(_state) != RAWSXP)
		verror("Invalid bigWig writer state");
	state.assign((const char *)RAW(_state), (const char *)RAW(_state) + Rf_xlength(_state));
}

SEXP state_value(const vector<char> &state)
{
	SEXP answer = rprotect_ptr(RSaneAllocVector(RAWSXP, state.size()));

	if (!state.empty())
		memcpy(RAW(answer), state.data(), state.size());
	return answer;
}

// A numeric vector of whole numbers as int64
void int64_arg(SEXP _v, const char *name, vector<int64_t> &res)
{
	R_xlen_t n = Rf_xlength(_v);

	res.resize(n);
	if (Rf_isInteger(_v)) {
		const int *p = INTEGER(_v);
		for (R_xlen_t i = 0; i < n; ++i) {
			if (p[i] == NA_INTEGER)
				verror("%s holds NA", name);
			res[i] = p[i];
		}
	} else if (Rf_isReal(_v)) {
		const double *p = REAL(_v);
		for (R_xlen_t i = 0; i < n; ++i) {
			if (!std::isfinite(p[i]) || p[i] != floor(p[i]) || fabs(p[i]) > 9.0e15)
				verror("%s must hold whole numbers, not %g", name, p[i]);
			res[i] = (int64_t)p[i];
		}
	} else
		verror("%s must be numeric", name);
}

// The threads that compress blocks: at most gmax.processes and 8
unsigned num_threads(SEXP _envir)
{
	IntervUtils iu(_envir);
	unsigned hw = std::thread::hardware_concurrency();

	return (unsigned)min<uint64_t>({ hw ? hw : 1, max<uint64_t>(iu.get_max_processes(), 1), 8 });
}

}

extern "C" {

// Creates the file: the chromosomes are names (sorted byte by byte) and sizes. Returns the state.
SEXP gtrack_bigwig_create(SEXP _file, SEXP _names, SEXP _sizes, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		string file = file_arg(_file);

		if (!Rf_isString(_names))
			verror("names must be a character vector");

		vector<int64_t> sizes;

		int64_arg(_sizes, "sizes", sizes);
		if ((R_xlen_t)sizes.size() != Rf_xlength(_names))
			verror("names and sizes must have the same length");

		vector<BigWigWriter::Chrom> chroms(sizes.size());

		for (uint64_t i = 0; i < sizes.size(); ++i) {
			if (STRING_ELT(_names, i) == NA_STRING)
				verror("names holds NA");
			chroms[i].name = CHAR(STRING_ELT(_names, i));
			chroms[i].size = sizes[i] < 0 ? 0 : (uint64_t)sizes[i];
		}

		vector<char> state;

		BigWigWriter::create(file, chroms, state);
		return state_value(state);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

// Appends records: chromids (0-based, into the names of gtrack_bigwig_create), starts, ends and
// values, in the order of the file. Returns the new state.
SEXP gtrack_bigwig_append(SEXP _file, SEXP _state, SEXP _chromids, SEXP _starts, SEXP _ends, SEXP _values, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		string file = file_arg(_file);
		vector<char> state;
		vector<int64_t> chromids, starts, ends;

		state_arg(_state, state);
		int64_arg(_chromids, "chromids", chromids);
		int64_arg(_starts, "starts", starts);
		int64_arg(_ends, "ends", ends);
		if (!Rf_isReal(_values))
			verror("values must be a double vector");

		uint64_t n = chromids.size();

		if (starts.size() != n || ends.size() != n || (uint64_t)Rf_xlength(_values) != n)
			verror("chromids, starts, ends and values must have the same length");

		BigWigWriter::append(file, state, n, chromids.data(), starts.data(), ends.data(), REAL(_values), num_threads(_envir));
		return state_value(state);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

// Writes the index, the zoom levels, the summary and the header. Returns the number of records.
SEXP gtrack_bigwig_finish(SEXP _file, SEXP _state, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		string file = file_arg(_file);
		vector<char> state;

		state_arg(_state, state);

		uint64_t num_items = BigWigWriter::finish(file, state, num_threads(_envir));
		SEXP answer = rprotect_ptr(RSaneAllocVector(REALSXP, 1));

		REAL(answer)[0] = (double)num_items;
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
