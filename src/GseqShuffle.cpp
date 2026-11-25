// C_gseq_shuffle: gseq.shuffle()'s engine. The shuffle itself is SeqShuffle.cpp,
// shared with pymisha, so the two give the same sequences for the same seed.

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "SeqShuffle.h"
#include "rdbutils.h"

using namespace rdb;
using namespace std;

extern "C" {

// r_seqs: character vector (NA allowed); r_k, r_n: integers; r_seed: a whole number
// with |seed| <= 2^53. Returns a character vector of length(seqs) * n, the n
// shuffles of the first sequence first.
SEXP C_gseq_shuffle(SEXP r_seqs, SEXP r_k, SEXP r_n, SEXP r_seed)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(r_seqs))
			verror("seqs must be a character vector");
		if (!Rf_isInteger(r_k) || Rf_length(r_k) != 1 || INTEGER(r_k)[0] == NA_INTEGER)
			verror("k must be a single integer");
		if (!Rf_isInteger(r_n) || Rf_length(r_n) != 1 || INTEGER(r_n)[0] == NA_INTEGER)
			verror("n must be a single integer");
		if (!Rf_isReal(r_seed) || Rf_length(r_seed) != 1)
			verror("seed must be a single number");

		const int k = INTEGER(r_k)[0];
		const int n = INTEGER(r_n)[0];
		const double seed_d = REAL(r_seed)[0];
		if (k < 1 || k > seqshuffle::MAX_K)
			verror("k must be between 1 and %d", seqshuffle::MAX_K);
		if (n < 1)
			verror("n must be a positive integer");
		if (!std::isfinite(seed_d) || seed_d != std::floor(seed_d) || std::fabs(seed_d) > 9007199254740992.0)
			verror("seed must be a whole number between -2^53 and 2^53");
		const int64_t seed = (int64_t)seed_d;

		const R_xlen_t num_seqs = Rf_xlength(r_seqs);
		if ((double)num_seqs * n > (double)R_XLEN_T_MAX)
			verror("length(seqs) * n is too large");

		SEXP res = rprotect_ptr(RSaneAllocVector(STRSXP, num_seqs * n));
		string buf;

		for (R_xlen_t i = 0; i < num_seqs; ++i) {
			check_interrupt();
			SEXP el = STRING_ELT(r_seqs, i);
			if (el == NA_STRING) {
				for (int j = 0; j < n; ++j)
					SET_STRING_ELT(res, i * n + j, NA_STRING);
				continue;
			}
			const char *seq = CHAR(el);
			const R_xlen_t len = Rf_xlength(el);
			if (len >= (R_xlen_t)INT32_MAX)
				verror("Sequence %lld is too long to shuffle (%lld bases)", (long long)(i + 1), (long long)len);
			const cetype_t enc = Rf_getCharCE(el);

			seqshuffle::Shuffler shuffler(seq, (size_t)len, k);
			seqshuffle::Rng rng(seed, (uint64_t)i);
			buf.resize((size_t)len);
			for (int j = 0; j < n; ++j) {
				if (j)
					check_interrupt();
				shuffler.shuffle(rng, len ? &buf[0] : nullptr);
				SET_STRING_ELT(res, i * n + j, RSaneMkCharLenCE(buf.data(), (int)len, enc));
			}
		}

		runprotect(1);
		return res;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const std::bad_alloc &) {
		rerror("Out of memory in gseq.shuffle");
	} catch (const std::exception &e) {
		rerror("Error in C_gseq_shuffle: %s", e.what());
	}
	rerror("Unknown error in C_gseq_shuffle");
	return R_NilValue;
}

} // extern "C"
