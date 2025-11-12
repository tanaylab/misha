#include <cstdint>
#include "GIntervalsBigSet1D.h"
#include "GenomeSeqFetch.h"
#include "GenomeIndex.h"
#include "rdbinterval.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gseqread(SEXP _intervals, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		IntervUtils iu(_envir);
		GIntervalsFetcher1D *intervals = NULL;
		iu.convert_rintervs(_intervals, &intervals, NULL);
		unique_ptr<GIntervalsFetcher1D> intervals_guard(intervals);
		intervals->sort();

		if (!intervals->size())
			return R_NilValue;

        vector<char> buf;
        // Reserve based on a heuristic: start with 256KB and grow as needed
        buf.reserve(262144);
		SEXP answer;
		rprotect(answer = RSaneAllocVector(STRSXP, intervals->size()));

		uint64_t seqlen = 0;
		GenomeSeqFetch seqfetch;
		seqfetch.set_seqdir(string(rdb::get_groot(_envir)) + "/seq");

		for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
			seqfetch.read_interval(intervals->cur_interval(), iu.get_chromkey(), buf);
            seqlen += buf.size();
			iu.verify_max_data_size(seqlen, "Result sequence");
            // Avoid strlen by constructing string with known length
            SET_STRING_ELT(answer,
                           iu.get_orig_interv_idx(intervals->cur_interval()),
                           Rf_mkCharLenCE(&*buf.begin(), (int)buf.size(), CE_NATIVE));
			check_interrupt();
		}
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

// Validate genome index file (called during gdb.init)
// Throws error if index is corrupt or has checksum mismatch
SEXP gseq_validate_index(SEXP _seqdir, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_seqdir) || Rf_length(_seqdir) != 1)
			verror("seqdir argument must be a string");

		const char *seqdir = CHAR(STRING_ELT(_seqdir, 0));
		string idx_path = string(seqdir) + "/genome.idx";

		// Try to load index - will throw TGLException if corrupt
		GenomeIndex index;
		index.load(idx_path);

		return R_NilValue;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
