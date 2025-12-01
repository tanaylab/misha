#include "ChainIntervalConverter.h"
#include "rdbinterval.h" // For ChainInterval, ChainIntervals, IntervUtils
#include "rdbutils.h" // For verror, rprotect, runprotect, RSaneAllocVector, TGLError
#include "GIntervals.h"
#include <cstring> // For strcmp
#include <unordered_map>
#include <cmath> // For isnan

using namespace std;
using namespace rdb;

ChainIntervalConverter::ChainIntervalConverter(rdb::IntervUtils &iu) :
	m_iu(iu)
{
}

void ChainIntervalConverter::convert_rchain_intervs(SEXP rchain, rdb::ChainIntervals &chain_intervs, std::vector<std::string> &src_id2chrom)
{
	if (!Rf_isVector(rchain) || Rf_length(rchain) != ChainInterval::NUM_COLS)
		TGLError("Invalid format of chain argument");

	SEXP colnames = Rf_getAttrib(rchain, R_NamesSymbol);

	if (!Rf_isString(colnames) || Rf_length(colnames) != ChainInterval::NUM_COLS)
		verror("Invalid format of chain argument");

	for (unsigned i = 0; i < ChainInterval::NUM_COLS; i++) {
		if (strcmp(CHAR(STRING_ELT(colnames, i)), ChainInterval::COL_NAMES[i]))
			verror("Invalid format of chain argument");
	}

	// convert the first 3 columns of the data frame
	GIntervals intervs;
	m_iu.convert_rintervs(rchain, &intervs, NULL);

	// convert the rest of the columns
	SEXP src_chroms = VECTOR_ELT(rchain, ChainInterval::CHROM_SRC);
	SEXP src_chrom_levels = Rf_getAttrib(src_chroms, R_LevelsSymbol);
	SEXP src_starts = VECTOR_ELT(rchain, ChainInterval::START_SRC);
	// Note: src_ends (END_SRC) is not used as it's calculated from start_src + (end - start) in constructor
	SEXP src_strands = VECTOR_ELT(rchain, ChainInterval::STRAND_SRC);
	SEXP tgt_strands = VECTOR_ELT(rchain, ChainInterval::STRAND);
	SEXP chain_ids = VECTOR_ELT(rchain, ChainInterval::CHAIN_ID);
	SEXP scores = VECTOR_ELT(rchain, ChainInterval::SCORE);

	for (unsigned i = 0; i < ChainInterval::NUM_COLS; i++) {
		if (i != 0 && Rf_length(VECTOR_ELT(rchain, i)) != Rf_length(VECTOR_ELT(rchain, i - 1)))
			verror("Number of rows in column %s differs than the number of rows in column %s", ChainInterval::COL_NAMES[i - 1], ChainInterval::COL_NAMES[i]);
	}

	if (!Rf_isReal(src_starts) && !Rf_isInteger(src_starts))
		verror("Invalid format of intervals argument");

	unordered_map<string, int> src_chrom2id;

	for (unsigned i = 0; i < intervs.size(); i++) {
		if ((Rf_isFactor(src_chroms) && INTEGER(src_chroms)[i] < 0) || (Rf_isReal(src_starts) && std::isnan(REAL(src_starts)[i])))
			verror("Invalid format of interval at index %d", i + 1);

		const char *src_chrom = Rf_isString(src_chroms) ? CHAR(STRING_ELT(src_chroms, i)) : CHAR(STRING_ELT(src_chrom_levels, INTEGER(src_chroms)[i] - 1));
		unordered_map<string, int>::const_iterator isrc_chrom2id = src_chrom2id.find(src_chrom);
		int src_chromid;

		if (isrc_chrom2id == src_chrom2id.end()) {
			src_chromid = src_id2chrom.size();
			src_id2chrom.push_back(src_chrom);
			src_chrom2id[src_chrom] = src_chromid;
		} else
			src_chromid = isrc_chrom2id->second;

		int64_t src_start = (int64_t)(Rf_isReal(src_starts) ? REAL(src_starts)[i] : INTEGER(src_starts)[i]);

		// Read strands and convert +1/-1 to 0/1 for internal storage
		int tgt_strand_val = Rf_isReal(tgt_strands) ? (int)REAL(tgt_strands)[i] : INTEGER(tgt_strands)[i];
		int src_strand_val = Rf_isReal(src_strands) ? (int)REAL(src_strands)[i] : INTEGER(src_strands)[i];
		int tgt_strand = (tgt_strand_val == 1) ? 0 : 1;
		int src_strand = (src_strand_val == 1) ? 0 : 1;

		ChainInterval interval(intervs[i].chromid, intervs[i].start, intervs[i].end, tgt_strand, src_chromid, src_start, src_strand);

		// Read chain_id and score from R object if available
		if (Rf_isReal(chain_ids))
			interval.chain_id = (int64_t)REAL(chain_ids)[i];
		else if (Rf_isInteger(chain_ids))
			interval.chain_id = INTEGER(chain_ids)[i];
		else
			interval.chain_id = i;  // Fallback: assign monotonically increasing ID for stable sorting

		if (Rf_isReal(scores))
			interval.score = REAL(scores)[i];
		else if (Rf_isInteger(scores))
			interval.score = (double)INTEGER(scores)[i];
		else
			interval.score = 0.0;  // Fallback: default score

		// Validate that chain intervals are non-zero length (strict validation)
		if (interval.start >= interval.end)
			verror("Chain file contains zero-length or invalid interval at row %d (start=%lld, end=%lld)",
			       i + 1, (long long)interval.start, (long long)interval.end);

		// Don't check chromosome boundaries - chain files may have coordinates that exceed database sizes
		// (matching Kent's liftOver behavior which is lenient about target coordinates)
		interval.verify(m_iu.get_chromkey(), src_id2chrom, false);
		chain_intervs.push_back(interval);
	}
}

SEXP ChainIntervalConverter::convert_chain_intervs(const rdb::ChainIntervals &chain_intervs, std::vector<std::string> &src_id2chrom)
{
	GIntervals tmp_intervals;
	tmp_intervals.reserve(chain_intervs.size());
	for (ChainIntervals::const_iterator iinterval = chain_intervs.begin(); iinterval != chain_intervs.end(); ++iinterval)
		tmp_intervals.push_back((GInterval)*iinterval);

    SEXP answer = m_iu.convert_intervs(&tmp_intervals, ChainInterval::NUM_COLS);
	SEXP src_chroms, src_chroms_idx, src_starts, src_ends, src_strands, tgt_strands, chain_ids, scores;
    SEXP col_names = Rf_getAttrib(answer, R_NamesSymbol);
    rprotect(col_names);
	unsigned num_src_chroms = src_id2chrom.size();

    rprotect(src_chroms_idx = RSaneAllocVector(INTSXP, chain_intervs.size()));
    rprotect(src_starts = RSaneAllocVector(REALSXP, chain_intervs.size()));
    rprotect(src_ends = RSaneAllocVector(REALSXP, chain_intervs.size()));
    rprotect(src_chroms = RSaneAllocVector(STRSXP, num_src_chroms));
    rprotect(src_strands = RSaneAllocVector(INTSXP, chain_intervs.size()));
    rprotect(tgt_strands = RSaneAllocVector(INTSXP, chain_intervs.size()));
    rprotect(chain_ids = RSaneAllocVector(REALSXP, chain_intervs.size()));
    rprotect(scores = RSaneAllocVector(REALSXP, chain_intervs.size()));

	for (ChainIntervals::const_iterator iinterval = chain_intervs.begin(); iinterval != chain_intervs.end(); ++iinterval) {
		INTEGER(src_chroms_idx)[iinterval - chain_intervs.begin()] = iinterval->chromid_src + 1;
		REAL(src_starts)[iinterval - chain_intervs.begin()] = iinterval->start_src;
		REAL(src_ends)[iinterval - chain_intervs.begin()] = iinterval->end_src;
		// Convert 0/1 to +1/-1 for R output
		INTEGER(src_strands)[iinterval - chain_intervs.begin()] = (iinterval->strand_src == 0) ? 1 : -1;
		INTEGER(tgt_strands)[iinterval - chain_intervs.begin()] = (iinterval->strand == 0) ? 1 : -1;
		REAL(chain_ids)[iinterval - chain_intervs.begin()] = iinterval->chain_id;
		REAL(scores)[iinterval - chain_intervs.begin()] = iinterval->score;
	}

	for (unsigned id = 0; id < num_src_chroms; ++id)
		SET_STRING_ELT(src_chroms, id, Rf_mkChar(src_id2chrom[id].c_str()));

    for (int i = 0; i < ChainInterval::NUM_COLS; i++)
		SET_STRING_ELT(col_names, i, Rf_mkChar(ChainInterval::COL_NAMES[i]));

    Rf_setAttrib(src_chroms_idx, R_LevelsSymbol, src_chroms);
    Rf_setAttrib(src_chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(answer, ChainInterval::STRAND, tgt_strands);
    SET_VECTOR_ELT(answer, ChainInterval::CHROM_SRC, src_chroms_idx);
    SET_VECTOR_ELT(answer, ChainInterval::START_SRC, src_starts);
    SET_VECTOR_ELT(answer, ChainInterval::END_SRC, src_ends);
    SET_VECTOR_ELT(answer, ChainInterval::STRAND_SRC, src_strands);
    SET_VECTOR_ELT(answer, ChainInterval::CHAIN_ID, chain_ids);
    SET_VECTOR_ELT(answer, ChainInterval::SCORE, scores);

    runprotect(7);
    return answer;
}

