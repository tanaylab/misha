#include <cstdint>
#include "GIntervalsMeta1D.h"
#include "rdbutils.h"

const char *GIntervalsMeta1D::STAT_COL_NAMES[NUM_STAT_COLS] = {
	"chrom", "contains_overlaps", "size", "unified_overlap_size", "unified_touching_size",
	"range", "unified_overlap_range"
};

void GIntervalsMeta1D::init(const char *name, SEXP meta, const GenomeChromKey &chromkey)
{
	if (!is1d(meta) || !Rf_isVector(meta) || Rf_length(meta) < 1) 
		verror("%s: Invalid format of .meta file", name);

	m_chromkey = (GenomeChromKey *)&chromkey;
	m_size = 0;
	m_range = 0;
	m_contains_overlaps = false;
	m_user_chrom2size = &m_chrom2size;
	m_chrom2size.clear();
	m_chrom2unified_overlap_size.clear();
	m_chrom2unified_touching_size.clear();
	m_chrom2range.clear();
	m_chrom2unified_overlap_range.clear();
	m_chrom2size.resize(m_chromkey->get_num_chroms(), 0);
	m_chrom2unified_overlap_size.resize(m_chromkey->get_num_chroms(), 0);
	m_chrom2unified_touching_size.resize(m_chromkey->get_num_chroms(), 0);
	m_chrom2range.resize(m_chromkey->get_num_chroms(), 0);
	m_chrom2unified_overlap_range.resize(m_chromkey->get_num_chroms(), 0);

	SEXP stat = VECTOR_ELT(meta, 0);
	SEXP colnames = Rf_getAttrib(stat, R_NamesSymbol);

	if (Rf_length(stat) != NUM_STAT_COLS || !Rf_isString(colnames) || Rf_length(colnames) != NUM_STAT_COLS || strcmp(CHAR(STRING_ELT(colnames, 0)), STAT_COL_NAMES[0]))
		verror("%s: Invalid format of .meta file", name);

	for (int i = 1; i < NUM_STAT_COLS; ++i) {
		if (Rf_length(VECTOR_ELT(stat, i - 1)) != Rf_length(VECTOR_ELT(stat, i)) || strcmp(CHAR(STRING_ELT(colnames, i)), STAT_COL_NAMES[i]))
			verror("%s: Invalid format of .meta file", name);
	}

	SEXP chroms = VECTOR_ELT(stat, CHROM_COL);
	SEXP chrom_levels = Rf_getAttrib(chroms, R_LevelsSymbol);
	SEXP sizes = VECTOR_ELT(stat, SIZE_COL);
	SEXP unified_overlap_sizes = VECTOR_ELT(stat, UNIFIED_OVERLAP_SIZE_COL);
	SEXP unified_touching_sizes = VECTOR_ELT(stat, UNIFIED_TOUCHING_SIZE_COL);
	SEXP ranges = VECTOR_ELT(stat, RANGE_COL);
	SEXP unified_overlap_ranges = VECTOR_ELT(stat, UNIFIED_OVERLAP_RANGE_COL);
	SEXP contains_overlaps = VECTOR_ELT(stat, CONTAINS_OVERLAPS_COL);

	for (int i = 0; i < Rf_length(sizes); ++i) {
		const char *chrom = Rf_isString(chroms) ? CHAR(STRING_ELT(chroms, i)) : CHAR(STRING_ELT(chrom_levels, INTEGER(chroms)[i] - 1));
		int chromid = m_chromkey->chrom2id(chrom);
		int64_t size = (int64_t)(Rf_isReal(sizes) ? REAL(sizes)[i] : INTEGER(sizes)[i]);
		int64_t unified_overlap_size = (int64_t)(Rf_isReal(unified_overlap_sizes) ? REAL(unified_overlap_sizes)[i] : INTEGER(unified_overlap_sizes)[i]);
		int64_t unified_touching_size = (int64_t)(Rf_isReal(unified_touching_sizes) ? REAL(unified_touching_sizes)[i] : INTEGER(unified_touching_sizes)[i]);
		int64_t range = (int64_t)(Rf_isReal(ranges) ? REAL(ranges)[i] : INTEGER(ranges)[i]);
		int64_t unified_overlap_range = (int64_t)(Rf_isReal(unified_overlap_ranges) ? REAL(unified_overlap_ranges)[i] : INTEGER(unified_overlap_ranges)[i]);

		m_chrom2size[chromid] = size;
		m_chrom2unified_overlap_size[chromid] = unified_overlap_size;
		m_chrom2unified_touching_size[chromid] = unified_touching_size;
		m_chrom2range[chromid] = range;
		m_chrom2unified_overlap_range[chromid] = unified_overlap_range;
		m_size += (uint64_t)size;
		m_range += (uint64_t)range;
		m_contains_overlaps |= LOGICAL(contains_overlaps)[i];
	}

	m_orig_chrom2size = m_chrom2size;
}

void GIntervalsMeta1D::init_masked_copy(GIntervalsMeta1D *obj, const set<int> &chromids_mask) const
{
	obj->m_chromkey = m_chromkey;
	obj->m_size = 0;
	obj->m_range = 0;
	obj->m_contains_overlaps = false;
	obj->m_user_chrom2size = &obj->m_chrom2size;
	obj->m_chrom2size.clear();
	obj->m_chrom2unified_overlap_size.clear();
	obj->m_chrom2unified_touching_size.clear();
	obj->m_chrom2range.clear();
	obj->m_chrom2unified_overlap_range.clear();
	obj->m_chrom2size.resize(m_chromkey->get_num_chroms(), 0);
	obj->m_chrom2unified_overlap_size.resize(m_chromkey->get_num_chroms(), 0);
	obj->m_chrom2unified_touching_size.resize(m_chromkey->get_num_chroms(), 0);
	obj->m_chrom2range.resize(m_chromkey->get_num_chroms(), 0);
	obj->m_chrom2unified_overlap_range.resize(m_chromkey->get_num_chroms(), 0);
	obj->m_orig_chrom2size = m_orig_chrom2size;

	for (int chromid = 0; chromid < (int)m_chromkey->get_num_chroms(); ++chromid) {
		if (chromids_mask.find(chromid) == chromids_mask.end()) 
			continue;

		obj->m_chrom2size[chromid] = m_chrom2size[chromid];
		obj->m_chrom2unified_overlap_size[chromid] = m_chrom2unified_overlap_size[chromid];
		obj->m_chrom2unified_touching_size[chromid] = m_chrom2unified_touching_size[chromid];
		obj->m_chrom2range[chromid] = m_chrom2range[chromid];
		obj->m_chrom2unified_overlap_range[chromid] = m_chrom2unified_overlap_range[chromid];
		obj->m_size += (uint64_t)m_chrom2size[chromid];
		obj->m_range += (uint64_t)m_chrom2range[chromid];
		obj->m_contains_overlaps |= m_chrom2size[chromid] != m_chrom2unified_overlap_size[chromid];
	}
}

pair<int, GIntervalsMeta1D::ChromStat> GIntervalsMeta1D::get_chrom_stat(GIntervalsFetcher1D *_intervals)
{
	pair<int, ChromStat> res;
	int &chromid = res.first;

	chromid = -1;

	if (_intervals->size()) {
		if (_intervals->num_chroms() > 1) 
			verror("get_chrom_stat found more than one chromosome in the intervals");

		// we are going to unify overlaps in our intervals => create a copy
		GIntervals intervals;

		intervals.reserve(_intervals->size());
		for (_intervals->begin_iter(); !_intervals->isend(); _intervals->next())
			intervals.push_back(_intervals->cur_interval());

		chromid = intervals.front().chromid;

		ChromStat &chromstat = res.second;
		chromstat.size = intervals.size();
		chromstat.range = intervals.range();
		intervals.unify_overlaps(false);
		chromstat.unified_overlap_size = intervals.size();
		chromstat.unified_overlap_range = intervals.range();
		intervals.unify_overlaps(true);
		chromstat.unified_touching_size = intervals.size();
		chromstat.contains_overlaps = chromstat.size != chromstat.unified_overlap_size;
	}
	return res;
}

void GIntervalsMeta1D::init_chromstats(vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	chromstats.clear();
	chromstats.resize(iu.get_chromkey().get_num_chroms());
}

void GIntervalsMeta1D::save_plain_intervals_meta(const char *path, const vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	GIntervals intervals;
	SEXP zeroline = iu.convert_intervs(&intervals, GInterval::NUM_COLS, false);
	save_meta(path, zeroline, chromstats, iu);
}

void GIntervalsMeta1D::save_meta(const char *path, SEXP zeroline, const vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	// uint64_t num_intervals = 0;
	// for (vector<ChromStat>::const_iterator istat = chromstats.begin(); istat < chromstats.end(); ++istat){
	// 	num_intervals += istat->size;
	// }

	SEXP rstat;
	SEXP colnames;
	SEXP rownames;
	SEXP chroms, chroms_idx;
    SEXP rsize, roverlap_size, rtouching_size, rrange, roverlap_range, roverlaps;

    rstat = rprotect_ptr(RSaneAllocVector(VECSXP, NUM_STAT_COLS));

    Rf_setAttrib(rstat, R_NamesSymbol, (colnames = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_STAT_COLS))));
	Rf_setAttrib(rstat, R_ClassSymbol, Rf_mkString("data.frame"));

	for (int i = 0; i < NUM_STAT_COLS; i++)
		SET_STRING_ELT(colnames, i, Rf_mkChar(STAT_COL_NAMES[i]));

	int num_nonempty_chroms = 0;
	for (vector<ChromStat>::const_iterator ichromstat = chromstats.begin(); ichromstat != chromstats.end(); ++ichromstat) {
		if (ichromstat->size) 
			++num_nonempty_chroms;
	}

    chroms_idx = rprotect_ptr(RSaneAllocVector(INTSXP, num_nonempty_chroms));
    rsize = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    roverlap_size = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    rtouching_size = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    rrange = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    roverlap_range = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    roverlaps = rprotect_ptr(RSaneAllocVector(LGLSXP, num_nonempty_chroms));
    rownames = rprotect_ptr(RSaneAllocVector(INTSXP, num_nonempty_chroms));
    chroms = rprotect_ptr(RSaneAllocVector(STRSXP, iu.get_chromkey().get_num_chroms()));

	for (unsigned id = 0; id < (unsigned)iu.get_chromkey().get_num_chroms(); ++id)
		SET_STRING_ELT(chroms, id, Rf_mkChar(iu.id2chrom(id).c_str()));

	int res_index = 0;
	for (vector<ChromStat>::const_iterator ichromstat = chromstats.begin(); ichromstat != chromstats.end(); ++ichromstat) {
		if (!ichromstat->size) 
			continue;

		INTEGER(chroms_idx)[res_index] = ichromstat - chromstats.begin() + 1;
		REAL(rsize)[res_index] = ichromstat->size;
		REAL(roverlap_size)[res_index] = ichromstat->unified_overlap_size;
		REAL(rtouching_size)[res_index] = ichromstat->unified_touching_size;
		REAL(rrange)[res_index] = ichromstat->range;
		REAL(roverlap_range)[res_index] = ichromstat->unified_overlap_range;
		LOGICAL(roverlaps)[res_index] = ichromstat->contains_overlaps;
		INTEGER(rownames)[res_index] = res_index + 1;
		++res_index;
	}

    Rf_setAttrib(rstat, R_RowNamesSymbol, rownames);
    Rf_setAttrib(chroms_idx, R_LevelsSymbol, chroms);
    Rf_setAttrib(chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(rstat, CHROM_COL, chroms_idx);
    SET_VECTOR_ELT(rstat, SIZE_COL, rsize);
    SET_VECTOR_ELT(rstat, UNIFIED_OVERLAP_SIZE_COL, roverlap_size);
    SET_VECTOR_ELT(rstat, UNIFIED_TOUCHING_SIZE_COL, rtouching_size);
    SET_VECTOR_ELT(rstat, RANGE_COL, rrange);
    SET_VECTOR_ELT(rstat, UNIFIED_OVERLAP_RANGE_COL, roverlap_range);
    SET_VECTOR_ELT(rstat, CONTAINS_OVERLAPS_COL, roverlaps);

	GIntervalsMeta::save_meta(path, rstat, zeroline);
}

