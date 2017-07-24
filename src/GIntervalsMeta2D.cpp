#include "GIntervalsMeta2D.h"
#include "rdbutils.h"

const char *GIntervalsMeta2D::STAT_COL_NAMES[NUM_STAT_COLS] = {
	"chrom1", "chrom2", "contains_overlaps", "size", "surface"
};

void GIntervalsMeta2D::init(const char *name, SEXP meta, const GenomeChromKey &chromkey)
{
	if (!is2d(meta) || !isVector(meta) || length(meta) < 1) 
		verror("%s: Invalid format of .meta file", name);

	m_chromkey = (GenomeChromKey *)&chromkey;
	m_size = 0;
	m_surface = 0;
	m_chroms2size.clear();
	m_contains_overlaps.clear();
	m_surfaces.clear();
	m_chroms2size.resize(m_chromkey->get_num_chroms() * m_chromkey->get_num_chroms(), 0);
	m_contains_overlaps.resize(m_chromkey->get_num_chroms() * m_chromkey->get_num_chroms(), false);
	m_surfaces.resize(m_chromkey->get_num_chroms() * m_chromkey->get_num_chroms(), 0);

	SEXP stat = VECTOR_ELT(meta, 0);
	SEXP colnames = getAttrib(stat, R_NamesSymbol);

	if (length(stat) != NUM_STAT_COLS || !isString(colnames) || length(colnames) != NUM_STAT_COLS || strcmp(CHAR(STRING_ELT(colnames, 0)), STAT_COL_NAMES[0]))
		verror("%s: Invalid format of .meta file", name);

	for (int i = 1; i < NUM_STAT_COLS; ++i) {
		if (length(VECTOR_ELT(stat, i - 1)) != length(VECTOR_ELT(stat, i)) || strcmp(CHAR(STRING_ELT(colnames, i)), STAT_COL_NAMES[i]))
			verror("%s: Invalid format of .meta file", name);
	}

	SEXP chroms1 = VECTOR_ELT(stat, CHROM1_COL);
	SEXP chrom_levels1 = getAttrib(chroms1, R_LevelsSymbol);
	SEXP chroms2 = VECTOR_ELT(stat, CHROM2_COL);
	SEXP chrom_levels2 = getAttrib(chroms2, R_LevelsSymbol);
	SEXP sizes = VECTOR_ELT(stat, SIZE_COL);
	SEXP surfaces = VECTOR_ELT(stat, SURFACE_COL);
	SEXP contains_overlaps = VECTOR_ELT(stat, CONTAINS_OVERLAPS_COL);

	for (int i = 0; i < length(sizes); ++i) {
		const char *chrom1 = isString(chroms1) ? CHAR(STRING_ELT(chroms1, i)) : CHAR(STRING_ELT(chrom_levels1, INTEGER(chroms1)[i] - 1));
		const char *chrom2 = isString(chroms2) ? CHAR(STRING_ELT(chroms2, i)) : CHAR(STRING_ELT(chrom_levels2, INTEGER(chroms2)[i] - 1));
		int chromid1 = m_chromkey->chrom2id(chrom1);
		int chromid2 = m_chromkey->chrom2id(chrom2);
		int64_t size = (int64_t)(isReal(sizes) ? REAL(sizes)[i] : INTEGER(sizes)[i]);
		double surface = REAL(surfaces)[i];
		int idx = chroms2idx(chromid1, chromid2);

		m_chroms2size[idx] = size;
		m_surfaces[idx] = surface;
		m_contains_overlaps[idx] = LOGICAL(contains_overlaps)[i];
		m_size += (size_t)size;
		m_surface += surface;
	}

	m_orig_chroms2size = m_chroms2size;
}

void GIntervalsMeta2D::init_masked_copy(GIntervalsMeta2D *obj, const set<ChromPair> &chrompairs_mask) const
{
	obj->m_chromkey = m_chromkey;
	obj->m_size = 0;
	obj->m_surface = 0;
	obj->m_chroms2size.clear();
	obj->m_contains_overlaps.clear();
	obj->m_surfaces.clear();
	obj->m_chroms2size.resize(m_chroms2size.size(), 0);
	obj->m_contains_overlaps.resize(m_contains_overlaps.size(), false);
	obj->m_surfaces.resize(m_surfaces.size(), 0);
	obj->m_orig_chroms2size = m_orig_chroms2size;

	for (int chromid = 0; chromid < obj->m_chroms2size.size(); ++chromid) {
		int chromid1 = idx2chrom1(chromid);
		int chromid2 = idx2chrom2(chromid);

		if (chrompairs_mask.find(ChromPair(chromid1, chromid2)) == chrompairs_mask.end())
			continue;

		obj->m_chroms2size[chromid] = m_chroms2size[chromid];
		obj->m_contains_overlaps[chromid] = m_contains_overlaps[chromid];
		obj->m_surfaces[chromid] = m_surfaces[chromid];
		obj->m_size += (size_t)m_chroms2size[chromid];
		obj->m_surface += m_surfaces[chromid];
	}
}

void GIntervalsMeta2D::init_chromstats(vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	chromstats.clear();
	chromstats.resize(iu.get_chromkey().get_num_chroms() * iu.get_chromkey().get_num_chroms());
}

void GIntervalsMeta2D::save_plain_intervals_meta(const char *path, const vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	GIntervals2D intervals;
	SEXP zeroline = iu.convert_intervs(&intervals, GInterval2D::NUM_COLS, false);
	save_meta(path, zeroline, chromstats, iu);
}

void GIntervalsMeta2D::save_meta(const char *path, SEXP zeroline, const vector<ChromStat> &chromstats, const IntervUtils &iu)
{
	SEXP rstat;
	SEXP colnames;
	SEXP rownames;
	SEXP chroms1, chroms2, chroms_idx1, chroms_idx2;

	rprotect(rstat = allocVector(VECSXP, NUM_STAT_COLS));

	setAttrib(rstat, R_NamesSymbol, (colnames = allocVector(STRSXP, NUM_STAT_COLS)));
	setAttrib(rstat, R_ClassSymbol, mkString("data.frame"));

	for (int i = 0; i < NUM_STAT_COLS; i++)
		SET_STRING_ELT(colnames, i, mkChar(STAT_COL_NAMES[i]));

	int num_nonempty_chroms = 0;
	for (vector<ChromStat>::const_iterator ichromstat = chromstats.begin(); ichromstat != chromstats.end(); ++ichromstat) {
		if (ichromstat->size) 
			++num_nonempty_chroms;
	}

	SET_VECTOR_ELT(rstat, CHROM1_COL, (chroms_idx1 = allocVector(INTSXP, num_nonempty_chroms)));
	SET_VECTOR_ELT(rstat, CHROM2_COL, (chroms_idx2 = allocVector(INTSXP, num_nonempty_chroms)));
	SET_VECTOR_ELT(rstat, SIZE_COL, allocVector(REALSXP, num_nonempty_chroms));
	SET_VECTOR_ELT(rstat, SURFACE_COL, allocVector(REALSXP, num_nonempty_chroms));
	SET_VECTOR_ELT(rstat, CONTAINS_OVERLAPS_COL, allocVector(LGLSXP, num_nonempty_chroms));

	setAttrib(rstat, R_RowNamesSymbol, (rownames = allocVector(INTSXP, num_nonempty_chroms)));
	setAttrib(chroms_idx1, R_LevelsSymbol, (chroms1 = allocVector(STRSXP, iu.get_chromkey().get_num_chroms())));
	setAttrib(chroms_idx2, R_LevelsSymbol, (chroms2 = allocVector(STRSXP, iu.get_chromkey().get_num_chroms())));
	setAttrib(chroms_idx1, R_ClassSymbol, mkString("factor"));
	setAttrib(chroms_idx2, R_ClassSymbol, mkString("factor"));

	for (unsigned id = 0; id < (unsigned)iu.get_chromkey().get_num_chroms(); ++id) {
		SET_STRING_ELT(chroms1, id, mkChar(iu.id2chrom(id).c_str()));
		SET_STRING_ELT(chroms2, id, mkChar(iu.id2chrom(id).c_str()));
	}

	int res_index = 0;
	for (int chromid1 = 0; chromid1 < iu.get_chromkey().get_num_chroms(); ++chromid1) {
		for (int chromid2 = 0; chromid2 < iu.get_chromkey().get_num_chroms(); ++chromid2) {
			const ChromStat &chromstat = chromstats[chromid1 * iu.get_chromkey().get_num_chroms() + chromid2];

			if (!chromstat.size) 
				continue;

			INTEGER(chroms_idx1)[res_index] = chromid1 + 1;
			INTEGER(chroms_idx2)[res_index] = chromid2 + 1;
			REAL(VECTOR_ELT(rstat, SIZE_COL))[res_index] = chromstat.size;
			REAL(VECTOR_ELT(rstat, SURFACE_COL))[res_index] = chromstat.surface;
			LOGICAL(VECTOR_ELT(rstat, CONTAINS_OVERLAPS_COL))[res_index] = chromstat.contains_overlaps;
			INTEGER(rownames)[res_index] = res_index + 1;
			++res_index;
		}
	}

	GIntervalsMeta::save_meta(path, rstat, zeroline);
}

