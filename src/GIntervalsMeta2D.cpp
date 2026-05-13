#include <algorithm>
#include <cstdint>
#include "GIntervalsMeta2D.h"
#include "rdbutils.h"

const char *GIntervalsMeta2D::STAT_COL_NAMES[NUM_STAT_COLS] = {
	"chrom1", "chrom2", "contains_overlaps", "size", "surface"
};

void GIntervalsMeta2D::rebuild_pair_index()
{
	m_pair_keys_sorted.clear();
	m_pair_keys_sorted.reserve(m_pair_stats.size());
	for (const auto &kv : m_pair_stats)
		m_pair_keys_sorted.push_back(kv.first);
	std::sort(m_pair_keys_sorted.begin(), m_pair_keys_sorted.end());

	m_orig_size_prefix.assign(m_pair_keys_sorted.size() + 1, 0);
	for (size_t i = 0; i < m_pair_keys_sorted.size(); ++i) {
		const PairStat &ps = m_pair_stats[m_pair_keys_sorted[i]];
		m_orig_size_prefix[i + 1] = m_orig_size_prefix[i] + ps.orig_size;
	}
}

void GIntervalsMeta2D::init(const char *name, SEXP meta, const GenomeChromKey &chromkey)
{
	if (!is2d(meta) || !Rf_isVector(meta) || Rf_length(meta) < 1) {
		verror("%s: Invalid format of .meta file", name);
	}

	m_chromkey = (GenomeChromKey *)&chromkey;
	m_size = 0;
	m_surface = 0;
	m_pair_stats.clear();
	m_pair_keys_sorted.clear();
	m_orig_size_prefix.clear();

	SEXP stat = VECTOR_ELT(meta, 0);
	SEXP colnames = Rf_getAttrib(stat, R_NamesSymbol);

	if (Rf_length(stat) != NUM_STAT_COLS || !Rf_isString(colnames) || Rf_length(colnames) != NUM_STAT_COLS || strcmp(CHAR(STRING_ELT(colnames, 0)), STAT_COL_NAMES[0]))
		verror("%s: Invalid format of .meta file", name);

	for (int i = 1; i < NUM_STAT_COLS; ++i) {
		if (Rf_length(VECTOR_ELT(stat, i - 1)) != Rf_length(VECTOR_ELT(stat, i)) || strcmp(CHAR(STRING_ELT(colnames, i)), STAT_COL_NAMES[i]))
			verror("%s: Invalid format of .meta file", name);
	}

	SEXP chroms1 = VECTOR_ELT(stat, CHROM1_COL);
	SEXP chrom_levels1 = Rf_getAttrib(chroms1, R_LevelsSymbol);
	SEXP chroms2 = VECTOR_ELT(stat, CHROM2_COL);
	SEXP chrom_levels2 = Rf_getAttrib(chroms2, R_LevelsSymbol);
	SEXP sizes = VECTOR_ELT(stat, SIZE_COL);
	SEXP surfaces = VECTOR_ELT(stat, SURFACE_COL);
	SEXP contains_overlaps = VECTOR_ELT(stat, CONTAINS_OVERLAPS_COL);

	m_pair_stats.reserve(Rf_length(sizes));

	for (int i = 0; i < Rf_length(sizes); ++i) {
		const char *chrom1 = Rf_isString(chroms1) ? CHAR(STRING_ELT(chroms1, i)) : CHAR(STRING_ELT(chrom_levels1, INTEGER(chroms1)[i] - 1));
		const char *chrom2 = Rf_isString(chroms2) ? CHAR(STRING_ELT(chroms2, i)) : CHAR(STRING_ELT(chrom_levels2, INTEGER(chroms2)[i] - 1));
		int chromid1 = m_chromkey->chrom2id(chrom1);
		int chromid2 = m_chromkey->chrom2id(chrom2);
		int64_t size = (int64_t)(Rf_isReal(sizes) ? REAL(sizes)[i] : INTEGER(sizes)[i]);
		double surface = REAL(surfaces)[i];

		if (size == 0)
			continue; // skip phantom empty rows; sparse store keeps only populated pairs

		PairStat &ps = m_pair_stats[pair_key(chromid1, chromid2)];
		ps.size = size;
		ps.orig_size = size;
		ps.surface = surface;
		ps.contains_overlaps = LOGICAL(contains_overlaps)[i];

		m_size += (uint64_t)size;
		m_surface += surface;
	}

	rebuild_pair_index();
}

void GIntervalsMeta2D::init_masked_copy(GIntervalsMeta2D *obj, const set<ChromPair> &chrompairs_mask) const
{
	obj->m_chromkey = m_chromkey;
	obj->m_size = 0;
	obj->m_surface = 0;
	obj->m_pair_stats.clear();
	obj->m_pair_keys_sorted.clear();
	obj->m_orig_size_prefix.clear();

	// Preserve orig_size for every populated pair (masked or not) so udata
	// offset computations remain consistent with the unmasked parent.
	obj->m_pair_stats.reserve(m_pair_stats.size());
	for (const auto &kv : m_pair_stats) {
		PairStat ps_copy;
		ps_copy.orig_size = kv.second.orig_size;

		ChromPair cp(key_chrom1(kv.first), key_chrom2(kv.first));
		if (chrompairs_mask.find(cp) != chrompairs_mask.end()) {
			ps_copy.size = kv.second.size;
			ps_copy.surface = kv.second.surface;
			ps_copy.contains_overlaps = kv.second.contains_overlaps;
			obj->m_size += (uint64_t)kv.second.size;
			obj->m_surface += kv.second.surface;
		}
		obj->m_pair_stats[kv.first] = ps_copy;
	}

	obj->rebuild_pair_index();
}

void GIntervalsMeta2D::ChromStats2D::pack(void *&ptr) const
{
	// Layout: [uint32_t count][PackedEntry * count]. count must not exceed
	// m_max_pairs (set by begin_save) - the caller pre-allocated the buffer.
	uint32_t count = (uint32_t)m_stats.size();
	if (m_max_pairs && count > m_max_pairs)
		verror("ChromStats2D: kid produced %u populated chrom-pairs but only %zu were reserved", (unsigned)count, m_max_pairs);
	memcpy(ptr, &count, sizeof(count));
	ptr = (char *)ptr + sizeof(count);
	for (const auto &kv : m_stats) {
		PackedEntry e;
		e.chromid1 = (uint32_t)key_chrom1(kv.first);
		e.chromid2 = (uint32_t)key_chrom2(kv.first);
		e.stat = kv.second;
		memcpy(ptr, &e, sizeof(e));
		ptr = (char *)ptr + sizeof(e);
	}
}

void GIntervalsMeta2D::ChromStats2D::unpack(void *&ptr)
{
	uint32_t count = 0;
	memcpy(&count, ptr, sizeof(count));
	ptr = (char *)ptr + sizeof(count);
	for (uint32_t i = 0; i < count; ++i) {
		PackedEntry e;
		memcpy(&e, ptr, sizeof(e));
		ptr = (char *)ptr + sizeof(e);
		// Skip empty entries (defensive); set() also skips them.
		if (!e.stat.size) continue;
		m_stats[pack_key((int)e.chromid1, (int)e.chromid2)] = e.stat;
	}
}

void GIntervalsMeta2D::save_plain_intervals_meta(const char *path, const ChromStats2D &chromstats, const IntervUtils &iu)
{
	GIntervals2D intervals;
	SEXP zeroline = iu.convert_intervs(&intervals, GInterval2D::NUM_COLS, false);
	save_meta(path, zeroline, chromstats, iu);
}

void GIntervalsMeta2D::save_meta(const char *path, SEXP zeroline, const ChromStats2D &chromstats, const IntervUtils &iu)
{
	SEXP rstat;
	SEXP colnames;
	SEXP rownames;
	SEXP chroms1, chroms2, chroms_idx1, chroms_idx2, rsize, rsurface, roverlaps;

    rstat = rprotect_ptr(RSaneAllocVector(VECSXP, NUM_STAT_COLS));
    colnames = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_STAT_COLS));
    chroms1 = rprotect_ptr(RSaneAllocVector(STRSXP, iu.get_chromkey().get_num_chroms()));
    chroms2 = rprotect_ptr(RSaneAllocVector(STRSXP, iu.get_chromkey().get_num_chroms()));

	for (int i = 0; i < NUM_STAT_COLS; i++)
		SET_STRING_ELT(colnames, i, Rf_mkChar(STAT_COL_NAMES[i]));

	// Collect populated (chromid1, chromid2, stat) tuples and sort by
	// (chromid1, chromid2) so the saved .meta has stable row order matching
	// the prior dense-walk implementation.
	struct Row { int c1; int c2; ChromStat stat; };
	vector<Row> rows;
	rows.reserve(chromstats.size());
	for (auto it = chromstats.begin(); it != chromstats.end(); ++it) {
		if (!it->second.size)
			continue;
		Row r;
		r.c1 = ChromStats2D::key_chrom1(it->first);
		r.c2 = ChromStats2D::key_chrom2(it->first);
		r.stat = it->second;
		rows.push_back(r);
	}
	std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
		if (a.c1 != b.c1) return a.c1 < b.c1;
		return a.c2 < b.c2;
	});

	int num_nonempty_chroms = (int)rows.size();

    chroms_idx1 = rprotect_ptr(RSaneAllocVector(INTSXP, num_nonempty_chroms));
    chroms_idx2 = rprotect_ptr(RSaneAllocVector(INTSXP, num_nonempty_chroms));
    rsize = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    rsurface = rprotect_ptr(RSaneAllocVector(REALSXP, num_nonempty_chroms));
    roverlaps = rprotect_ptr(RSaneAllocVector(LGLSXP, num_nonempty_chroms));
    rownames = rprotect_ptr(RSaneAllocVector(INTSXP, num_nonempty_chroms));

    for (unsigned id = 0; id < (unsigned)iu.get_chromkey().get_num_chroms(); ++id) {
		SET_STRING_ELT(chroms1, id, Rf_mkChar(iu.id2chrom(id).c_str()));
		SET_STRING_ELT(chroms2, id, Rf_mkChar(iu.id2chrom(id).c_str()));
	}

	for (int res_index = 0; res_index < num_nonempty_chroms; ++res_index) {
		const Row &r = rows[res_index];
		INTEGER(chroms_idx1)[res_index] = r.c1 + 1;
		INTEGER(chroms_idx2)[res_index] = r.c2 + 1;
		REAL(rsize)[res_index] = r.stat.size;
		REAL(rsurface)[res_index] = r.stat.surface;
		LOGICAL(roverlaps)[res_index] = r.stat.contains_overlaps;
		INTEGER(rownames)[res_index] = res_index + 1;
	}

    Rf_setAttrib(rstat, R_RowNamesSymbol, rownames);
    Rf_setAttrib(chroms_idx1, R_LevelsSymbol, chroms1);
    Rf_setAttrib(chroms_idx2, R_LevelsSymbol, chroms2);
    Rf_setAttrib(chroms_idx1, R_ClassSymbol, Rf_mkString("factor"));
    Rf_setAttrib(chroms_idx2, R_ClassSymbol, Rf_mkString("factor"));

    SET_VECTOR_ELT(rstat, CHROM1_COL, chroms_idx1);
    SET_VECTOR_ELT(rstat, CHROM2_COL, chroms_idx2);
    SET_VECTOR_ELT(rstat, SIZE_COL, rsize);
    SET_VECTOR_ELT(rstat, SURFACE_COL, rsurface);
    SET_VECTOR_ELT(rstat, CONTAINS_OVERLAPS_COL, roverlaps);

    Rf_setAttrib(rstat, R_NamesSymbol, colnames);
    Rf_setAttrib(rstat, R_ClassSymbol, Rf_mkString("data.frame"));

	GIntervalsMeta::save_meta(path, rstat, zeroline);
}
