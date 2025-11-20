#include <cstdint>
#include <sys/stat.h>
#include <sys/types.h>

#include "GIntervalsBigSet2D.h"
#include "IntervalsIndex2D.h"
#include "rdbutils.h"

//------------------------------------- GIntervalsBigSet2D --------------------------------------

// Static members initialization
std::map<std::string, std::shared_ptr<IntervalsIndex2D>> GIntervalsBigSet2D::s_index_cache;
std::mutex GIntervalsBigSet2D::s_cache_mutex;

std::shared_ptr<IntervalsIndex2D> GIntervalsBigSet2D::get_intervals_index(const std::string &intervset_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	auto it = s_index_cache.find(intervset_dir);
	if (it != s_index_cache.end()) {
		return it->second;
	}

	// Load index
	auto idx = std::make_shared<IntervalsIndex2D>();
	std::string idx_path = intervset_dir + "/intervals2d.idx";
	if (!idx->load(idx_path)) {
		// Index doesn't exist - return nullptr
		return nullptr;
	}

	// Cache it
	s_index_cache[intervset_dir] = idx;
	return idx;
}

std::string GIntervalsBigSet2D::get_intervset_dir(const std::string &intervset, SEXP envir) {
	return interv2path(envir, intervset.c_str());
}

void GIntervalsBigSet2D::init(const char *intervset, SEXP meta, const IntervUtils &iu)
{
	GIntervalsBigSet::init(intervset, iu);
	GIntervalsMeta2D::init(intervset, meta, iu.get_chromkey());

	if (!is2d(meta))
		verror("Intervals set %s: expecting 1D intervals", intervset);

	m_cur_chromid = m_chroms2size.size();
	m_iter_chromid = -1;
	m_iter_index = -1;
	m_iter_chrom_index = 0;
	m_do_sort = false;
	m_iinterval = m_intervals.end();

	// Initialize smart handle state
	m_dat2d_fp = nullptr;
	m_dat2d_open = false;
}

void GIntervalsBigSet2D::load_chrom(int chromid1, int chromid2)
{
	m_iter_chrom_index = 0;
	if (get_num_intervals(chromid1, chromid2)) {
		if (m_intervals.empty() || m_intervals.front().chromid1() != chromid1 || m_intervals.front().chromid2() != chromid2) {
			// Construct per-chromosome filename
			string intervset_dir = interv2path(m_iu->get_env(), m_intervset);
			string filename = intervset_dir + "/" + m_iu->id2chrom(chromid1) + "-" + m_iu->id2chrom(chromid2);

			// SURGICAL PIVOT: Check for per-chromosome file first
			struct stat st;
			SEXP rintervals = R_NilValue;

			if (stat(filename.c_str(), &st) == 0) {
				// PER-CHROMOSOME PATH: Per-pair file exists
				rintervals = rprotect_ptr(RSaneUnserialize(filename.c_str()));
			} else {
				// INDEXED PATH: Use intervals2d.dat with index
				string dat_path = intervset_dir + "/intervals2d.dat";

				// Smart handle: open intervals2d.dat once, reuse across pairs
				if (!m_dat2d_open || m_dat2d_path != dat_path) {
					// Close previous file if open
					if (m_dat2d_open && m_dat2d_fp) {
						fclose(m_dat2d_fp);
						m_dat2d_fp = nullptr;
					}

					// Open new file
					m_dat2d_fp = fopen(dat_path.c_str(), "rb");
					if (!m_dat2d_fp) {
						verror("Cannot open indexed 2D intervals file %s: %s",
							   dat_path.c_str(), strerror(errno));
					}
					m_dat2d_open = true;
					m_dat2d_path = dat_path;
				}

				// Get index and lookup entry
				auto idx = get_intervals_index(intervset_dir);
				if (!idx) {
					verror("2D intervals index not found for %s", m_intervset.c_str());
				}

				const IntervalsPairEntry* entry = idx->get_entry(chromid1, chromid2);
				if (!entry || entry->length == 0) {
					// Empty chromosome pair
					m_intervals.clear();
					return;
				}

				// Seek to the pair's data in intervals2d.dat
				if (fseek(m_dat2d_fp, entry->offset, SEEK_SET) != 0) {
					verror("Failed to seek in %s: %s", dat_path.c_str(), strerror(errno));
				}

				// Read serialized R object from current position
				rintervals = rprotect_ptr(RSaneUnserialize(m_dat2d_fp));
			}

			// Convert R intervals to C++ intervals
			m_iu->convert_rintervs(rintervals, NULL, &m_intervals);
			runprotect(1);

			// set udata
			uint64_t offset = 0;
			int idx = chroms2idx(chromid1, chromid2);
			for (int i = 0; i < idx; ++i)
				offset += m_orig_chroms2size[i];
			for (GIntervals2D::iterator iinterval = m_intervals.begin(); iinterval < m_intervals.end(); ++iinterval)
				iinterval->udata() = (void *)(intptr_t)(iinterval - m_intervals.begin() + offset);

			if (m_do_sort)
				m_intervals.sort(m_compare);
		}
	} else
		m_intervals.clear();
}

pair<ChromPair, GIntervalsBigSet2D::ChromStat> GIntervalsBigSet2D::get_chrom_stat(GIntervalsFetcher2D *intervals, const IntervUtils &iu)
{
	pair<ChromPair, ChromStat> res(ChromPair(-1, -1), ChromStat());

	if (intervals->size()) {
		if (intervals->num_chrom_pairs() > 1) 
			verror("get_chrom_stat found more than one chromosome pair in the intervals");

		ChromPair &chrompair = res.first;
		intervals->begin_iter();
		chrompair.chromid1 = intervals->cur_interval().chromid1();
		chrompair.chromid2 = intervals->cur_interval().chromid2();

		ChromStat &chromstat = res.second;
		chromstat.size = intervals->size();
		chromstat.surface = intervals->surface();
		try {
			intervals->verify_no_overlaps(iu.get_chromkey());
			chromstat.contains_overlaps = false;
		} catch (TGLException &e) {
			if (e.code() == GIntervalsFetcher2D::OVERLAPPING_INTERVAL) 
				chromstat.contains_overlaps = true;
			else
				throw;
		}
	}
	return res;
}

void GIntervalsBigSet2D::begin_save(const char *intervset, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	string path = interv2path(iu.get_env(), intervset);
	if (mkdir(path.c_str(), 0777))
		verror("Cannot create intervals directory at %s: %s", path.c_str(), strerror(errno));

	init_chromstats(chromstats, iu);
}

void GIntervalsBigSet2D::save_chrom_plain_intervals(const char *intervset, GIntervals2D &intervals, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	if (intervals.size()) {
		SEXP rintervals = iu.convert_intervs(&intervals);
		save_chrom(intervset, &intervals, rintervals, iu, chromstats);
		intervals.clear();
	}
}

void GIntervalsBigSet2D::save_chrom(const char *intervset, GIntervalsFetcher2D *intervals, SEXP rintervals, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	if (!intervals->size()) 
		return;

	pair<ChromPair, ChromStat> res = get_chrom_stat(intervals, iu);
	ChromPair &chrompair = res.first;
	ChromStat &chromstat = res.second;
	chromstats[chrompair.chromid1 * iu.get_chromkey().get_num_chroms() + chrompair.chromid2] = chromstat;

	string filename = interv2path(iu.get_env(), intervset);
	filename += "/";
	filename += iu.id2chrom(chrompair.chromid1);
	filename += "-";
	filename += iu.id2chrom(chrompair.chromid2);
	RSaneSerialize(rintervals, filename.c_str());
}

void GIntervalsBigSet2D::end_save_plain_intervals(const char *intervset, const IntervUtils &iu, const vector<ChromStat> &chromstats)
{
	save_plain_intervals_meta(interv2path(iu.get_env(), intervset).c_str(), chromstats, iu);
}

void GIntervalsBigSet2D::end_save(const char *intervset, SEXP zeroline, const IntervUtils &iu, const vector<ChromStat> &chromstats)
{
	save_meta(interv2path(iu.get_env(), intervset).c_str(), zeroline, chromstats, iu);
}

GIntervalsFetcher2D *GIntervalsBigSet2D::create_masked_copy(const set<ChromPair> &chrompairs_mask) const
{
	GIntervalsBigSet2D *obj = new GIntervalsBigSet2D();

	init_masked_copy(obj, chrompairs_mask);

	obj->m_intervset = m_intervset;
	obj->m_iu = m_iu;
	obj->m_cur_chromid = obj->m_chroms2size.size();
	obj->m_iter_chromid = -1;
	obj->m_iter_index = -1;
	obj->m_iter_chrom_index = 0;
	obj->m_do_sort = false;
	obj->m_intervals.clear();
	obj->m_iinterval = obj->m_intervals.end();
	obj->m_orig_chroms2size = m_orig_chroms2size;

	if (m_do_sort)
		obj->sort(m_compare);

	return obj;
}

void GIntervalsBigSet2D::begin_iter()
{
	m_iter_chromid = -1;
	m_iter_index = 0;
	m_iter_chrom_index = 0;
	m_intervals.clear();
	m_iinterval = m_intervals.end();
	for (m_cur_chromid = 0; m_cur_chromid < (int)m_chroms2size.size(); ++m_cur_chromid) {
		if (m_chroms2size[m_cur_chromid]) {
			int chromid1 = idx2chrom1(m_cur_chromid);
			int chromid2 = idx2chrom2(m_cur_chromid);
			load_chrom(chromid1, chromid2);
			m_iinterval = m_intervals.begin();
			return;
		}
	}
}

void GIntervalsBigSet2D::begin_chrom_iter(int chromid1, int chromid2)
{
	int target_chromid = chroms2idx(chromid1, chromid2);
	m_iter_chromid = target_chromid;
	m_iter_index = 0;
	m_iter_chrom_index = 0;
	for (m_cur_chromid = 0; m_cur_chromid < (int)m_chroms2size.size(); ++m_cur_chromid) {
		if (m_cur_chromid == target_chromid) {
			if (m_chroms2size[m_cur_chromid]) {
				load_chrom(chromid1, chromid2);
				m_iinterval = m_intervals.begin();
			} else {
				m_intervals.clear();
				m_iinterval = m_intervals.end();
			}
			return;
		}
		m_iter_index += m_chroms2size[m_cur_chromid];
	}
	m_intervals.clear();
	m_iinterval = m_intervals.end();
}

void GIntervalsBigSet2D::sort(Compare_t compare)
{
	m_do_sort = true;
	m_compare = compare;
	if (m_intervals.size()) 
		m_intervals.sort(m_compare);
}

void GIntervalsBigSet2D::verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix) const
{
	for (vector<bool>::const_iterator icontains_overlaps = m_contains_overlaps.begin(); icontains_overlaps < m_contains_overlaps.end(); ++icontains_overlaps)  {
		if (*icontains_overlaps) 
			TGLError<GIntervalsFetcher2D>(OVERLAPPING_INTERVAL, "%sIntervals set %s contains overlapping intervals", error_prefix, m_intervset.c_str());
	}
}
