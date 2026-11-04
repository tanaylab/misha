#include <string>
#include <unistd.h>
#include <fcntl.h>
#include <cstdint>
#include <sys/stat.h>
#include <sys/types.h>

#include "GIntervalsBigSet1D.h"
#include "IntervalsIndex1D.h"
#include "rdbutils.h"

//------------------------------------- GIntervalsBigSet1D --------------------------------------

// Static members initialization
std::map<std::string, std::shared_ptr<IntervalsIndex1D>> GIntervalsBigSet1D::s_index_cache;
std::map<std::string, std::string> GIntervalsBigSet1D::s_index_stamp;
std::mutex GIntervalsBigSet1D::s_cache_mutex;

// The stamp of an index file: its size and its header, which holds the checksum of its entries (the
// offsets and lengths read through it); "" if it cannot be read. One that another session wrote
// since it was cached with another layout differs. (Not the file's inode and mtime: a file deleted
// and written again within the same second can get the same ones.)
static std::string index_file_stamp(const std::string &path)
{
	char header[64];
	struct stat st;
	int fd = open(path.c_str(), O_RDONLY);
	if (fd < 0)
		return std::string();
	ssize_t n = fstat(fd, &st) ? -1 : read(fd, header, sizeof(header));
	close(fd);
	if (n < 0)
		return std::string();
	return std::to_string((long long)st.st_size) + ":" + std::string(header, n);
}

std::shared_ptr<IntervalsIndex1D> GIntervalsBigSet1D::get_intervals_index(const std::string &intervset_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	auto it = s_index_cache.find(intervset_dir);
	if (it != s_index_cache.end()) {
		return it->second;
	}

	// Load index
	auto idx = std::make_shared<IntervalsIndex1D>();
	std::string idx_path = intervset_dir + "/intervals.idx";
	// stamped before the load: a file replaced in between is seen as replaced
	std::string stamp = index_file_stamp(idx_path);
	if (!idx->load(idx_path)) {
		// Index doesn't exist - return nullptr
		return nullptr;
	}

	// Cache it
	s_index_cache[intervset_dir] = idx;
	if (!stamp.empty())
		s_index_stamp[intervset_dir] = stamp;
	return idx;
}

std::string GIntervalsBigSet1D::get_intervset_dir(const std::string &intervset, SEXP envir) {
	return interv2path(envir, intervset.c_str());
}

void GIntervalsBigSet1D::invalidate_index_cache(const std::string &intervset_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	s_index_cache.erase(intervset_dir);
	s_index_stamp.erase(intervset_dir);
}

void GIntervalsBigSet1D::clear_index_cache() {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	s_index_cache.clear();
	s_index_stamp.clear();
}

bool GIntervalsBigSet1D::drop_stale_index(const std::string &intervset_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	if (s_index_cache.find(intervset_dir) == s_index_cache.end())
		return false;
	auto it = s_index_stamp.find(intervset_dir);
	if (it != s_index_stamp.end() && index_file_stamp(intervset_dir + "/intervals.idx") == it->second)
		return false;
	s_index_cache.erase(intervset_dir);
	s_index_stamp.erase(intervset_dir);
	return true;
}

void GIntervalsBigSet1D::init(const char *intervset, SEXP meta, const IntervUtils &iu)
{
	GIntervalsBigSet::init(intervset, iu);
	GIntervalsMeta1D::init(intervset, meta, iu.get_chromkey());

	m_cur_chromid = m_chrom2size.size();
	m_iter_chrom = -1;
	m_iter_index = -1;
	m_iter_chrom_index = 0;
	m_do_sort = false;
	m_do_unify_overlaps = false;
	m_iinterval = m_intervals.end();

	// Initialize smart handle state
	m_dat_fp = nullptr;
	m_dat_open = false;
}

void GIntervalsBigSet1D::load_chrom(int chromid)
{
	m_iter_chrom_index = 0;
	if (get_num_intervals(chromid)) {
		if (m_intervals.empty() || m_intervals.front().chromid != chromid) {
			// Construct per-chromosome filename
			string intervset_dir = interv2path(m_iu->get_env(), m_intervset);
			const GenomeTrack::ChromFiles &files = m_iu->get_chrom_files(intervset_dir);

			// SURGICAL PIVOT: Check for per-chromosome file first
			SEXP rintervals = R_NilValue;

			if (files.exists(chromid)) {
				// PER-CHROMOSOME PATH: Per-chromosome file exists
				rintervals = rprotect_ptr(RSaneUnserialize(files.path(chromid).c_str()));
			} else {
				// INDEXED PATH: Use intervals.dat with index
				// An index that another session replaced since this process cached it would be read at
				// the old offsets: checked when the data file is opened, and after a size or unserialize
				// error the read is done once more, with the data file reopened (as
				// GenomeTrackSparse::init_read)
				auto read_indexed = [&]() -> SEXP {
					string dat_path = intervset_dir + "/intervals.dat";

					// Smart handle: open intervals.dat once, reuse across chromosomes
					if (!m_dat_open || m_dat_path != dat_path) {
						// Close previous file if open
						if (m_dat_open && m_dat_fp) {
							fclose(m_dat_fp);
							m_dat_fp = nullptr;
						}

						// Open new file (and the cached index checked against its file: one stat, once per call)
						drop_stale_index(intervset_dir);
						m_dat_fp = fopen(dat_path.c_str(), "rb");
						if (!m_dat_fp) {
							verror("Cannot open indexed intervals file %s: %s",
								   dat_path.c_str(), strerror(errno));
						}
						m_dat_open = true;
						m_dat_path = dat_path;
					}

					// Get index and lookup entry
					auto idx = get_intervals_index(intervset_dir);
					if (!idx) {
						verror("Intervals index not found for %s", m_intervset.c_str());
					}

					const IntervalsContigEntry* entry = idx->get_entry(chromid);
					if (!entry || entry->length == 0) {
						// Empty chromosome
						return R_NilValue;
					}

					struct stat dat_st;
					if (fstat(fileno(m_dat_fp), &dat_st) || entry->offset + entry->length > (uint64_t)dat_st.st_size)
						verror("Index entry of %s extends past the end of %s", m_intervset.c_str(), dat_path.c_str());
					// Seek to the chromosome's data in intervals.dat
					if (fseek(m_dat_fp, entry->offset, SEEK_SET) != 0) {
						verror("Failed to seek in %s: %s", dat_path.c_str(), strerror(errno));
					}

					// Read serialized R object from current position
					return RSaneUnserialize(m_dat_fp);
				};
				SEXP res;
				try {
					res = read_indexed();
				} catch (TGLException &) {
					if (!drop_stale_index(intervset_dir) && !m_dat_open)
						throw;
					if (m_dat_fp)
						fclose(m_dat_fp);
					m_dat_fp = nullptr;
					m_dat_open = false;
					res = read_indexed();
				}
				if (res == R_NilValue) {
					// Empty chromosome
					m_intervals.clear();
					return;
				}
				rintervals = rprotect_ptr(res);
			}

			// Convert R intervals to C++ intervals
			m_iu->convert_rintervs(rintervals, &m_intervals, NULL);
			// Name it: convert_rintervs protects on some of its input shapes, so
			// "the last one protected" is not necessarily rintervals.
			runprotect(rintervals);

			// set udata
			uint64_t offset = m_orig_chrom2offset[chromid];
			for (GIntervals::iterator iinterval = m_intervals.begin(); iinterval < m_intervals.end(); ++iinterval)
				iinterval->udata = (void *)(intptr_t)(iinterval - m_intervals.begin() + offset);

			if (m_do_sort)
				m_intervals.sort(m_compare);

			if (m_do_unify_overlaps)
				m_intervals.unify_overlaps(m_unify_touching_intervals);
		}
	} else
		m_intervals.clear();
}

void GIntervalsBigSet1D::begin_save(const char *intervset, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	string path = interv2path(iu.get_env(), intervset);
	if (mkdir(path.c_str(), 0777))
		verror("Cannot create intervals directory at %s: %s", path.c_str(), strerror(errno));

	init_chromstats(chromstats, iu);
}

void GIntervalsBigSet1D::save_chrom_plain_intervals(const char *intervset, GIntervals &intervals, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	if (intervals.size()) {
		SEXP rintervals = iu.convert_intervs(&intervals);
		save_chrom(intervset, &intervals, rintervals, iu, chromstats);
		// Every caller runs this once per contig. convert_intervs' result has to be
		// given back or the protect stack grows without bound across the genome.
		runprotect(rintervals);
		intervals.clear();
	}
}

void GIntervalsBigSet1D::save_chrom(const char *intervset, GIntervalsFetcher1D *intervals, SEXP rintervals, const IntervUtils &iu, vector<ChromStat> &chromstats)
{
	if (!intervals->size()) 
		return;

	pair<int, ChromStat> res = get_chrom_stat(intervals);
	int &chromid = res.first;
	ChromStat &chromstat = res.second;
	chromstats[chromid] = chromstat;

	string filename = interv2path(iu.get_env(), intervset);
	filename += "/";
	filename += iu.id2chrom(chromid);
	RSaneSerialize(rintervals, filename.c_str());
}

void GIntervalsBigSet1D::end_save_plain_intervals(const char *intervset, const IntervUtils &iu, const vector<ChromStat> &chromstats)
{
	save_plain_intervals_meta(interv2path(iu.get_env(), intervset).c_str(), chromstats, iu);
}

void GIntervalsBigSet1D::end_save(const char *intervset, SEXP zeroline, const IntervUtils &iu, const vector<ChromStat> &chromstats)
{
	save_meta(interv2path(iu.get_env(), intervset).c_str(), zeroline, chromstats, iu);
}

GIntervalsFetcher1D *GIntervalsBigSet1D::create_masked_copy(const set<int> &chromids_mask) const
{
	GIntervalsBigSet1D *obj = new GIntervalsBigSet1D();

	init_masked_copy(obj, chromids_mask);

	obj->m_intervset = m_intervset;
	obj->m_iu = m_iu;
	obj->m_cur_chromid = obj->m_chrom2size.size();
	obj->m_iter_chrom = -1;
	obj->m_iter_index = -1;
	obj->m_iter_chrom_index = 0;
	obj->m_do_sort = false;
	obj->m_do_unify_overlaps = false;
	obj->m_iinterval = obj->m_intervals.end();

	if (m_do_sort)
		obj->sort(m_compare);

	if (m_do_unify_overlaps) 
		obj->unify_overlaps(m_unify_touching_intervals);

	return obj;
}

void GIntervalsBigSet1D::begin_iter()
{
	m_iter_chrom = -1;
	m_iter_index = 0;
	m_iter_chrom_index = 0;
	m_intervals.clear();
	for (m_cur_chromid = 0; m_cur_chromid < (int)m_chrom2size.size(); ++m_cur_chromid) {
		if (get_num_intervals(m_cur_chromid)) {
			load_chrom(m_cur_chromid);
			m_iinterval = m_intervals.begin();
			return;
		}
	}
}

void GIntervalsBigSet1D::begin_chrom_iter(int chromid)
{
	m_iter_chrom = chromid;
	m_iter_index = 0;
	m_iter_chrom_index = 0;
	for (m_cur_chromid = 0; m_cur_chromid < (int)m_chrom2size.size(); ++m_cur_chromid) {
		if (m_cur_chromid == chromid) {
			if (get_num_intervals(m_cur_chromid)) {
				load_chrom(m_cur_chromid);
				m_iinterval = m_intervals.begin();
				return;
			}
			break;
		}
		m_iter_index += get_num_intervals(m_cur_chromid);
	}
	m_intervals.clear();
	m_iinterval = m_intervals.end();
}

void GIntervalsBigSet1D::sort(Compare_t compare)
{
	m_do_sort = true;
	m_compare = compare;
	if (m_intervals.size())
		m_intervals.sort(m_compare); 
}

void GIntervalsBigSet1D::unify_overlaps(bool unify_touching_intervals)
{
	if (m_do_unify_overlaps && m_unify_touching_intervals == unify_touching_intervals) 
		return;

	m_do_unify_overlaps = true;
	m_unify_touching_intervals = unify_touching_intervals;

	m_size = 0;
	m_range = 0;
	if (m_unify_touching_intervals) {
		m_user_chrom2size = &m_chrom2unified_touching_size;
		for (vector<int64_t>::const_iterator isize = m_chrom2unified_touching_size.begin(); isize < m_chrom2unified_touching_size.end(); ++isize)
			m_size += *isize;
	} else {
		m_user_chrom2size = &m_chrom2unified_overlap_size;
		for (vector<int64_t>::const_iterator isize = m_chrom2unified_overlap_size.begin(); isize < m_chrom2unified_overlap_size.end(); ++isize) 
			m_size += *isize;
	}
	for (vector<int64_t>::const_iterator irange = m_chrom2unified_overlap_range.begin(); irange < m_chrom2unified_overlap_range.end(); ++irange) 
		m_range += *irange;

	if (m_intervals.size())
		m_intervals.unify_overlaps(m_unify_touching_intervals);
}

void GIntervalsBigSet1D::verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix) const
{
	if (m_contains_overlaps)
		TGLError<GIntervalsFetcher1D>(OVERLAPPING_INTERVAL, "%sIntervals set %s contains overlapping intervals", error_prefix, m_intervset.c_str());
}

