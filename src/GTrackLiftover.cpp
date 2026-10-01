#include <cstdint>
#include "port.h"

#include <cmath>
#include <limits>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "GenomeTrackFixedBin.h"
#include "GenomeTrackIndexedWriter.h"
#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "AggregationHelpers.h"

using namespace std;
using namespace rdb;

//----------------------------------------- BufferedIntervals ---------------------------------------------

class BufferedIntervals {
public:
	BufferedIntervals() : m_last_val(numeric_limits<float>::quiet_NaN()) {}
	BufferedIntervals(const BufferedIntervals &) : m_last_val(numeric_limits<float>::quiet_NaN()) {}

	~BufferedIntervals() { close(); }

	void open(const char *path, const char *mode, int chromid) {
		if (m_bfile.open(path, mode))
			TGLError("Opening file %s: %s", path, strerror(errno));
		m_last_interval.chromid = chromid;
	}

	void close() {
		flush();
		m_bfile.close();
	}

	void write_interval(const GInterval &interval, float val) {
		if (((std::isnan(val) && std::isnan(m_last_val)) || val == m_last_val) && m_last_interval.end == interval.start)
			m_last_interval.end = interval.end;
		else {
			flush();
			m_last_interval = interval;
			m_last_val = val;
		}
	}

	bool read_interval() {
		if (m_bfile.read(&m_last_interval.start, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_interval.end, sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_val, sizeof(float)) != sizeof(float))
		{
			if (m_bfile.eof())
				return false;
			if (m_bfile.error())
				TGLError("Failed to read a file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			TGLError("Invalid format of a file %s", m_bfile.file_name().c_str());
		}

		if (isinf(m_last_val))
			m_last_val = numeric_limits<float>::quiet_NaN();
		return true;
	}

	void flush() {
		if (m_last_interval.start != -1)
			write_last_interval();
	}

	int seek(long offset, int whence) { return m_bfile.seek(offset, whence); }

	const string &file_name() { return m_bfile.file_name(); }

	bool opened() const { return m_bfile.opened(); }

	int chromid() const { return m_last_interval.chromid; }

	const GInterval &last_interval() const { return m_last_interval; }
	float last_val() const { return m_last_val; }

private:
	BufferedFile m_bfile;
	GInterval    m_last_interval;
	float        m_last_val;

	void write_last_interval() {
		static const int RECORD_SIZE = sizeof(m_last_interval.start) + sizeof(m_last_interval.end) + sizeof(m_last_val);

		uint64_t size = 0;
		size += m_bfile.write(&m_last_interval.start, sizeof(m_last_interval.start));
		size += m_bfile.write(&m_last_interval.end, sizeof(m_last_interval.end));
		size += m_bfile.write(&m_last_val, sizeof(m_last_val));

		if ((int)size != RECORD_SIZE) {
			if (m_bfile.error())
				TGLError("Failed to write intervals to file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			TGLError("Failed to write intervals to file %s", m_bfile.file_name().c_str());
		}

		m_last_interval.start = -1;
	}
};

//----------------------------------------- BufferedIntervals2D ---------------------------------------------

class BufferedIntervals2D {
public:
	BufferedIntervals2D() : m_last_val(numeric_limits<float>::quiet_NaN()) {}
	BufferedIntervals2D(const BufferedIntervals2D &) : m_last_val(numeric_limits<float>::quiet_NaN()) {}

	~BufferedIntervals2D() { close(); }

	void open(const char *path, const char *mode, int chromid1, int chromid2) {
		if (m_bfile.open(path, mode))
			TGLError("Opening file %s: %s", path, strerror(errno));
		m_last_interval.chromid1() = chromid1;
		m_last_interval.chromid2() = chromid2;
	}

	void close() {
		m_bfile.close();
	}

	void write_interval(const GInterval &interval1, const GInterval &interval2, float val) {
		static const int RECORD_SIZE = sizeof(interval1.start) + sizeof(interval2.start) + sizeof(interval1.end) + sizeof(interval2.end) + sizeof(val);

		uint64_t size = 0;
		size += m_bfile.write(&interval1.start, sizeof(interval1.start));
		size += m_bfile.write(&interval1.end, sizeof(interval1.end));
		size += m_bfile.write(&interval2.start, sizeof(interval2.start));
		size += m_bfile.write(&interval2.end, sizeof(interval2.end));
		size += m_bfile.write(&val, sizeof(val));

		if ((int)size != RECORD_SIZE) {
			if (m_bfile.error())
				TGLError("Failed to write intervals to file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			TGLError("Failed to write intervals to file %s", m_bfile.file_name().c_str());
		}
	}

	bool read_interval() {
		if (m_bfile.read(&m_last_interval.start1(), sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_interval.end1(), sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_interval.start2(), sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_interval.end2(), sizeof(int64_t)) != sizeof(int64_t) ||
				m_bfile.read(&m_last_val, sizeof(float)) != sizeof(float))
		{
			if (m_bfile.eof())
				return false;
			if (m_bfile.error())
				TGLError("Failed to read a file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
			TGLError("Invalid format of a file %s", m_bfile.file_name().c_str());
		}

		if (isinf(m_last_val))
			m_last_val = numeric_limits<float>::quiet_NaN();
		return true;
	}

	void flush() {}

	int seek(long offset, int whence) { return m_bfile.seek(offset, whence); }

	const string &file_name() { return m_bfile.file_name(); }

	bool opened() const { return m_bfile.opened(); }

	BufferedFile &bfile() { return m_bfile; }

	int chromid1() const { return m_last_interval.chromid1(); }
	int chromid2() const { return m_last_interval.chromid2(); }

	const GInterval2D &last_interval() const { return m_last_interval; }
	float last_val() const { return m_last_val; }

private:
	BufferedFile m_bfile;
	GInterval2D  m_last_interval;
	float        m_last_val;
};

//----------------------------------------- BufferedIntervalsCache ---------------------------------------------

// Cache manager for BufferedIntervals to avoid having too many open files
class BufferedIntervalsCache {
public:
	BufferedIntervalsCache(const string &dir, const GenomeChromKey &chromkey, size_t max_open = 256)
		: m_dir(dir), m_chromkey(chromkey), m_max_open(max_open) {}

	~BufferedIntervalsCache() {
		// Close all open files
		for (map<int, BufferedIntervals*>::iterator it = m_cache.begin(); it != m_cache.end(); ++it) {
			if (it->second) {
				it->second->close();
				delete it->second;
			}
		}
	}

	BufferedIntervals* get(int chromid) {
		// Check if already in cache
		map<int, BufferedIntervals*>::iterator it = m_cache.find(chromid);
		if (it != m_cache.end()) {
			// Move to end of LRU list
			m_lru.remove(chromid);
			m_lru.push_back(chromid);
			return it->second;
		}

		// Need to open new file
		// First check if cache is full
		if (m_cache.size() >= m_max_open) {
			// Remove least recently used
			int lru_chromid = m_lru.front();
			m_lru.pop_front();

			BufferedIntervals *lru = m_cache[lru_chromid];
			lru->close();
			delete lru;
			m_cache.erase(lru_chromid);
		}

		// Open new file
		// Check if this chromid was opened before (use "a+" for append, "w+" for new)
		bool previously_opened = (m_opened_chroms.find(chromid) != m_opened_chroms.end());
		const char *mode = previously_opened ? "a+" : "w+";

		BufferedIntervals *bi = new BufferedIntervals();
		char filename[FILENAME_MAX];
		snprintf(filename, sizeof(filename), "%s/_%s", m_dir.c_str(), m_chromkey.id2chrom(chromid).c_str());
		bi->open(filename, mode, chromid);

		m_cache[chromid] = bi;
		m_lru.push_back(chromid);
		m_opened_chroms.insert(chromid);
		return bi;
	}

	void flush_all() {
		for (map<int, BufferedIntervals*>::iterator it = m_cache.begin(); it != m_cache.end(); ++it) {
			if (it->second) {
				it->second->flush();
			}
		}
	}

	void close_all() {
		for (map<int, BufferedIntervals*>::iterator it = m_cache.begin(); it != m_cache.end(); ++it) {
			if (it->second) {
				it->second->close();
				delete it->second;
			}
		}
		m_cache.clear();
		m_lru.clear();
	}

	// Get a BufferedIntervals for reading (not from cache)
	BufferedIntervals* get_for_reading(int chromid) {
		BufferedIntervals *bi = new BufferedIntervals();
		char filename[FILENAME_MAX];
		snprintf(filename, sizeof(filename), "%s/_%s", m_dir.c_str(), m_chromkey.id2chrom(chromid).c_str());
		bi->open(filename, "r", chromid);
		return bi;
	}

private:
	string m_dir;
	const GenomeChromKey &m_chromkey;
	size_t m_max_open;
	map<int, BufferedIntervals*> m_cache;
	list<int> m_lru;  // Least recently used list
	set<int> m_opened_chroms;  // Track which chroms have been opened before
};

struct GIntervalVal {
	GInterval interval;
	float     val;
	int64_t   chain_id;

	GIntervalVal(const GInterval &_interval, float _val, int64_t _chain_id) :
		interval(_interval), val(_val), chain_id(_chain_id) {}
	bool operator<(const GIntervalVal &obj) const {
		if (interval.start != obj.interval.start)
			return interval.start < obj.interval.start;
		if (interval.end != obj.interval.end)
			return interval.end < obj.interval.end;
		if (chain_id != obj.chain_id)
			return chain_id < obj.chain_id;
		return val < obj.val;
	}
};

// Aggregate possibly-overlapping mapped rectangles into DISJOINT rectangles before
// inserting them into a StatQuadTree (which requires non-overlapping objects -
// overlapping inserts corrupt read-back and double-count). Disjoint source rects can
// map onto overlapping target rects because the chain shifts x and y independently, so
// the 2D path needs the same collect -> segment -> aggregate treatment the 1D path has.
//
// We coordinate-compress the x and y boundaries into a grid and, for every grid cell,
// aggregate (via agg_cfg) the values of all source rects covering it. Each cell carries
// a unique contribution id so aggregate_values never folds distinct sources together.
// Cells are emitted one rect each (no run merging).
//
// ponytail: per-cell cost is O(active rects per x-slab); for grid-aligned data
// (uniform bins, each rect == one cell) that is ~O(N), but it is O(N^2)
// worst case for pathological nested/offset rectangles, and a non-grid-aligned disjoint
// rect gets split at a neighbour's boundary (more rects, identical values). Upgrade path
// if it ever bites: decompose per overlap-cluster and merge equal-value runs.
struct RectVal2D { int64_t x1, y1, x2, y2; float v; };

static void aggregate_2d_rects(vector<RectVal2D> &rects, const AggregationConfig &agg_cfg, RectsQuadTree &qtree)
{
	if (rects.empty())
		return;

	vector<int64_t> xs, ys;
	xs.reserve(rects.size() * 2);
	ys.reserve(rects.size() * 2);
	for (vector<RectVal2D>::const_iterator r = rects.begin(); r != rects.end(); ++r) {
		xs.push_back(r->x1); xs.push_back(r->x2);
		ys.push_back(r->y1); ys.push_back(r->y2);
	}
	sort(xs.begin(), xs.end()); xs.erase(unique(xs.begin(), xs.end()), xs.end());
	sort(ys.begin(), ys.end()); ys.erase(unique(ys.begin(), ys.end()), ys.end());

	vector<size_t> by_x1(rects.size()), by_x2(rects.size());
	for (size_t i = 0; i < rects.size(); ++i) { by_x1[i] = i; by_x2[i] = i; }
	sort(by_x1.begin(), by_x1.end(), [&](size_t a, size_t b) { return rects[a].x1 < rects[b].x1; });
	sort(by_x2.begin(), by_x2.end(), [&](size_t a, size_t b) { return rects[a].x2 < rects[b].x2; });

	set<size_t> active;
	size_t ia = 0, ir = 0;
	int64_t encounter = 0; // unique id per contribution -> aggregate_values never merges them

	for (size_t xi = 0; xi + 1 < xs.size(); ++xi) {
		int64_t xa = xs[xi], xb = xs[xi + 1];

		while (ir < by_x2.size() && rects[by_x2[ir]].x2 <= xa) { active.erase(by_x2[ir]); ++ir; }
		while (ia < by_x1.size() && rects[by_x1[ia]].x1 <= xa) { active.insert(by_x1[ia]); ++ia; }

		if (active.empty())
			continue;

		// Per y-band aggregation over the rects active in this x-slab. Every active rect
		// fully covers each cell it touches, so all contributors weigh equally (len = 1).
		map<int, AggregationState> band_states;
		for (set<size_t>::const_iterator it = active.begin(); it != active.end(); ++it) {
			const RectVal2D &r = rects[*it];
			int b_lo = (int)(lower_bound(ys.begin(), ys.end(), r.y1) - ys.begin());
			int b_hi = (int)(lower_bound(ys.begin(), ys.end(), r.y2) - ys.begin());
			for (int b = b_lo; b < b_hi; ++b) {
				aggregation_state_add(band_states[b], (double)r.v, 1.0, 1.0, encounter, encounter, encounter);
				++encounter;
			}
		}

		for (map<int, AggregationState>::iterator kv = band_states.begin(); kv != band_states.end(); ++kv) {
			double aggregated = aggregate_values(agg_cfg, kv->second);
			float out_val = numeric_limits<float>::quiet_NaN();
			if (!std::isnan(aggregated))
				out_val = (float)aggregated;
			qtree.insert(RectsQuadTree::ValueType(xa, ys[kv->first], xb, ys[kv->first + 1], out_val));
		}
		check_interrupt();
	}
}

// Lifted points of a POINTS track are buffered as (x, y, val) records. In memory a point
// also carries seq, its position in the buffer, which is its contribution id: it orders
// first / last / nth like the encounter id above and keeps aggregate_values from merging
// distinct points.
struct PointVal2D { int64_t x, y; float v; uint32_t seq; };

static const int POINT_RECORD_SIZE = 2 * sizeof(int64_t) + sizeof(float);

static void write_point(BufferedFile &file, int64_t x, int64_t y, float v)
{
	char buf[POINT_RECORD_SIZE];
	memcpy(buf, &x, sizeof(x));
	memcpy(buf + sizeof(x), &y, sizeof(y));
	memcpy(buf + sizeof(x) + sizeof(y), &v, sizeof(v));
	if (file.write(buf, POINT_RECORD_SIZE) != POINT_RECORD_SIZE) {
		if (file.error())
			TGLError("Failed to write points to file %s: %s", file.file_name().c_str(), strerror(errno));
		TGLError("Failed to write points to file %s", file.file_name().c_str());
	}
}

static bool read_point(BufferedFile &file, int64_t &x, int64_t &y, float &v)
{
	char buf[POINT_RECORD_SIZE];
	uint64_t size = file.read(buf, POINT_RECORD_SIZE);
	if (size != POINT_RECORD_SIZE) {
		if (!size && file.eof())
			return false;
		if (file.error())
			TGLError("Failed to read a file %s: %s", file.file_name().c_str(), strerror(errno));
		TGLError("Invalid format of a file %s", file.file_name().c_str());
	}
	memcpy(&x, buf, sizeof(x));
	memcpy(&y, buf + sizeof(x), sizeof(y));
	memcpy(&v, buf + sizeof(x) + sizeof(y), sizeof(v));

	// an infinite source value becomes NaN, as in BufferedIntervals2D::read_interval
	if (isinf(v))
		v = numeric_limits<float>::quiet_NaN();
	return true;
}

// Writes the points lifted into one target chromosome pair (buffered in the file infilename)
// as a POINTS track file. Points that landed on the same target point are merged by agg_cfg.
// As in gtrack.2d.import_contacts the pair is split into subtrees by gmax.data.size, and only
// one subtree is held in memory at a time. The buffer file is removed once it is read: after
// the split, so that the pair is not on disk twice.
static void write_2d_points(const string &infilename, const char *filename, const string &dirname, int chromid1, int chromid2,
                            IntervUtils &iu, const AggregationConfig &agg_cfg)
{
	BufferedFile infile;
	if (infile.open(infilename.c_str(), "r"))
		TGLError("Opening file %s: %s", infilename.c_str(), strerror(errno));

	int64_t num_points = infile.file_size() / POINT_RECORD_SIZE;
	// Unsigned: an unset gmax.data.size reads as UINT64_MAX. Below 1 (0, or 0.5 read as 0)
	// it counts as 1 rather than dividing by zero.
	uint64_t max_data_size = max(iu.get_max_data_size(), (uint64_t)1);
	int64_t num_subtrees = max((int64_t)((uint64_t)num_points / max_data_size), (int64_t)1);
	num_subtrees = (int64_t)1 << (2 * (int)(log2(num_subtrees) / 2));  // round down to a power of 4

	GenomeTrackRectsPoints gtrack(iu.get_track_chunk_size(), iu.get_track_num_chunks());
	PointsQuadTreeCachedSerializer serializer;
	gtrack.init_write(filename, chromid1, chromid2);
	gtrack.init_serializer(serializer, 0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2), num_subtrees, false);

	BufferedFiles subtree_files(num_subtrees > 1 ? num_subtrees : 0);
	int64_t x, y;
	float v;

	if (num_subtrees > 1) {
		const Rectangles &subarenas = serializer.get_subarenas();
		char subtree_filename[FILENAME_MAX];

		for (int64_t i = 0; i < num_subtrees; ++i) {
			snprintf(subtree_filename, sizeof(subtree_filename), "%s/_subtree%ld", dirname.c_str(), (long)i);
			subtree_files[i] = new BufferedFile();
			if (subtree_files[i]->open(subtree_filename, "w+"))
				TGLError("Opening file %s: %s", subtree_filename, strerror(errno));
		}

		while (read_point(infile, x, y, v)) {
			Point point(x, y);
			for (Rectangles::const_iterator isubarena = subarenas.begin(); isubarena != subarenas.end(); ++isubarena) {
				if (point.do_intersect(*isubarena)) {
					write_point(*subtree_files[isubarena - subarenas.begin()], x, y, v);
					break;
				}
			}
			check_interrupt();
		}

		for (int64_t i = 0; i < num_subtrees; ++i) {
			if (subtree_files[i]->flush())
				TGLError("Failed to write points to file %s: %s", subtree_files[i]->file_name().c_str(), strerror(errno));
			subtree_files[i]->seek(0, SEEK_SET);
		}

		infile.close();
		unlink(infilename.c_str());
	}

	vector<PointVal2D> points;
	AggregationState agg_state;

	for (int64_t i = 0; i < num_subtrees; ++i) {
		BufferedFile &file = num_subtrees > 1 ? *subtree_files[i] : infile;

		// seq is 32-bit; check the count before reserving, or a huge subtree fails as "Out of memory"
		int64_t num_file_points = file.file_size() / POINT_RECORD_SIZE;
		if (num_file_points > (int64_t)numeric_limits<uint32_t>::max())
			TGLError("Too many points lifted into chromosome pair (%s, %s); lower gmax.data.size",
					 iu.id2chrom(chromid1).c_str(), iu.id2chrom(chromid2).c_str());

		points.clear();
		points.reserve(num_file_points);
		while (read_point(file, x, y, v)) {
			PointVal2D p = { x, y, v, (uint32_t)points.size() };
			points.push_back(p);
			check_interrupt();
		}

		string read_filename = file.file_name();
		file.close();
		unlink(read_filename.c_str());

		sort(points.begin(), points.end(), [](const PointVal2D &a, const PointVal2D &b) {
			if (a.x != b.x)
				return a.x < b.x;
			if (a.y != b.y)
				return a.y < b.y;
			return a.seq < b.seq;
		});

		for (size_t first = 0, last; first < points.size(); first = last) {
			agg_state.reset();
			for (last = first; last < points.size() && points[last].x == points[first].x && points[last].y == points[first].y; ++last)
				aggregation_state_add(agg_state, (double)points[last].v, 1.0, 1.0, points[last].seq, points[last].seq, points[last].seq);

			double aggregated = aggregate_values(agg_cfg, agg_state);
			float out_val = numeric_limits<float>::quiet_NaN();
			if (!std::isnan(aggregated))
				out_val = (float)aggregated;
			serializer.insert(PointsQuadTree::ValueType(points[first].x, points[first].y, out_val));
			check_interrupt();
		}
	}

	serializer.end();
}

// The chrom id by which the source track is read, per chain source chromosome (-1: the source
// genome lacks it). The files of a per-chromosome track are found by name, so the chain's ids
// serve. An indexed track (track.idx) is keyed by the chrom ids of its own database, which
// gtrack.liftover passes as src_chroms (that database's chromosomes in chrom id order); a
// chain chromosome is matched there by name, as is, then without or with the "chr" prefix.
static vector<int> map_src_chainid2genomeid(const vector<string> &src_id2chrom, bool src_indexed, SEXP _src_chroms, const char *src_track_dir)
{
	vector<int> src_chainid2genomeid(src_id2chrom.size(), -1);

	if (!src_indexed) {
		for (size_t i = 0; i < src_id2chrom.size(); ++i)
			src_chainid2genomeid[i] = (int)i;
		return src_chainid2genomeid;
	}

	if (!Rf_isString(_src_chroms))
		TGLError("Source track %s is in indexed format, but no database holds it (a directory with chrom_sizes.txt whose "
				 "\"tracks\" directory contains the track); an indexed track is read by the chromosome ids of its database",
				 src_track_dir);

	map<string, int> chrom2id;
	for (int i = 0; i < Rf_length(_src_chroms); ++i)
		chrom2id.insert(make_pair(string(CHAR(STRING_ELT(_src_chroms, i))), i));

	for (size_t i = 0; i < src_id2chrom.size(); ++i) {
		const string &name = src_id2chrom[i];
		map<string, int>::const_iterator ichrom = chrom2id.find(name);
		if (ichrom == chrom2id.end())
			ichrom = name.compare(0, 3, "chr") ? chrom2id.find("chr" + name) : chrom2id.find(name.substr(3));
		if (ichrom != chrom2id.end())
			src_chainid2genomeid[i] = ichrom->second;
	}
	return src_chainid2genomeid;
}

extern "C" {

SEXP gtrack_liftover(SEXP _track,
                     SEXP _src_track_dir,
                     SEXP _src_chroms,
                     SEXP _chain,
                     SEXP _src_overlap_policy,
                     SEXP _tgt_overlap_policy,
                     SEXP _multi_target_agg,
                     SEXP _multi_target_params,
                     SEXP _na_rm,
                     SEXP _min_n,
                     SEXP _min_score,
                     SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_src_track_dir) || Rf_length(_src_track_dir) != 1)
			verror("Track source directory argument is not a string");

		if (!Rf_isNull(_src_chroms) && !Rf_isString(_src_chroms))
			verror("Source chromosomes argument is not a character vector");

		if (!Rf_isString(_src_overlap_policy) || Rf_length(_src_overlap_policy) != 1)
			verror("Source overlap policy argument is not a string");

		if (!Rf_isString(_tgt_overlap_policy) || Rf_length(_tgt_overlap_policy) != 1)
			verror("Target overlap policy argument is not a string");

		// Validate min_score (optional, currently unused)
		if (!Rf_isNull(_min_score)) {
			if (!Rf_isReal(_min_score) || Rf_length(_min_score) != 1)
				verror("min_score must be a single numeric value");
		}

		if (!Rf_isString(_multi_target_agg) || Rf_length(_multi_target_agg) != 1)
			verror("multi_target_agg argument is not a string");

		if (!Rf_isInteger(_multi_target_params) || Rf_length(_multi_target_params) != 1)
			verror("params argument must be an integer scalar");

		if (!Rf_isLogical(_na_rm) || Rf_length(_na_rm) != 1)
			verror("na_rm argument must be a logical scalar");

		if (!Rf_isInteger(_min_n) || Rf_length(_min_n) != 1)
			verror("min_n argument must be an integer scalar");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *src_track_dir = CHAR(STRING_ELT(_src_track_dir, 0));
		const char *src_overlap_policy = CHAR(STRING_ELT(_src_overlap_policy, 0));
		const char *tgt_overlap_policy = CHAR(STRING_ELT(_tgt_overlap_policy, 0));
		// Convert "auto" to "auto_score" alias
		std::string effective_tgt_policy = tgt_overlap_policy;
		if (!strcmp(tgt_overlap_policy, "auto"))
			effective_tgt_policy = "auto_score";
		const char *multi_target_agg_str = CHAR(STRING_ELT(_multi_target_agg, 0));

		int na_rm_int = Rf_asLogical(_na_rm);
		if (na_rm_int == NA_LOGICAL)
			verror("na_rm must not be NA");

		int min_n_int = Rf_asInteger(_min_n);
		if (min_n_int == NA_INTEGER)
			min_n_int = -1;
		else if (min_n_int < 0)
			verror("min_n must be non-negative");

		int nth_index = Rf_asInteger(_multi_target_params);

		AggregationConfig agg_cfg;
		agg_cfg.type = parse_aggregation_type(multi_target_agg_str);
		agg_cfg.na_rm = (na_rm_int != 0);
		agg_cfg.min_n = min_n_int;
		agg_cfg.nth_index = -1;

		if (agg_cfg.type == AggregationType::NTH) {
			if (nth_index == NA_INTEGER || nth_index <= 0)
				verror("params must be a positive integer for 'nth' aggregation");
			agg_cfg.nth_index = nth_index;
		}

		IntervUtils iu(_envir);
		ChainIntervals chain_intervs;
		vector<string> src_id2chrom;

		iu.convert_rchain_intervs(_chain, chain_intervs, src_id2chrom);

		// Handle target overlaps first
		chain_intervs.sort_by_tgt();
		chain_intervs.handle_tgt_overlaps(effective_tgt_policy, iu.get_chromkey(), src_id2chrom);

		// Handle source overlaps
		chain_intervs.sort_by_src();
		chain_intervs.handle_src_overlaps(src_overlap_policy, iu.get_chromkey(), src_id2chrom);

		// Build auxiliary structures for efficient source interval mapping
		chain_intervs.buildSrcAux();

		// Set the target overlap policy for score-based selection during liftover
		chain_intervs.set_tgt_overlap_policy(effective_tgt_policy);

		GenomeChromKey src_chromkey;
		for (vector<string>::const_iterator ichrom = src_id2chrom.begin(); ichrom != src_id2chrom.end(); ++ichrom)
			src_chromkey.add_chrom(*ichrom, numeric_limits<int64_t>::max());

		struct stat idx_st;
		bool src_indexed = stat((string(src_track_dir) + "/track.idx").c_str(), &idx_st) == 0;
		vector<int> src_chainid2genomeid = map_src_chainid2genomeid(src_id2chrom, src_indexed, _src_chroms, src_track_dir);

		string dirname = create_track_dir(_envir, track);

		GenomeTrack::Type src_track_type = GenomeTrack::get_type(src_track_dir, src_chromkey);

		if (GenomeTrack::is_1d(src_track_type)) {
			GIntervals all_genome_intervs;
			iu.get_all_genome_intervs(all_genome_intervs);
			// Collect intervals in memory instead of writing to temp files immediately
			map<int, vector<GIntervalVal> > chrom_intervals;
			char filename[FILENAME_MAX];
			GIntervals tgt_intervals;
			vector<ChainMappingMetadata> mapping_meta;
			unsigned binsize = 0;

			Progress_reporter progress;
			progress.init(src_id2chrom.size() + iu.get_chromkey().get_num_chroms(), 1);

			// write the target intervals + values to the temporary files
			if (src_track_type == GenomeTrack::FIXED_BIN) {
				GenomeTrackFixedBin src_track;
				for (vector<string>::const_iterator ichrom = src_id2chrom.begin(); ichrom != src_id2chrom.end(); ++ichrom) {
					int src_chromid_in_chain = ichrom - src_id2chrom.begin();  // chromid in the chain's coordinate system
					int chromid_to_use = src_chainid2genomeid[src_chromid_in_chain];
					float val;

					if (chromid_to_use < 0) {  // not in the source genome
						progress.report(1);
						continue;
					}

					try {
						snprintf(filename, sizeof(filename), "%s/%s", src_track_dir, ichrom->c_str());
						src_track.init_read(filename, chromid_to_use);
						if (binsize > 0 && binsize != src_track.get_bin_size()) {
							char filename2[FILENAME_MAX];
							snprintf(filename2, sizeof(filename2), "%s/%s", src_track_dir, (ichrom - 1)->c_str());
							TGLError("Binsize of track file %s differs from the binsize of track file %s (%d vs. %d)",
									filename, filename2, src_track.get_bin_size(), binsize);
						} else
							binsize = src_track.get_bin_size();
					} catch (TGLException &) {  // some of source chroms might be missing, this is normal
						progress.report(1);
						continue;
					}

					GInterval src_interval(src_chromid_in_chain, 0, src_track.get_bin_size(), 0);
					ChainIntervals::const_iterator hint = chain_intervs.begin();

					for (int64_t i = 0; i < src_track.get_num_samples(); ++i) {
						src_track.read_next_bin(val);

						mapping_meta.clear();
						hint = chain_intervs.map_interval(src_interval, tgt_intervals, hint, &mapping_meta);
						if (!tgt_intervals.empty() && mapping_meta.size() != tgt_intervals.size())
							TGLError("Metadata size mismatch: %zu vs %zu", mapping_meta.size(), tgt_intervals.size());
						for (size_t idx = 0; idx < tgt_intervals.size(); ++idx) {
							// Aggregation dedup key = (chain id, source bin index i). Pieces of
							// the SAME source bin split by "agg" target segmentation share the key
							// (counted once), while DIFFERENT source bins of the same chain that
							// land in one target bin stay distinct and are aggregated. Keying by
							// chain_id alone collapsed them, returning only the first bin's value.
							// (Mirrors IntervalsLiftover's (chain_id, interval_id) key.)
							int64_t agg_key = ((int64_t)mapping_meta[idx].chain_id << 32) ^ (int64_t)i;
							chrom_intervals[tgt_intervals[idx].chromid].push_back(GIntervalVal(tgt_intervals[idx], val, agg_key));
						}

						src_interval.start += src_track.get_bin_size();
						src_interval.end += src_track.get_bin_size();
						check_interrupt();
					}

					progress.report(1);
				}
			} else if (src_track_type == GenomeTrack::SPARSE) {
				GenomeTrackSparse src_track;
				for (vector<string>::const_iterator ichrom = src_id2chrom.begin(); ichrom != src_id2chrom.end(); ++ichrom) {
					int src_chromid_in_chain = ichrom - src_id2chrom.begin();  // chromid in the chain's coordinate system
					int chromid_to_use = src_chainid2genomeid[src_chromid_in_chain];

					if (chromid_to_use < 0) {  // not in the source genome
						progress.report(1);
						continue;
					}

					try {
						snprintf(filename, sizeof(filename), "%s/%s", src_track_dir, ichrom->c_str());
						src_track.init_read(filename, chromid_to_use);
					} catch (TGLException &) {  // some of source chroms might be missing, this is normal
						progress.report(1);
						continue;
					}

					const GIntervals &src_intervals = src_track.get_intervals();
					const vector<float> &vals = src_track.get_vals();

					// The intervals from the track have chromid=src_chromid_in_genome,
					// but the chain expects chromid=src_chromid_in_chain
					for (uint64_t i = 0; i < src_intervals.size(); ++i) {
						GInterval remapped_interval = src_intervals[i];
						remapped_interval.chromid = src_chromid_in_chain;
						// Reset hint for each interval to avoid missing overlaps with non-consecutive chains
						ChainIntervals::const_iterator hint = chain_intervs.begin();
						mapping_meta.clear();
						hint = chain_intervs.map_interval(remapped_interval, tgt_intervals, hint, &mapping_meta);
						for (size_t idx = 0; idx < tgt_intervals.size(); ++idx) {
							// Dedup key = (chain id, source interval index i); see the FIXED_BIN
							// path above. Keeps distinct source intervals of one chain separate.
							int64_t agg_key = ((int64_t)mapping_meta[idx].chain_id << 32) ^ (int64_t)i;
							chrom_intervals[tgt_intervals[idx].chromid].push_back(GIntervalVal(tgt_intervals[idx], vals[i], agg_key));
						}
						check_interrupt();
					}

					progress.report(1);
				}
			} else
				TGLError("Source track type %s is currently not supported in liftover", GenomeTrack::TYPE_NAMES[src_track_type]);

			// Process collected intervals for each chromosome, sort and save to track files.
			//
			// On indexed-format target DBs (large-contig genomes), stream the per-chrom
			// bytes directly into track.dat + track.idx via GenomeTrackIndexedWriter
			// instead of writing one file per chrom. The legacy per-chrom path created
			// N files even for chroms not covered by the chain (mostly empty
			// signature-only sparse files, or NaN-filled fixed-bin files); on a 1M-contig
			// target genome that is ~1M open+close syscalls. Byte layout is bit-for-bit
			// identical to the convert-after-create path (see test-track-indexed-direct.R
			// for the sparse and fixed-bin equivalence tests).
			//
			// The chrom loop already iterates 0..num_chroms in strictly increasing order
			// which satisfies GenomeTrackIndexedWriter::begin_chrom's monotonicity check.
			const bool indexed = rdb::is_db_indexed(_envir);
			GenomeTrackIndexedWriter iw;
			vector<char> sparse_header_bytes;
			vector<char> fixedbin_header_bytes;
			if (indexed) {
				GenomeTrack::Type out_type = (src_track_type == GenomeTrack::FIXED_BIN) ? GenomeTrack::FIXED_BIN : GenomeTrack::SPARSE;
				iw.init(dirname, out_type, (uint32_t)iu.get_chromkey().get_num_chroms());
				if (out_type == GenomeTrack::SPARSE)
					GenomeTrackSparse::pack_header(sparse_header_bytes);
				else if (binsize > 0)
					GenomeTrackFixedBin::pack_header(fixedbin_header_bytes, binsize);
			}

			for (int chromid = 0; chromid < (int)iu.get_chromkey().get_num_chroms(); ++chromid) {
				// Create empty file if this chromosome has no data
				if (chrom_intervals.find(chromid) == chrom_intervals.end()) {
					// Always create empty chromosome files so downstream readers and indexing
					// see the expected number of bins/intervals, even when the database is indexed.
					if (src_track_type == GenomeTrack::FIXED_BIN) {
						// Only create empty fixed bin tracks if we know the binsize (i.e., at least one chromosome was processed)
						if (binsize > 0) {
							int64_t chrom_size = iu.get_chromkey().get_chrom_size(chromid);
							int64_t end_bin = (int64_t)ceil(chrom_size / (double)binsize);
							if (indexed) {
								iw.begin_chrom(chromid);
								iw.write_bytes(fixedbin_header_bytes.data(), fixedbin_header_bytes.size());
								if (end_bin > 0) {
									const int64_t chunk_size = 65536;
									vector<float> na_chunk((size_t)min<int64_t>(end_bin, chunk_size), numeric_limits<float>::quiet_NaN());
									int64_t remaining = end_bin;
									while (remaining > 0) {
										uint64_t to_write = (uint64_t)min<int64_t>(remaining, (int64_t)na_chunk.size());
										iw.write_bytes(&na_chunk[0], to_write * sizeof(float));
										remaining -= to_write;
									}
								}
								iw.end_chrom();
							} else {
								snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());
								GenomeTrackFixedBin gtrack;
								gtrack.init_write(filename, binsize, chromid);
								// Fill the chromosome with NaN values so the bin count matches the chromosome size
								if (end_bin > 0) {
									const int64_t chunk_size = 65536;
									vector<float> na_chunk((size_t)min<int64_t>(end_bin, chunk_size), numeric_limits<float>::quiet_NaN());
									int64_t remaining = end_bin;
									while (remaining > 0) {
										uint64_t to_write = (uint64_t)min<int64_t>(remaining, (int64_t)na_chunk.size());
										gtrack.write_next_bins(&na_chunk[0], to_write);
										remaining -= to_write;
									}
								}
							}
						}
						// indexed && binsize==0: skip entirely (no source data at all,
						// matches legacy per-chrom behaviour of writing no file).
					} else if (src_track_type == GenomeTrack::SPARSE) {
						if (indexed) {
							iw.begin_chrom(chromid);
							iw.write_bytes(sparse_header_bytes.data(), sparse_header_bytes.size());
							iw.end_chrom();
						} else {
							snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());
							GenomeTrackSparse gtrack;
							gtrack.init_write(filename, chromid);
						}
					}
					progress.report(1);
					continue;
				}

				vector<GIntervalVal> &interv_vals = chrom_intervals[chromid];

				sort(interv_vals.begin(), interv_vals.end());

				if (src_track_type == GenomeTrack::FIXED_BIN) {
					GenomeTrackFixedBin gtrack;
					if (indexed) {
						iw.begin_chrom(chromid);
						iw.write_bytes(fixedbin_header_bytes.data(), fixedbin_header_bytes.size());
					} else {
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());
						gtrack.init_write(filename, binsize, chromid);
					}
					int64_t chrom_size = iu.get_chromkey().get_chrom_size(chromid);
					int64_t end_bin = (int64_t)ceil(chrom_size / (double)binsize);
					int64_t coord1 = 0;
					int64_t coord2 = coord1 + binsize;
					vector<GIntervalVal>::const_iterator iinterv_val = interv_vals.begin();
					AggregationState agg_state;
					agg_state.contributions.reserve(8);

					for (int64_t bin = 0; bin < end_bin; ++bin) {
						agg_state.reset();

						// Advance the persistent cursor past contributions that end at or before this
						// bin's start. interv_vals is sorted by start, so a contribution with end <=
						// coord1 cannot overlap this bin or any later one (coord1 only increases).
						// Contributions with end > coord1 may still overlap a later bin, so we must NOT
						// advance past them here: doing so dropped a contribution that spanned a bin
						// boundary whenever an overlapping sibling (e.g. a parallel chain kept by the
						// "agg" target-overlap policy) shared its interval.
						while (iinterv_val != interv_vals.end() && iinterv_val->interval.end <= coord1)
							++iinterv_val;

						// Collect every contribution overlapping [coord1, coord2). Sorted by start, so
						// stop as soon as a contribution starts at or after coord2.
						for (vector<GIntervalVal>::const_iterator iter = iinterv_val; iter != interv_vals.end(); ++iter) {
							if (iter->interval.start >= coord2)
								break;
							int64_t overlap_start = max(coord1, iter->interval.start);
							int64_t overlap_end = min(coord2, iter->interval.end);
							if (overlap_start < overlap_end) {
								double overlap_len = static_cast<double>(overlap_end - overlap_start);
								int64_t bin_end_clamped = std::min<int64_t>(coord2, chrom_size);
								double locus_len = static_cast<double>(std::max<int64_t>(0, bin_end_clamped - coord1));
								if (locus_len == 0.0)
									locus_len = static_cast<double>(coord2 - coord1);
								aggregation_state_add(
									agg_state,
									static_cast<double>(iter->val),
									overlap_len,
									locus_len,
									overlap_start,
									overlap_end,
									iter->chain_id
								);
							}
							check_interrupt();
						}

						double aggregated = aggregate_values(agg_cfg, agg_state);
						float out_val = numeric_limits<float>::quiet_NaN();
						if (!std::isnan(aggregated))
							out_val = static_cast<float>(aggregated);

						if (indexed)
							iw.write_bytes(&out_val, sizeof(out_val));
						else
							gtrack.write_next_bin(out_val);

						coord1 = coord2;
						coord2 += binsize;
						check_interrupt();
					}
					if (indexed)
						iw.end_chrom();
				} else if (src_track_type == GenomeTrack::SPARSE) {
					GenomeTrackSparse gtrack;
					if (indexed) {
						iw.begin_chrom(chromid);
						iw.write_bytes(sparse_header_bytes.data(), sparse_header_bytes.size());
					} else {
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());
						gtrack.init_write(filename, chromid);
					}
					AggregationState agg_state;
					agg_state.contributions.reserve(4);
					vector<char> rec_buf;

					// Merge overlapping intervals before writing
					// Intervals are sorted by start, then end
					size_t idx = 0;
					while (idx < interv_vals.size()) {
						// Start a new merged interval
						int64_t merged_start = interv_vals[idx].interval.start;
						int64_t merged_end = interv_vals[idx].interval.end;
						agg_state.reset();

						// Collect all intervals that overlap with the current merged region
						while (idx < interv_vals.size() &&
								interv_vals[idx].interval.start < merged_end) {
							// This interval overlaps with merged region
							// Extend merged_end if needed
							if (interv_vals[idx].interval.end > merged_end)
								merged_end = interv_vals[idx].interval.end;

							double contrib_len = static_cast<double>(
								std::max<int64_t>(0, interv_vals[idx].interval.end - interv_vals[idx].interval.start));
							if (contrib_len == 0.0)
								contrib_len = 1.0;

							aggregation_state_add(
								agg_state,
								static_cast<double>(interv_vals[idx].val),
								contrib_len,
								contrib_len,
								interv_vals[idx].interval.start,
								interv_vals[idx].interval.end,
								interv_vals[idx].chain_id
							);
							++idx;
							check_interrupt();
						}

						double aggregated = aggregate_values(agg_cfg, agg_state);
						float out_val = numeric_limits<float>::quiet_NaN();
						if (!std::isnan(aggregated))
							out_val = static_cast<float>(aggregated);

						GInterval merged_interval(chromid, merged_start, merged_end, 0);
						if (indexed) {
							rec_buf.clear();
							GenomeTrackSparse::pack_record(rec_buf, merged_interval.start, merged_interval.end, out_val);
							iw.write_bytes(rec_buf.data(), rec_buf.size());
						} else {
							gtrack.write_next_interval(merged_interval, out_val);
						}
						check_interrupt();
					}
					if (indexed)
						iw.end_chrom();
				}

				progress.report(1);
			}

			if (indexed)
				iw.finalize();

			progress.report_last();
		} else {   // 2D tracks
			if (src_track_type != GenomeTrack::RECTS && src_track_type != GenomeTrack::POINTS)
				TGLError("Source track type %s is currently not supported in liftover", GenomeTrack::TYPE_NAMES[src_track_type]);

			// The buffers are keyed by target chrom-pair. Sizing this container from the
			// whole-genome 2D intervals set was wrong twice over: that set is a deferred
			// placeholder (length 0) on databases above gmulticontig.2d.threshold contigs,
			// making every buffered_intervs[chromid1 * num_chroms + chromid2] access run off
			// the end of an empty vector; and even when it is materialised it allocates
			// num_chroms^2 objects (207,025 on hg38 with scaffolds) to open a handful of
			// files. A sparse map creates a buffer only for the pairs that are written to.
			map<ChromPair, BufferedIntervals2D> buffered_intervs;
			char filename[FILENAME_MAX];
			GIntervals tgt_intervals[2];

			Progress_reporter progress;
			progress.init(src_id2chrom.size() * src_id2chrom.size(), 1);

			// The source is read in one pass, so it needs only the chunks on the iterator's stack.
			// gtrack.num.chunks = 0 (the default) caches every chunk it reads, i.e. the whole
			// chromosome-pair file (tens of GB for a Hi-C pair); bound the cache instead.
			int64_t src_num_chunks = iu.get_track_num_chunks() ? iu.get_track_num_chunks() : 100;

			// convert source intervals and write them to files
			for (vector<string>::const_iterator ichrom1 = src_id2chrom.begin(); ichrom1 != src_id2chrom.end(); ++ichrom1) {
				for (vector<string>::const_iterator ichrom2 = src_id2chrom.begin(); ichrom2 != src_id2chrom.end(); ++ichrom2) {
					GenomeTrackRectsRects src_track_rects(iu.get_track_chunk_size(), src_num_chunks);
					GenomeTrackRectsPoints src_track_points(iu.get_track_chunk_size(), src_num_chunks);
					GenomeTrack2D *src_track;

					if (src_track_type == GenomeTrack::RECTS)
						 src_track = &src_track_rects;
					else
						src_track = &src_track_points;

					int chromid1 = ichrom1 - src_id2chrom.begin();
					int chromid2 = ichrom2 - src_id2chrom.begin();

					// An indexed source is looked up by the source genome's chrom ids, not the
					// chain's (as in the 1D path); with the chain's ids it read another pair.
					int read_chromid1 = src_chainid2genomeid[chromid1];
					int read_chromid2 = src_chainid2genomeid[chromid2];
					if (read_chromid1 < 0 || read_chromid2 < 0) {  // not in the source genome
						progress.report(1);
						continue;
					}

					try {
						snprintf(filename, sizeof(filename), "%s/%s-%s", src_track_dir, ichrom1->c_str(), ichrom2->c_str());
						src_track->init_read(filename, read_chromid1, read_chromid2);
						if (!src_track->has_data_for_pair()) {
							progress.report(1);
							continue;
						}
					} catch (TGLException &) {  // some of source chroms might be missing, this is normal
						progress.report(1);
						continue;
					}

					// The quadtree iterates rects/points in spatial (not coordinate) order,
					// so a carried map_interval hint is not a monotonically advancing lower
					// bound and would drop mappings. Map each axis from begin() every time
					// (the same per-interval reset the 1D sparse path uses).
					GInterval src_intervals[2];

					src_intervals[0].chromid = chromid1;
					src_intervals[1].chromid = chromid2;

					if (src_track_type == GenomeTrack::RECTS) {
						src_track_rects.load();
						RectsQuadTreeCached &qtree = src_track_rects.get_qtree();
						RectsQuadTreeCached::Iterator iqtree(&qtree);

						for (iqtree.begin(); !iqtree.is_end(); iqtree.next()) {
							src_intervals[0].start = iqtree->x1;
							src_intervals[0].end = iqtree->x2;
							src_intervals[1].start = iqtree->y1;
							src_intervals[1].end = iqtree->y2;

							chain_intervs.map_interval(src_intervals[0], tgt_intervals[0], chain_intervs.begin());
							chain_intervs.map_interval(src_intervals[1], tgt_intervals[1], chain_intervs.begin());

							for (GIntervals::const_iterator iinterv1 = tgt_intervals[0].begin(); iinterv1 != tgt_intervals[0].end(); ++iinterv1) {
								for (GIntervals::const_iterator iinterv2 = tgt_intervals[1].begin(); iinterv2 != tgt_intervals[1].end(); ++iinterv2) {
									BufferedIntervals2D &buffered_interv = buffered_intervs[ChromPair(iinterv1->chromid, iinterv2->chromid)];
									if (!buffered_interv.opened()) {
										snprintf(filename, sizeof(filename), "%s/_%s-%s", dirname.c_str(), iu.get_chromkey().id2chrom(iinterv1->chromid).c_str(), iu.get_chromkey().id2chrom(iinterv2->chromid).c_str());
										buffered_interv.open(filename, "w+", iinterv1->chromid, iinterv2->chromid);
									}
									buffered_interv.write_interval(*iinterv1, *iinterv2, iqtree->v);
								}
							}
							check_interrupt();
						}
					} else {
						src_track_points.load();
						PointsQuadTreeCached &qtree = src_track_points.get_qtree();
						PointsQuadTreeCached::Iterator iqtree(&qtree);

						for (iqtree.begin(); !iqtree.is_end(); iqtree.next()) {
							src_intervals[0].start = iqtree->x;
							src_intervals[0].end = iqtree->x + 1;
							src_intervals[1].start = iqtree->y;
							src_intervals[1].end = iqtree->y + 1;

							chain_intervs.map_interval(src_intervals[0], tgt_intervals[0], chain_intervs.begin());
							chain_intervs.map_interval(src_intervals[1], tgt_intervals[1], chain_intervs.begin());

							for (GIntervals::const_iterator iinterv1 = tgt_intervals[0].begin(); iinterv1 != tgt_intervals[0].end(); ++iinterv1) {
								for (GIntervals::const_iterator iinterv2 = tgt_intervals[1].begin(); iinterv2 != tgt_intervals[1].end(); ++iinterv2) {
									BufferedIntervals2D &buffered_interv = buffered_intervs[ChromPair(iinterv1->chromid, iinterv2->chromid)];
									if (!buffered_interv.opened()) {
										snprintf(filename, sizeof(filename), "%s/_%s-%s", dirname.c_str(), iu.get_chromkey().id2chrom(iinterv1->chromid).c_str(), iu.get_chromkey().id2chrom(iinterv2->chromid).c_str());
										buffered_interv.open(filename, "w+", iinterv1->chromid, iinterv2->chromid);
									}
									write_point(buffered_interv.bfile(), iinterv1->start, iinterv2->start, iqtree->v);
								}
							}
							check_interrupt();
						}
					}

					progress.report(1);
				}
			}

			progress.report_last();

			// Close every pair buffer before the pairs are written, and reopen each when its pair
			// is written. A pair split into subtrees opens a file per subtree; with all the pair
			// buffers (and their stdio buffers) still open, hundreds of target pairs could exceed
			// the limit on open files. gtrack.2d.import_contacts closes its pair files the same way.
			for (map<ChromPair, BufferedIntervals2D>::iterator ibuffered_interv = buffered_intervs.begin(); ibuffered_interv != buffered_intervs.end(); ++ibuffered_interv) {
				if (ibuffered_interv->second.bfile().flush())
					TGLError("Failed to write intervals to file %s: %s", ibuffered_interv->second.file_name().c_str(), strerror(errno));
				ibuffered_interv->second.close();
			}

			// read the temporary files one by one into memory, build the quad tree and save it in corresponding track;
			// each buffer file is removed once it is read (on a Hi-C track these are the size of the track)
			progress.init(buffered_intervs.size(), 1);
			for (map<ChromPair, BufferedIntervals2D>::iterator ibuffered_interv = buffered_intervs.begin(); ibuffered_interv != buffered_intervs.end(); ++ibuffered_interv) {
				BufferedIntervals2D &buffered_interv = ibuffered_interv->second;
				int chromid1 = buffered_interv.chromid1();
				int chromid2 = buffered_interv.chromid2();
				string buffer_filename = buffered_interv.file_name();

				snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_2d_filename(iu.get_chromkey(), chromid1, chromid2).c_str());

				// A POINTS source gives a POINTS track
				if (src_track_type == GenomeTrack::POINTS) {
					write_2d_points(buffer_filename, filename, dirname, chromid1, chromid2, iu, agg_cfg);
					progress.report(1);
					continue;
				}

				// Collect the mapped rects, then aggregate overlaps into disjoint
				// rects before inserting (the quadtree forbids overlapping objects).
				vector<RectVal2D> rects;
				buffered_interv.open(buffer_filename.c_str(), "r", chromid1, chromid2);
				while (buffered_interv.read_interval()) {
					const GInterval2D &gi = buffered_interv.last_interval();
					RectVal2D rv;
					rv.x1 = gi.start1(); rv.x2 = gi.end1();
					rv.y1 = gi.start2(); rv.y2 = gi.end2();
					rv.v = buffered_interv.last_val();
					rects.push_back(rv);
				}
				buffered_interv.close();
				unlink(buffer_filename.c_str());

				RectsQuadTree qtree;

				qtree.reset(0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2));
				aggregate_2d_rects(rects, agg_cfg, qtree);

				if (!qtree.empty()) {
					GenomeTrackRectsRects gtrack(iu.get_track_chunk_size(), iu.get_track_num_chunks());

					gtrack.init_write(filename, chromid1, chromid2);
					gtrack.write(qtree);
				}

				progress.report(1);
			}

			progress.report_last();
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
