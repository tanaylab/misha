#ifndef GINTERVALSMETA2D_H_INCLUDED
#define GINTERVALSMETA2D_H_INCLUDED

#include <cstdint>
#include <unordered_map>
#include "GIntervalsFetcher2D.h"
#include "GIntervalsMeta.h"

class GIntervalsMeta2D : public GIntervalsMeta, public GIntervalsFetcher2D {
public:
	enum StatCols {
		CHROM1_COL, CHROM2_COL, CONTAINS_OVERLAPS_COL, SIZE_COL, SURFACE_COL, NUM_STAT_COLS
	};

	static const char *STAT_COL_NAMES[NUM_STAT_COLS];

	struct ChromStat {
		bool    contains_overlaps;
		uint64_t  size;
		double  surface;

		ChromStat() : contains_overlaps(false), size(0), surface(0.) {}
	};

	GIntervalsMeta2D() : GIntervalsFetcher2D(BIGSET2D) {}

	virtual ~GIntervalsMeta2D() {}

	static bool is2d(SEXP meta) { return Rf_length(VECTOR_ELT(meta, 0)) == NUM_STAT_COLS; }

	static void init_chromstats(vector<ChromStat> &chromstats, const IntervUtils &iu);
	static void save_plain_intervals_meta(const char *path, const vector<ChromStat> &chromstats, const IntervUtils &iu);
	static void save_meta(const char *path, SEXP zeroline, const vector<ChromStat> &chromstats, const IntervUtils &iu);

	//-------------------------------- GIntervalsFetcher2D interface -----------------------------------

	virtual uint64_t size() const { return m_size; }

	virtual uint64_t size(int chromid1, int chromid2) const {
		auto it = m_pair_stats.find(pair_key(chromid1, chromid2));
		return it == m_pair_stats.end() ? 0 : (uint64_t)it->second.size;
	}

	virtual int num_chrom_pairs() const;

	virtual double surface() const { return m_surface; } // complexity: O(1)
	virtual double surface(int chromid1, int chromid2) const {
		auto it = m_pair_stats.find(pair_key(chromid1, chromid2));
		return it == m_pair_stats.end() ? 0. : it->second.surface;
	}

	virtual bool get_next_chroms(int *chromid1, int *chromid2);

protected:
	// Sparse per-pair stats. The dense vector<...>(N*N) representation that
	// was used here previously OOMs on databases with 1M+ contigs (Phase 7a).
	// Only populated chrom-pairs appear in m_pair_stats.
	struct PairStat {
		int64_t size;
		int64_t orig_size;
		double  surface;
		bool    contains_overlaps;

		PairStat() : size(0), orig_size(0), surface(0.), contains_overlaps(false) {}
	};

	// Pack a chrom-pair into a 64-bit key for unordered_map.
	static uint64_t pair_key(int chromid1, int chromid2) {
		return (static_cast<uint64_t>(static_cast<uint32_t>(chromid1)) << 32) |
		        static_cast<uint32_t>(chromid2);
	}
	static int key_chrom1(uint64_t k) { return (int)(k >> 32); }
	static int key_chrom2(uint64_t k) { return (int)(k & 0xFFFFFFFFu); }

	std::unordered_map<uint64_t, PairStat> m_pair_stats;

	// Sorted list of populated pair keys; the order matches a row-major walk
	// of the dense (chromid1, chromid2) index, so iteration order is stable
	// across the dense->sparse switch.
	std::vector<uint64_t> m_pair_keys_sorted;

	// Prefix sum of orig_size across m_pair_keys_sorted (size = N_populated + 1).
	// Used by derived classes to compute global udata offsets without an
	// O(N*N) walk.
	std::vector<int64_t> m_orig_size_prefix;

	uint64_t                       m_size;
	double                       m_surface;
	GenomeChromKey              *m_chromkey;

	void init(const char *name, SEXP meta, const GenomeChromKey &chromkey);
	void init_masked_copy(GIntervalsMeta2D *obj, const set<ChromPair> &chrompairs_mask) const;

	// Rebuild m_pair_keys_sorted (sorted) and m_orig_size_prefix from
	// m_pair_stats. O(P log P) in the number of populated pairs.
	void rebuild_pair_index();

	// True if a populated entry exists for (c1, c2).
	bool has_pair(int chromid1, int chromid2) const {
		return m_pair_stats.find(pair_key(chromid1, chromid2)) != m_pair_stats.end();
	}

	// Lookup size (returns 0 if pair missing).
	int64_t pair_size(int chromid1, int chromid2) const {
		auto it = m_pair_stats.find(pair_key(chromid1, chromid2));
		return it == m_pair_stats.end() ? 0 : it->second.size;
	}

	int64_t pair_orig_size(int chromid1, int chromid2) const {
		auto it = m_pair_stats.find(pair_key(chromid1, chromid2));
		return it == m_pair_stats.end() ? 0 : it->second.orig_size;
	}

	// Legacy dense (chromid1*N + chromid2) <-> chrom-pair helpers, retained
	// because the static write-side API (vector<ChromStat> sized N*N) still
	// uses them. The arithmetic itself is cheap; only the dense allocation
	// was the OOM risk and that has been removed from the in-memory layout.
	int chroms2idx(int chromid1, int chromid2) const { return chromid1 * m_chromkey->get_num_chroms() + chromid2; }
	int idx2chrom1(int idx) const { return idx / m_chromkey->get_num_chroms(); }
	int idx2chrom2(int idx) const { return idx % m_chromkey->get_num_chroms(); }
};

//------------------------------------- IMPLEMENTATION ------------------------------------------

inline int GIntervalsMeta2D::num_chrom_pairs() const
{
	int res = 0;
	for (const auto &kv : m_pair_stats) {
		if (kv.second.size)
			++res;
	}
	return res;
}

inline bool GIntervalsMeta2D::get_next_chroms(int *chromid1, int *chromid2)
{
	if ((uint64_t)*chromid2 < m_chromkey->get_num_chroms() - 1)
		++*chromid2;
	else {
		++*chromid1;
		*chromid2 = 0;
	}
	return (uint64_t)*chromid1 < m_chromkey->get_num_chroms() && (uint64_t)*chromid2 < m_chromkey->get_num_chroms();
}

#endif

