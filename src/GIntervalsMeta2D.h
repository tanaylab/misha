#ifndef GINTERVALSMETA2D_H_INCLUDED
#define GINTERVALSMETA2D_H_INCLUDED

#include <algorithm>
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

	// Sparse container for write-side ChromStat accumulation.
	//
	// Phase 7b: replaces the prior dense vector<ChromStat> sized N*N that
	// OOMed on databases with O(1M) contigs (N*N would be 10^12 entries).
	// Only populated chrom-pairs occupy memory. Used by 2D big-set save
	// paths (BinnedTransform / Extract / Quantiles / Partition / Screener /
	// Summary / Iterator) and by gtrack.create_meta for 2D tracks.
	class ChromStats2D {
	public:
		ChromStats2D() : m_max_pairs(0) {}
		void clear() { m_stats.clear(); }

		// Number of populated chrom-pairs (NOT N*N).
		size_t size() const { return m_stats.size(); }
		bool   empty() const { return m_stats.empty(); }

		// Upper bound on populated pairs used to size shared-memory result
		// buffers for child processes. Set by begin_save based on the kid
		// scope's num_chrom_pairs().
		void   set_max_pairs(size_t n) { m_max_pairs = n; }
		size_t max_pairs() const { return m_max_pairs; }

		// Insert / overwrite the stat for (chromid1, chromid2). Empty stats
		// (size == 0) are skipped to match the dense semantics (empty cells
		// were never serialized).
		void set(int chromid1, int chromid2, const ChromStat &stat) {
			if (!stat.size) return;
			m_stats[pack_key(chromid1, chromid2)] = stat;
		}

		// Lookup by chrom-pair. Returns nullptr if pair is not populated.
		const ChromStat *find(int chromid1, int chromid2) const {
			auto it = m_stats.find(pack_key(chromid1, chromid2));
			return it == m_stats.end() ? nullptr : &it->second;
		}

		// Merge populated entries from `other` into this. Later wins on key
		// collision (mirrors the previous dense overwrite semantics).
		void merge(const ChromStats2D &other) {
			for (const auto &kv : other.m_stats)
				m_stats[kv.first] = kv.second;
		}

		// Iteration interface.
		using const_iterator = std::unordered_map<uint64_t, ChromStat>::const_iterator;
		const_iterator begin() const { return m_stats.begin(); }
		const_iterator end()   const { return m_stats.end(); }

		// Decode (chromid1, chromid2) from a packed key returned by begin/end.
		static int key_chrom1(uint64_t k) { return (int)(k >> 32); }
		static int key_chrom2(uint64_t k) { return (int)(k & 0xFFFFFFFFu); }

		// Shared-memory packing helpers. Layout written to shmem:
		//   uint32_t count                    (number of populated entries)
		//   PackedEntry entries[count]        (chrom1, chrom2, ChromStat)
		// Buffer size is bounded by max_pairs_bytes(max_pairs) and the actual
		// kid output must fit within that bound.
		struct PackedEntry {
			uint32_t chromid1;
			uint32_t chromid2;
			ChromStat stat;
		};
		static size_t max_pairs_bytes(size_t max_pairs) {
			return sizeof(uint32_t) + max_pairs * sizeof(PackedEntry);
		}
		void pack(void *&ptr) const;
		void unpack(void *&ptr);

	private:
		std::unordered_map<uint64_t, ChromStat> m_stats;
		size_t m_max_pairs;

		static uint64_t pack_key(int chromid1, int chromid2) {
			return (static_cast<uint64_t>(static_cast<uint32_t>(chromid1)) << 32) |
			       static_cast<uint32_t>(chromid2);
		}
	};

	GIntervalsMeta2D() : GIntervalsFetcher2D(BIGSET2D) {}

	virtual ~GIntervalsMeta2D() {}

	static bool is2d(SEXP meta) { return Rf_length(VECTOR_ELT(meta, 0)) == NUM_STAT_COLS; }

	// Sparse save-meta API (Phase 7b). The dense vector<ChromStat>(N*N) API
	// has been removed because it OOMed at O(1M) contigs.
	static void save_plain_intervals_meta(const char *path, const ChromStats2D &chromstats, const IntervUtils &iu);
	static void save_meta(const char *path, SEXP zeroline, const ChromStats2D &chromstats, const IntervUtils &iu);

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

	// Dense (chromid1*N + chromid2) <-> chrom-pair helpers, retained for
	// backward compatibility (callers may pre-compute a packed index in some
	// internal paths). With both the in-memory representation (Phase 7a) and
	// the write-side API (Phase 7b) sparse, none of the meta/big-set 2D code
	// allocates anything of size N*N anymore.
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
	// Phase 7b: walk only populated chrom-pairs.
	//
	// The previous implementation incremented chromid2 (then chromid1) by 1
	// and was O(N*N) when callers used "advance one, ask size(c1,c2), skip
	// if empty" to iterate the sparse population. On a 1M-contig DB that
	// loop never terminated.
	//
	// We use the sorted populated-pair list maintained by GIntervalsMeta2D
	// (m_pair_keys_sorted) and seek to the first key strictly greater than
	// the caller-supplied (chromid1, chromid2). End-of-iteration is reported
	// by setting chromid1/chromid2 to num_chroms (matches the dense walker's
	// post-end contract).
	const uint64_t num_chroms = m_chromkey->get_num_chroms();

	if (m_pair_keys_sorted.empty()) {
		*chromid1 = (int)num_chroms;
		*chromid2 = (int)num_chroms;
		return false;
	}

	uint64_t cur_key = (static_cast<uint64_t>(static_cast<uint32_t>(*chromid1)) << 32) |
	                   static_cast<uint32_t>(*chromid2);
	auto it = std::upper_bound(m_pair_keys_sorted.begin(), m_pair_keys_sorted.end(), cur_key);
	// Masked copies retain zero-size entries to preserve orig_size offsets;
	// skip them so the iterator only stops on truly populated pairs.
	while (it != m_pair_keys_sorted.end()) {
		auto sit = m_pair_stats.find(*it);
		if (sit != m_pair_stats.end() && sit->second.size)
			break;
		++it;
	}
	if (it == m_pair_keys_sorted.end()) {
		*chromid1 = (int)num_chroms;
		*chromid2 = (int)num_chroms;
		return false;
	}

	*chromid1 = key_chrom1(*it);
	*chromid2 = key_chrom2(*it);
	return (uint64_t)*chromid1 < num_chroms && (uint64_t)*chromid2 < num_chroms;
}

#endif

