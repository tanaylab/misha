// SeqShuffle.cpp - see SeqShuffle.h. Shared byte for byte between misha and pymisha.

#include <algorithm>

#include "SeqShuffle.h"

namespace seqshuffle {

namespace {

inline uint64_t splitmix64(uint64_t &x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

inline uint64_t rotl(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

// 0..3 for A, C, G, T in either case, -1 for anything else
inline int base_code(char c)
{
	switch (c) {
	case 'A': case 'a': return 0;
	case 'C': case 'c': return 1;
	case 'G': case 'g': return 2;
	case 'T': case 't': return 3;
	default: return -1;
	}
}

const char UPPER[4] = {'A', 'C', 'G', 'T'};
const char LOWER[4] = {'a', 'c', 'g', 't'};

// Fisher-Yates over v[0, n), drawing from the end down
template <typename T>
void shuffle_range(T *v, size_t n, Rng &rng)
{
	for (size_t i = n; i > 1; --i) {
		size_t j = (size_t)rng.below(i);
		std::swap(v[i - 1], v[j]);
	}
}

} // namespace

//------------------------------------------------------------------------------------------

Rng::Rng(int64_t seed, uint64_t stream)
{
	uint64_t x = (uint64_t)seed;
	uint64_t base = splitmix64(x);
	uint64_t y = base ^ (stream * 0xD1B54A32D192ED03ULL);
	for (int i = 0; i < 4; ++i)
		m_s[i] = splitmix64(y);
}

uint64_t Rng::next()
{
	const uint64_t result = rotl(m_s[1] * 5, 7) * 9;
	const uint64_t t = m_s[1] << 17;
	m_s[2] ^= m_s[0];
	m_s[3] ^= m_s[1];
	m_s[1] ^= m_s[2];
	m_s[0] ^= m_s[3];
	m_s[2] ^= t;
	m_s[3] = rotl(m_s[3], 45);
	return result;
}

uint64_t Rng::below(uint64_t n)
{
	// reject the lowest (2^64 mod n) values, so the rest divide evenly by n
	const uint64_t threshold = (0 - n) % n;
	for (;;) {
		uint64_t r = next();
		if (r >= threshold)
			return r % n;
	}
}

//------------------------------------------------------------------------------------------

Shuffler::Shuffler(const char *seq, size_t len, int k) : m_seq(seq), m_len(len), m_k(k)
{
	size_t max_len = 0;
	size_t i = 0;
	while (i < len) {
		if (base_code(seq[i]) < 0) {
			++i;
			continue;
		}
		size_t j = i;
		while (j < len && base_code(seq[j]) >= 0)
			++j;
		const size_t seg_len = j - i;
		// A segment of at most k bases holds at most one k-mer: it is its only
		// arrangement. (For k = 1 that is a single base.)
		if (seg_len > (size_t)k) {
			m_segments.push_back(Segment());
			Segment &seg = m_segments.back();
			seg.start = i;
			seg.codes.resize(seg_len);
			for (size_t p = 0; p < seg_len; ++p)
				seg.codes[p] = (uint8_t)base_code(seq[i + p]);
			seg.num_nodes = 0;
			seg.end_node = 0;
			if (k >= 2)
				build_graph(seg);
			max_len = std::max(max_len, seg_len);
		}
		i = j;
	}

	m_codes.resize(max_len);
	if (k >= 2) {
		m_in_tree.resize(max_len);
		m_tree_edge.resize(max_len);
		m_order.resize(max_len);
		m_next.resize(max_len);
	}
}

void Shuffler::build_graph(Segment &seg)
{
	const size_t L = seg.codes.size();
	const size_t num_edges = L - (size_t)m_k + 1;  // >= 2
	const size_t num_pos = num_edges + 1;          // (k-1)-mers
	const int w = m_k - 1;
	const uint64_t mask = (w >= 32) ? ~0ULL : ((1ULL << (2 * w)) - 1);

	// the (k-1)-mer starting at every position, packed two bits per base
	std::vector<uint64_t> kmer(num_pos);
	uint64_t code = 0;
	for (int p = 0; p < w - 1; ++p)
		code = (code << 2) | seg.codes[p];
	for (size_t p = 0; p < num_pos; ++p) {
		code = ((code << 2) | seg.codes[p + w - 1]) & mask;
		kmer[p] = code;
	}

	// vertices numbered in the order of their (k-1)-mers, which fixes the order
	// Wilson's algorithm visits them in
	std::vector<uint64_t> sorted(kmer);
	std::sort(sorted.begin(), sorted.end());
	sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
	seg.num_nodes = (int32_t)sorted.size();
	seg.node.resize(num_pos);
	for (size_t p = 0; p < num_pos; ++p)
		seg.node[p] = (int32_t)(std::lower_bound(sorted.begin(), sorted.end(), kmer[p]) - sorted.begin());
	seg.end_node = seg.node[num_edges];

	seg.out_start.assign(seg.num_nodes + 1, 0);
	for (size_t e = 0; e < num_edges; ++e)
		++seg.out_start[seg.node[e] + 1];
	for (int32_t v = 0; v < seg.num_nodes; ++v)
		seg.out_start[v + 1] += seg.out_start[v];
	seg.out_edges.resize(num_edges);
	std::vector<int32_t> fill(seg.out_start.begin(), seg.out_start.end() - 1);
	for (size_t e = 0; e < num_edges; ++e)
		seg.out_edges[fill[seg.node[e]]++] = (int32_t)e;
}

void Shuffler::shuffle(Rng &rng, char *out)
{
	std::copy(m_seq, m_seq + m_len, out);
	for (const Segment &seg : m_segments)
		shuffle_segment(seg, rng, out);
}

void Shuffler::shuffle_segment(const Segment &seg, Rng &rng, char *out)
{
	const size_t L = seg.codes.size();
	uint8_t *res = m_codes.data();

	if (m_k == 1) {
		std::copy(seg.codes.begin(), seg.codes.end(), res);
		shuffle_range(res, L, rng);
	} else {
		const int32_t V = seg.num_nodes;
		const int32_t root = seg.end_node;
		const int32_t *node = seg.node.data();
		const int32_t *out_start = seg.out_start.data();
		const int32_t *out_edges = seg.out_edges.data();
		char *in_tree = m_in_tree.data();
		int32_t *tree_edge = m_tree_edge.data();
		int32_t *order = m_order.data();
		int32_t *next = m_next.data();

		// Wilson's algorithm: a uniformly random spanning tree of last exits,
		// directed towards the end vertex. Every vertex but the end one has an
		// out-edge, and every vertex reaches the end one, so the walks stop.
		std::fill(in_tree, in_tree + V, 0);
		in_tree[root] = 1;
		for (int32_t v = 0; v < V; ++v) {
			for (int32_t u = v; !in_tree[u];) {
				const int32_t deg = out_start[u + 1] - out_start[u];
				const int32_t e = out_edges[out_start[u] + (deg > 1 ? (int32_t)rng.below(deg) : 0)];
				tree_edge[u] = e;  // a revisit overwrites it: that erases the loop
				u = node[e + 1];
			}
			for (int32_t u = v; !in_tree[u]; u = node[tree_edge[u] + 1])
				in_tree[u] = 1;
		}

		// the order each vertex is left in: the other out-edges at random, then
		// the tree edge
		for (int32_t v = 0; v < V; ++v) {
			const int32_t b = out_start[v];
			const int32_t deg = out_start[v + 1] - b;
			std::copy(out_edges + b, out_edges + b + deg, order + b);
			int32_t free_edges = deg;
			if (v != root) {
				int32_t *last = std::find(order + b, order + b + deg, tree_edge[v]);
				std::swap(*last, order[b + deg - 1]);
				--free_edges;
			}
			shuffle_range(order + b, (size_t)free_edges, rng);
			next[v] = b;
		}

		// the walk: the first k-1 bases stay, each edge adds its last base
		const size_t w = (size_t)m_k - 1;
		std::copy(seg.codes.begin(), seg.codes.begin() + w, res);
		int32_t cur = node[0];
		for (size_t p = w; p < L; ++p) {
			const int32_t e = order[next[cur]++];
			res[p] = seg.codes[e + w];
			cur = node[e + 1];
		}
	}

	const char *in = m_seq + seg.start;
	char *dst = out + seg.start;
	for (size_t p = 0; p < L; ++p) {
		const bool lower = in[p] >= 'a' && in[p] <= 'z';
		dst[p] = lower ? LOWER[res[p]] : UPPER[res[p]];
	}
}

} // namespace seqshuffle
