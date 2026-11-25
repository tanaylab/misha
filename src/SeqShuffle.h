// SeqShuffle.h - shuffling a DNA sequence while keeping its k-mer counts exactly.
//
// Shared byte for byte between misha (gseq.shuffle) and pymisha (gseq_shuffle), so
// that the two give the same output for the same seed. Nothing here touches the R
// or the Python API: the random numbers come from the self-contained generator
// below, never from R's or numpy's.
//
// What is kept. A/C/G/T, in either case, are shuffled; every other byte (N, IUPAC
// codes, gaps, ...) stays where it is and splits the sequence into segments, each
// shuffled on its own. Bases are compared ignoring case, and every position keeps
// the case it had, so a soft-masked stretch stays masked at the same coordinates.
// Within a segment the shuffle keeps the count of every k-mer and the first k-1
// bases - and with them, since the k-mers form a path through the (k-1)-mers, the
// last k-1 bases too. The k-mers that run into a fixed byte are unchanged, so the
// counts of the whole sequence are kept as well.
//
// The draw is uniform over all the sequences that satisfy those constraints, by
// the Eulerian-path construction of Altschul & Erickson (1985) and Kandel et al.
// (1996), as in Jiang et al.'s uShuffle (2008). The (k-1)-mers are the vertices
// of a multigraph and every k-mer of the segment is an edge from its prefix to
// its suffix, so a sequence with the same k-mer counts and the same first (k-1)-
// mer is an Eulerian path from the first (k-1)-mer to the last. A path is fixed
// by the last edge it leaves each vertex through - these form a spanning tree
// directed towards the end vertex - together with the order of the other edges
// out of each vertex (the BEST theorem). So: draw that tree uniformly, by
// Wilson's loop-erased random walk, put every vertex's other out-edges in a
// uniformly random order before its tree edge, and walk. Every edge-labelled
// path is equally likely, and every sequence has the same number of them (the
// product of the factorials of its k-mer multiplicities), so every sequence is
// too. Keeping each vertex's original last edge instead of drawing the tree - a
// common shortcut - gives a valid sequence but not a uniform one: sequences
// whose last exits differ from the input's are never drawn.

#ifndef SEQSHUFFLE_H_
#define SEQSHUFFLE_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seqshuffle {

// Largest k: a (k-1)-mer is packed into 64 bits, two per base.
const int MAX_K = 32;

// xoshiro256** (Blackman & Vigna), seeded through SplitMix64. Portable integer
// arithmetic only, so a seed gives the same stream on every platform and in
// both packages.
class Rng {
public:
	// One independent stream per (seed, stream) pair: gseq.shuffle uses the index
	// of the input sequence as the stream, so a sequence's shuffles do not depend
	// on the sequences around it.
	Rng(int64_t seed, uint64_t stream);

	uint64_t next();

	// Uniform integer in [0, n), without modulo bias. n must be > 0.
	uint64_t below(uint64_t n);

private:
	uint64_t m_s[4];
};

// The k-mer graph of one sequence, built once and shuffled any number of times.
class Shuffler {
public:
	// seq holds len bytes; 1 <= k <= MAX_K (not checked here).
	Shuffler(const char *seq, size_t len, int k);

	// Writes one shuffle of the sequence (len bytes) to out.
	void shuffle(Rng &rng, char *out);

private:
	struct Segment {
		size_t start;                  // offset in the sequence
		std::vector<uint8_t> codes;    // the bases, 0..3 for A, C, G, T
		// k >= 2 only: the (k-1)-mer graph. Edge e is the k-mer at codes[e]; it
		// leaves vertex node[e] and enters node[e + 1].
		std::vector<int32_t> node;
		std::vector<int32_t> out_start; // CSR of the out-edges of each vertex, by position
		std::vector<int32_t> out_edges;
		int32_t num_nodes;
		int32_t end_node;
	};

	void build_graph(Segment &seg);
	void shuffle_segment(const Segment &seg, Rng &rng, char *out);

	const char *m_seq;
	size_t m_len;
	int m_k;
	std::vector<Segment> m_segments;

	// work space, sized for the largest segment
	std::vector<uint8_t> m_codes;
	std::vector<char> m_in_tree;
	std::vector<int32_t> m_tree_edge;
	std::vector<int32_t> m_order;
	std::vector<int32_t> m_next;
};

} // namespace seqshuffle

#endif // SEQSHUFFLE_H_
