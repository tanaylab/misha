/*
 * SeqOneHot.h
 *
 * One-hot encoding of genomic windows read straight from the sequence store, for
 * gseq.onehot / gextract.windows (R misha) and gseq_onehot / gextract_windows (pymisha).
 *
 * THIS FILE IS SHARED VERBATIM BETWEEN misha AND pymisha: keep the two copies byte-identical.
 * It uses only GenomeSeqFetch::read_interval, GenomeChromKey and GInterval's 4-argument
 * constructor, which both repositories have, and reports interrupts and errors through the
 * callers (a functor and TGLException), so it builds against either entry-point layer.
 *
 * Encoding: channel 0 = A, 1 = C, 2 = G, 3 = T, case-insensitive. Any other byte (N, IUPAC
 * ambiguity codes, gaps) and every position outside the chromosome sets no channel, i.e. its
 * column is all zero. A reverse-complemented window reads the bases from its end and swaps
 * A <-> T and C <-> G, so it is the encoding of the reverse complement of the forward sequence.
 *
 * Layout: element (window n, channel c, position p) is written to out[n * s.n + c * s.c + p * s.p],
 * so one routine fills a C-order (N, 4, L) numpy array (s = {4L, L, 1}) and a column-major R
 * array of dim c(N, 4, L) (s = {1, N, 4N}) alike. Every element of a window is written, zeros
 * included, so the output need not be initialized.
 */

#ifndef SEQONEHOT_H_
#define SEQONEHOT_H_

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "GInterval.h"
#include "GenomeChromKey.h"
#include "GenomeSeqFetch.h"

namespace seqonehot {

// Channel of a byte: A/a 0, C/c 1, G/g 2, T/t 3, anything else -1.
inline int base2channel(unsigned char c)
{
	switch (c) {
	case 'A': case 'a': return 0;
	case 'C': case 'c': return 1;
	case 'G': case 'g': return 2;
	case 'T': case 't': return 3;
	default: return -1;
	}
}

// channel_table(c)[byte] is 1 if base2channel(byte) == c, else 0: the encoding is a table
// lookup per channel, with no branch on the base.
inline const unsigned char *channel_table(int c)
{
	struct Tables {
		unsigned char t[4][256];
		Tables() {
			for (int b = 0; b < 256; ++b)
				for (int ch = 0; ch < 4; ++ch)
					t[ch][b] = base2channel((unsigned char)b) == ch;
		}
	};
	static const Tables tables;
	return tables.t[c];
}

struct Strides {
	int64_t n;  // between windows
	int64_t c;  // between channels
	int64_t p;  // between positions
};

// Encodes seq[0, len) as window n: writes all four channels of positions [0, len).
template <typename T>
inline void encode(const char *seq, int64_t len, bool revcomp, int64_t n, T *out, const Strides &s)
{
	const unsigned char *useq = (const unsigned char *)seq;

	for (int c = 0; c < 4; ++c) {
		const unsigned char *tab = channel_table(revcomp ? 3 - c : c);
		T *row = out + n * s.n + c * s.c;
		if (revcomp) {
			for (int64_t p = 0; p < len; ++p)
				row[p * s.p] = (T)tab[useq[len - 1 - p]];
		} else {
			for (int64_t p = 0; p < len; ++p)
				row[p * s.p] = (T)tab[useq[p]];
		}
	}
}

// Reads [start, start + width) of chromid into buf[0, width). The part outside [0, chrom size)
// is 'N'. tmp is scratch space for GenomeSeqFetch. The read is always on the forward strand.
inline void read_padded(GenomeSeqFetch &seqfetch, const GenomeChromKey &chromkey, int chromid,
						int64_t start, int64_t width, std::vector<char> &tmp, char *buf)
{
	int64_t chromsize = (int64_t)chromkey.get_chrom_size(chromid);
	int64_t cs = std::max<int64_t>(start, 0);
	int64_t ce = std::min<int64_t>(start + width, chromsize);

	if (cs >= ce) {
		memset(buf, 'N', width);
		return;
	}

	seqfetch.read_interval(GInterval(chromid, cs, ce, 0), chromkey, tmp);

	int64_t left = cs - start;
	int64_t core = std::min<int64_t>((int64_t)tmp.size(), ce - cs);

	if (left > 0)
		memset(buf, 'N', left);
	memcpy(buf + left, tmp.data(), core);
	if (left + core < width)
		memset(buf + left + core, 'N', width - left - core);
}

// Number of windows read and encoded together: their sequences are held in one buffer of
// block * width bytes (at most ~16 MB). Within a block the reads go in genomic order, which
// keeps a per-chromosome sequence file open and the read-ahead window warm however the
// windows are ordered, while the output is written in the order of the windows.
inline int64_t block_size(int64_t width)
{
	return std::max<int64_t>(1, std::min<int64_t>(4096, (int64_t(1) << 24) / std::max<int64_t>(width, 1)));
}

// Encodes num windows of the given width: window i is [starts[i], starts[i] + width) of chromosome
// chromids[i], reverse complemented if revcomp[i] (revcomp may be NULL), written as window
// first_n + i. check_interrupt() is called between blocks; it may throw.
template <typename T, typename Interrupt>
void encode_windows(GenomeSeqFetch &seqfetch, const GenomeChromKey &chromkey,
					const int *chromids, const int64_t *starts, const unsigned char *revcomp,
					int64_t num, int64_t width, int64_t first_n, T *out, const Strides &s,
					Interrupt check_interrupt)
{
	if (num <= 0 || width <= 0)
		return;

	int64_t block = block_size(width);
	std::vector<char> buf(std::min(block, num) * width);
	std::vector<char> tmp;
	std::vector<int64_t> order;

	for (int64_t b0 = 0; b0 < num; b0 += block) {
		int64_t nb = std::min(block, num - b0);

		order.resize(nb);
		for (int64_t k = 0; k < nb; ++k)
			order[k] = k;
		std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
			if (chromids[b0 + a] != chromids[b0 + b])
				return chromids[b0 + a] < chromids[b0 + b];
			return starts[b0 + a] < starts[b0 + b];
		});

		for (int64_t k = 0; k < nb; ++k) {
			int64_t i = b0 + order[k];
			read_padded(seqfetch, chromkey, chromids[i], starts[i], width, tmp, &buf[order[k] * width]);
		}

		if (s.p == 1) {
			// positions contiguous (numpy): one window at a time
			for (int64_t k = 0; k < nb; ++k)
				encode(&buf[k * width], width, revcomp && revcomp[b0 + k], first_n + b0 + k, out, s);
		} else {
			// windows contiguous (R): one position and channel of every window of the block at a
			// time, so that consecutive writes land next to each other
			T *base = out + (first_n + b0) * s.n;
			for (int64_t p = 0; p < width; ++p) {
				for (int c = 0; c < 4; ++c) {
					const unsigned char *fwd = channel_table(c);
					const unsigned char *rev = channel_table(3 - c);
					T *dst = base + c * s.c + p * s.p;
					for (int64_t k = 0; k < nb; ++k) {
						const unsigned char *seq = (const unsigned char *)&buf[k * width];
						bool rc = revcomp && revcomp[b0 + k];
						dst[k * s.n] = (T)(rc ? rev[seq[width - 1 - p]] : fwd[seq[p]]);
					}
				}
			}
		}

		check_interrupt();
	}
}

}  // namespace seqonehot

#endif /* SEQONEHOT_H_ */
