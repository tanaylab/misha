#ifndef BIGWIGWRITER_H_
#define BIGWIGWRITER_H_

#include <cstdint>
#include <string>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------
// BigWigWriter writes a bigWig file (https://genome.ucsc.edu/goldenPath/help/bigWig.html; Kent et al.
// 2010, Bioinformatics 26:2204): header, chromosome B+ tree, zlib-compressed data blocks, R-tree
// index, total summary and zoom levels.
//
// The file is written in three steps, each of which opens the file, writes and closes it, so that
// each can be a call of its own from R or Python and nothing stays open between them. What a later
// step needs is carried in "state", bytes that the caller hands from one step to the next unchanged:
//
//   create()  writes everything that precedes the data blocks;
//   append()  writes data blocks, once per batch of records (as many batches as needed);
//   finish()  writes the index, the zoom levels, the total summary and the header.
//
// The chromosomes are given in the order of the file: sorted by name, byte by byte (C locale), as
// the B+ tree that maps names to chromosome ids is keyed. Chromosome id i is the i-th of them, and
// the records must come in the order of the ids, then of the starts, and must not overlap. A
// database whose chromosomes are 1, 2, ..., 10 is therefore written 1, 10, 2, ....
//
// Records hold float32 values; NaN values are not written. A run of at least MIN_FIXED_RUN records
// of equal width at a constant step (the bins of a fixed-bin track or iterator) is written as
// fixedStep sections, 4 bytes per record; other records as bedGraph sections, 12 bytes per record.
// Shorter runs are not worth it: each one would cost a data block of its own and split the
// bedGraph block around it (a dense track with scattered NaN bins came out larger with runs of 16).
//
// The zoom levels are libBigWig's choice of sizes: 16 times the mean record width (at least 10
// bases), then 4 times the previous one, up to the longest chromosome and at most MAX_ZOOM_LEVELS
// of them; a level that would hold no fewer records than the previous one is skipped. A zoom
// record summarizes the bases of the records inside one bin of the level's size, the bins being
// aligned to the chromosome start, so that each level is computed exactly from the one below it.
//
// The layout follows libBigWig's and UCSC's writers; the bigWig writer of figwig (Python, by Jacob
// Schreiber, MIT license) served as a reference for it.
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class BigWigWriter {
public:
	enum Errors { FILE_ERROR, BAD_ARG, BAD_RECORD, BAD_STATE };

	struct Chrom {
		string   name;
		uint64_t size;
	};

	static const unsigned MAX_ZOOM_LEVELS = 10;
	static const unsigned MIN_FIXED_RUN = 128;

	// Creates (or truncates) filename and writes the header placeholders and the chromosome tree.
	// chroms must be sorted by name, byte by byte, with unique non-empty names and sizes from 1 to
	// 2^32 - 1.
	static void create(const string &filename, const vector<Chrom> &chroms, vector<char> &state);

	// Writes records i = 0..n-1: chromosome id chromids[i] (an index into the chroms of create()),
	// [starts[i], ends[i]), values[i] (written as float32; NaN is skipped). Records must be sorted by
	// chromosome id, then start, must not overlap (records of earlier batches included), and must lie
	// within their chromosome. Data blocks are compressed on up to num_threads threads.
	static void append(const string &filename, vector<char> &state, uint64_t n, const int64_t *chromids,
	                   const int64_t *starts, const int64_t *ends, const double *values, unsigned num_threads);

	// Writes the index, the zoom levels, the total summary and the header. Zoom blocks are compressed
	// on up to num_threads threads. Returns the number of records written.
	static uint64_t finish(const string &filename, const vector<char> &state, unsigned num_threads);
};

#endif /* BIGWIGWRITER_H_ */
