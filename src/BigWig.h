#ifndef BIGWIG_H_
#define BIGWIG_H_

#include <cstdint>
#include <string>
#include <vector>
#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GIntervals.h"
#include "UnknownChroms.h"

using namespace std;

//----------------------------------------------------------------------------
// BigWig class reads a bigWig file directly (the full-resolution data blocks, never the zoom
// levels) and returns its records per chromosome, the same way Wig does for a text file.
//
// A partial read cannot pass silently: every data block the header counts must be reachable from
// the index, and every block must inflate (zlib checks it) to exactly the items its header declares.
// The header's total summary is not used as a check: UCSC's wigToBigWig writes a validCount larger
// than the bases its fixedStep records cover when span < step.
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class BigWig {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, BAD_COORD, INCOMPLETE };

	// true if the file starts with the bigWig magic number or is named *.bw / *.bigWig; for a
	// misnamed or empty file, init() then fails instead of the import yielding an empty track
	static bool is_bigwig(const char *filename);

	void init(const GenomeChromKey &chromkey, const string &filename);

	// Fills intervals with the records of data block iblock of the chromosome, (float)udata holding
	// the value. Blocks must be read in order from 0. Returns false past the chromosome's last block.
	bool get_data(int chromid, uint64_t iblock, GIntervals &intervals);

	// Chromosome names that carry data in the file but do not exist in the genome database.
	const UnknownChroms &get_unknown_chroms() const { return m_unknown_chroms; }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

private:
	struct Block {
		uint64_t offset;
		uint64_t size;
	};

	const GenomeChromKey  *m_chromkey;
	BufferedFile           m_bfile;
	uint32_t               m_uncompress_bufsize;
	vector<vector<Block> > m_blocks;          // data blocks of each bigWig chromosome, in file order
	vector<string>         m_bw_names;
	vector<int>            m_db2bw;           // database chromid -> bigWig chromosome index, or -1
	UnknownChroms          m_unknown_chroms;
	uint64_t               m_tree_budget;     // see init()
	vector<char>           m_cbuf;
	vector<char>           m_ubuf;
	int64_t                m_last_end{0};     // end of the previous record of the chromosome

	void read_at(uint64_t offset, void *buf, uint64_t size);
	void read_chrom_tree(uint64_t offset, uint32_t keysize, int depth, vector<pair<uint32_t, string> > &chroms);
	void read_rtree(uint64_t offset, int depth);
};

#endif /* BIGWIG_H_ */
