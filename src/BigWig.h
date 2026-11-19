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
// levels) and returns its records per chromosome, the same way Wig does for a text file. It also
// reads a bigBed, which keeps its records in the same container (header, chromosome tree, R-tree
// index, zlib-compressed data blocks).
//
// A partial read cannot pass silently: every data block the header counts must be reachable from
// the index (bigWig), all the records the header counts must be read (bigBed), and every block must
// inflate (zlib checks it) to exactly the items its header declares.
// The header's total summary is not used as a check: UCSC's wigToBigWig writes a validCount larger
// than the bases its fixedStep records cover when span < step.
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class BigWig {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, BAD_COORD, INCOMPLETE };

	// A bigBed record field past chrom, start and end: its name, and whether its autoSql type is a number
	struct Field {
		string name;
		bool   numeric;
	};

	// bigBed records: intervals (database chromid, start, end) and, for each field of get_fields(),
	// its column, in values if the field is numeric, in texts otherwise
	struct BedRecords {
		GIntervals              intervals;
		vector<vector<double> > values;
		vector<vector<string> > texts;
	};

	// true if the file starts with the bigWig magic number or is named *.bw / *.bigWig; for a
	// misnamed or empty file, init() then fails instead of the import yielding an empty track
	static bool is_bigwig(const char *filename);

	// Opens a bigWig, or a bigBed if bigbed is true
	void init(const GenomeChromKey &chromkey, const string &filename, bool bigbed = false);

	// bigWig: fills intervals with the records of data block iblock of the chromosome, (float)udata
	// holding the value. The blocks of a chromosome must be read in increasing order. Returns false
	// past the chromosome's last block. With query (sorted, non-overlapping intervals of the
	// chromosome), only the records that overlap it are kept.
	bool get_data(int chromid, uint64_t iblock, GIntervals &intervals, const GIntervals *query = NULL);

	// bigWig: fills intervals with the records of the chromosome that overlap query (sorted,
	// non-overlapping intervals of the chromosome), each once and whole. Only the data blocks the
	// index places over query are read.
	void get_data(int chromid, const GIntervals &query, GIntervals &intervals);

	// bigBed: the fields past chrom, start and end, named by the file's autoSql, or by the BED
	// standard names if the file has none
	const vector<Field> &get_fields() const { return m_fields; }

	// bigBed: reads all the records, in file order. Records on chromosomes the database lacks are
	// skipped (get_unknown_chroms()).
	void read_bed(BedRecords &records);

	// Chromosome names that carry data in the file but do not exist in the genome database.
	const UnknownChroms &get_unknown_chroms() const { return m_unknown_chroms; }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

	// Records read so far that ran past the end of their chromosome: clipped to it, or dropped
	// when they start at or past it; and the first of them in the file, "(chrom, start, end),
	// chromosome size N" ("" if none)
	uint64_t get_num_clipped() const { return m_num_clipped; }
	uint64_t get_num_dropped() const { return m_num_dropped; }
	const string &get_first_clipped() const { return m_first_clipped; }

private:
	struct Block {
		uint64_t offset;
		uint64_t size;
		int64_t  start;    // of the records in the block, on its first chromosome
		int64_t  end;
	};

	const GenomeChromKey  *m_chromkey;
	BufferedFile           m_bfile;
	bool                   m_bigbed{false};
	uint32_t               m_uncompress_bufsize;
	uint64_t               m_data_offset;
	vector<vector<Block> > m_blocks;          // data blocks of each file chromosome (the first, for a bigBed block), in file order
	vector<string>         m_bw_names;
	vector<int>            m_db2bw;           // database chromid -> file chromosome index, or -1
	vector<int>            m_bw2db;           // file chromosome index -> database chromid, or -1
	vector<Field>          m_fields;
	UnknownChroms          m_unknown_chroms;
	uint64_t               m_tree_budget;     // see init()
	vector<char>           m_cbuf;
	vector<char>           m_ubuf;
	int                    m_last_chromid{-1};
	int64_t                m_last_end{0};     // end of the previous record of the chromosome
	uint64_t               m_num_clipped{0};
	uint64_t               m_num_dropped{0};
	string                 m_first_clipped;
	uint64_t               m_first_clipped_offset{0};  // of its data block
	uint64_t               m_first_clipped_item{0};    // in the block

	const char *type_name() const { return m_bigbed ? "bigBed" : "bigWig"; }
	void read_at(uint64_t offset, void *buf, uint64_t size);
	void read_chrom_tree(uint64_t offset, uint32_t keysize, int depth, vector<pair<uint32_t, string> > &chroms);
	void read_rtree(uint64_t offset, int depth);
	void read_fields(uint64_t autosql_offset, uint16_t field_count, uint16_t defined_field_count);
	const char *read_block(const Block &block, uint64_t &size);

	// Applies the chromosome end to a record: returns false if it starts at or past it, clips its end otherwise
	bool clip(uint32_t bwid, int chromid, const Block &block, uint64_t item, int64_t start, int64_t &end);
};

#endif /* BIGWIG_H_ */
