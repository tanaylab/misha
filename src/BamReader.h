#ifndef BAMREADER_H_
#define BAMREADER_H_

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

//----------------------------------------------------------------------------
// BamReader reads the alignment records of a BAM file directly: the BGZF blocks, the header's
// reference names and lengths, then the fixed fields of each record and the reference span of its
// CIGAR. It replaces the `samtools view` pipe misha used to read BAM through.
//
// The BGZF blocks are inflated in batches, on up to nthreads threads (the threads run zlib only),
// and consumed in file order by the calling thread. Finding the block boundaries first and then
// inflating the blocks independently is possible because BGZF blocks are independent deflate
// streams by design (SAM/BAM specification, section 4.1), the same property htslib's threaded
// reading uses.
//
// A partial read cannot pass silently. The file must end with the BGZF end-of-file block, every
// block must inflate (zlib checks it) to exactly the length and CRC32 its trailer declares, the
// blocks must end exactly at the end of the file, and the records must end exactly at the end of
// the inflated stream. A record that htslib (and so samtools) refuses is refused too: a block
// shorter than the fixed fields, field lengths that do not fit in the record, a reference index
// outside the header, or a mapped record whose CIGAR does not consume as many bases as its
// sequence has.
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class BamReader {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, TRUNCATED };

	// The fields of a record that misha reads
	struct Record {
		int32_t  refid;      // index in ref_names(), -1 for none
		int32_t  pos;        // 0-based leftmost aligned base, -1 for none
		uint16_t flag;
		uint8_t  mapq;
		int32_t  next_refid; // index in ref_names() of the mate, -1 for none
		int32_t  next_pos;   // 0-based position of the mate, -1 for none
		int32_t  tlen;
		// Reference length of the CIGAR (its M, D, N, = and X operations): 0 for a record without a
		// CIGAR, -1 for a CIGAR with an operation code that SAM has no letter for (or an absurd length)
		int64_t  ref_span;
	};

	BamReader() {}
	~BamReader();

	BamReader(const BamReader &) = delete;
	BamReader &operator=(const BamReader &) = delete;

	// True if the file starts with the BGZF magic and its inflated stream with "BAM\1"
	static bool is_bam(const char *filename);

	// Opens the file and reads its header. nthreads is the number of threads that inflate blocks
	// (1: none but the calling one).
	void open(const std::string &filename, unsigned nthreads = 1);

	const std::vector<std::string> &ref_names() const { return m_ref_names; }
	const std::vector<int64_t>     &ref_lengths() const { return m_ref_lengths; }

	// Reads the next record. Returns false past the last one.
	bool next(Record &rec);

	// Appends to out the line `samtools view` prints for rec, reduced to the columns misha's SAM
	// parser reads, each with the value that parser takes from samtools' text: FLAG, RNAME, POS,
	// MAPQ, CIGAR, PNEXT and TLEN. QNAME and RNEXT are left empty, SEQ is "*", and the CIGAR is
	// written as "<span>M", which has the same reference span ("*" when that span is not positive).
	void append_sam_line(const Record &rec, std::string &out) const;

	uint64_t num_records() const { return m_num_records; }

private:
	struct Block {
		size_t   offset;     // in m_cbuf
		size_t   size;       // the whole block, header and trailer included
		size_t   out_offset; // in m_data
		uint32_t length;     // inflated length
	};

	std::string                m_filename;
	FILE                      *m_fp{nullptr};
	unsigned                   m_nthreads{1};
	uint64_t                   m_file_size{0};
	uint64_t                   m_cbuf_file_offset{0}; // offset in the file of m_cbuf[0]
	bool                       m_file_done{false};    // every byte of the file was read into m_cbuf
	bool                       m_blocks_done{false};  // every block was inflated
	std::vector<unsigned char> m_cbuf;                // compressed bytes read and not yet inflated
	size_t                     m_cbuf_len{0};
	std::vector<char>          m_data;                // inflated bytes
	size_t                     m_data_pos{0};         // the first byte not yet consumed
	size_t                     m_data_len{0};
	std::vector<Block>         m_blocks;
	std::vector<int>           m_block_status;
	std::vector<std::string>   m_ref_names;
	std::vector<int64_t>       m_ref_lengths;
	uint64_t                   m_num_records{0};

	// Makes n bytes available from m_data + m_data_pos on. Returns false if the stream ends first.
	bool fill(size_t n);
	// Inflates the next batch of blocks onto the end of m_data. Returns false at the end of the file.
	bool inflate_batch();
	void check_eof_block();
	[[noreturn]] void truncated(const char *what) const;
};

#endif /* BAMREADER_H_ */
