#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <new>
#include <system_error>
#include <thread>
#include <zlib.h>

#include "BamReader.h"
#include "TGLException.h"

using namespace std;

namespace {

// The BGZF end-of-file block: an empty block that every BGZF writer appends (SAM/BAM specification,
// section 4.1.2)
const unsigned char BGZF_EOF[28] = {
	0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
	0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

const size_t BGZF_MAX_BLOCK = 65536;  // the largest block, compressed or inflated
const size_t BATCH_PER_THREAD = 1 << 20; // compressed bytes read per batch, per inflating thread

inline uint16_t le16(const unsigned char *p) { return (uint16_t)(p[0] | (p[1] << 8)); }
inline uint32_t le32(const unsigned char *p) { return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); }
inline int32_t  le32s(const unsigned char *p) { return (int32_t)le32(p); }

// Size of the BGZF block at p, of which avail bytes are at hand: 0 if more bytes are needed to tell,
// -1 if p is not a BGZF block header (the gzip magic, deflate, FEXTRA, and a BC subfield), or if
// the size it declares cannot hold the header and the trailer
int64_t bgzf_block_size(const unsigned char *p, size_t avail)
{
	if (avail < 12)
		return avail >= 4 && (p[0] != 0x1f || p[1] != 0x8b || p[2] != 8 || !(p[3] & 4)) ? -1 : 0;
	if (p[0] != 0x1f || p[1] != 0x8b || p[2] != 8 || !(p[3] & 4))
		return -1;
	size_t xlen = le16(p + 10);
	if (avail < 12 + xlen)
		return 0;
	for (size_t i = 12; i + 4 <= 12 + xlen; ) {
		size_t slen = le16(p + i + 2);
		if (p[i] == 66 && p[i + 1] == 67 && slen == 2 && i + 6 <= 12 + xlen) {
			size_t size = (size_t)le16(p + i + 4) + 1;
			return size >= 12 + xlen + 8 ? (int64_t)size : -1;
		}
		i += 4 + slen;
	}
	return -1;
}

void append_uint(string &out, uint64_t v)
{
	char buf[24];
	char *p = buf + sizeof(buf);
	do {
		*--p = (char)('0' + v % 10);
		v /= 10;
	} while (v);
	out.append(p, buf + sizeof(buf) - p);
}

void append_int(string &out, int64_t v)
{
	if (v < 0) {
		out.push_back('-');
		append_uint(out, (uint64_t)0 - (uint64_t)v);
	} else
		append_uint(out, (uint64_t)v);
}

}  // namespace

BamReader::~BamReader()
{
	if (m_fp)
		fclose(m_fp);
}

bool BamReader::is_bam(const char *filename)
{
	FILE *fp = fopen(filename, "rb");
	if (!fp)
		return false;
	unsigned char magic[4];
	size_t n = fread(magic, 1, 4, fp);
	fclose(fp);
	if (n != 4 || magic[0] != 0x1f || magic[1] != 0x8b || magic[2] != 0x08 || magic[3] != 0x04)
		return false;
	gzFile gz = gzopen(filename, "rb");
	if (!gz)
		return false;
	char payload[4];
	int len = gzread(gz, payload, 4);
	gzclose(gz);
	return len == 4 && !memcmp(payload, "BAM\1", 4);
}

void BamReader::truncated(const char *what) const
{
	TGLError<BamReader>(TRUNCATED, "BAM file %s is truncated or corrupt: %s", m_filename.c_str(), what);
}

void BamReader::check_eof_block()
{
	unsigned char tail[sizeof(BGZF_EOF)];
	if (m_file_size < sizeof(BGZF_EOF) || fseeko(m_fp, (off_t)(m_file_size - sizeof(BGZF_EOF)), SEEK_SET) ||
	    fread(tail, 1, sizeof(tail), m_fp) != sizeof(tail) || memcmp(tail, BGZF_EOF, sizeof(tail)))
		truncated("it does not end with the BGZF end-of-file block. Rewrite it with `samtools view -b` if it is a whole file");
	if (fseeko(m_fp, 0, SEEK_SET))
		TGLError<BamReader>(FILE_ERROR, "Failed to read BAM file %s: %s", m_filename.c_str(), strerror(errno));
}

void BamReader::open(const string &filename, unsigned nthreads)
{
	m_filename = filename;
	m_nthreads = max(nthreads, 1u);

	m_fp = fopen(filename.c_str(), "rb");
	if (!m_fp)
		TGLError<BamReader>(FILE_ERROR, "Failed to open BAM file %s: %s", filename.c_str(), strerror(errno));

	struct stat st;
	if (fstat(fileno(m_fp), &st))
		TGLError<BamReader>(FILE_ERROR, "Failed to stat BAM file %s: %s", filename.c_str(), strerror(errno));
	if (!S_ISREG(st.st_mode))
		TGLError<BamReader>(FILE_ERROR, "BAM file %s is not a regular file", filename.c_str());
	m_file_size = (uint64_t)st.st_size;

	check_eof_block();

	m_cbuf.resize(max(BATCH_PER_THREAD * m_nthreads, 2 * BGZF_MAX_BLOCK));

	// header: magic, text, references
	if (!fill(12))
		truncated("its header is incomplete");
	const unsigned char *p = (const unsigned char *)&m_data[m_data_pos];
	if (memcmp(p, "BAM\1", 4))
		TGLError<BamReader>(INVALID_FORMAT, "%s is not a BAM file", filename.c_str());
	int32_t l_text = le32s(p + 4);
	if (l_text < 0)
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: invalid header text length %d", filename.c_str(), l_text);
	if (!fill(12 + (size_t)l_text))
		truncated("its header is incomplete");
	p = (const unsigned char *)&m_data[m_data_pos];
	int32_t n_ref = le32s(p + 8 + l_text);
	if (n_ref < 0)
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: invalid number of references %d", filename.c_str(), n_ref);
	m_data_pos += 12 + (size_t)l_text;

	// A reference takes at least 9 bytes (l_name, a name of at least its NUL, l_ref), so the space
	// reserved is bounded by what the bytes inflated so far can hold: a corrupt n_ref (up to 2^31 - 1)
	// would otherwise ask for some 80 GB up front, and fail as "out of memory" rather than as a
	// corrupt file. The loop below reports the header as incomplete once the data runs out.
	size_t refs_hint = min<size_t>((size_t)n_ref, (m_data_len - m_data_pos) / 9);
	m_ref_names.reserve(refs_hint);
	m_ref_lengths.reserve(refs_hint);
	for (int32_t i = 0; i < n_ref; ++i) {
		if (!fill(4))
			truncated("its header is incomplete");
		int32_t l_name = le32s((const unsigned char *)&m_data[m_data_pos]);
		if (l_name < 1)
			TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: invalid name length of reference %d", filename.c_str(), i + 1);
		if (!fill(8 + (size_t)l_name))
			truncated("its header is incomplete");
		const char *name = &m_data[m_data_pos + 4];
		m_ref_names.emplace_back(name, strnlen(name, l_name));
		m_ref_lengths.push_back(le32((const unsigned char *)name + l_name));
		m_data_pos += 8 + (size_t)l_name;
	}
}

bool BamReader::fill(size_t n)
{
	while (m_data_len - m_data_pos < n) {
		if (m_blocks_done || !inflate_batch())
			return false;
	}
	return true;
}

bool BamReader::inflate_batch()
{
	// move the bytes not yet consumed to the front
	if (m_data_pos) {
		memmove(m_data.data(), m_data.data() + m_data_pos, m_data_len - m_data_pos);
		m_data_len -= m_data_pos;
		m_data_pos = 0;
	}

	// read up to a full buffer
	while (!m_file_done && m_cbuf_len < m_cbuf.size()) {
		size_t n = fread(m_cbuf.data() + m_cbuf_len, 1, m_cbuf.size() - m_cbuf_len, m_fp);
		m_cbuf_len += n;
		if (!n) {
			if (ferror(m_fp))
				TGLError<BamReader>(FILE_ERROR, "Failed to read BAM file %s: %s", m_filename.c_str(), strerror(errno));
			m_file_done = true;
		}
	}

	// the blocks that lie whole in the buffer
	m_blocks.clear();
	size_t at = 0, out = m_data_len;
	while (at < m_cbuf_len) {
		int64_t size = bgzf_block_size(&m_cbuf[at], m_cbuf_len - at);
		if (size < 0) {
			char msg[200];
			snprintf(msg, sizeof(msg), "no valid BGZF block at offset %llu", (unsigned long long)(m_cbuf_file_offset + at));
			truncated(msg);
		}
		if (!size || at + (size_t)size > m_cbuf_len)
			break;
		uint32_t length = le32(&m_cbuf[at + size - 4]);
		if (length > BGZF_MAX_BLOCK) {
			char msg[200];
			snprintf(msg, sizeof(msg), "invalid length of the BGZF block at offset %llu", (unsigned long long)(m_cbuf_file_offset + at));
			truncated(msg);
		}
		m_blocks.push_back({ at, (size_t)size, out, length });
		at += size;
		out += length;
	}

	if (m_blocks.empty()) {
		if (m_cbuf_len)
			truncated("its last BGZF block is incomplete"); // the buffer holds more than any whole block
		m_blocks_done = true;
		return false;
	}

	if (m_data.size() < out)
		m_data.resize(max(out, m_data.size() + m_data.size() / 2));

	// inflate the blocks, on up to m_nthreads threads; each takes the next block not taken yet
	m_block_status.assign(m_blocks.size(), 0);
	atomic<size_t> next_block(0);
	auto work = [&]() {
		z_stream zs;
		memset(&zs, 0, sizeof(zs));
		if (inflateInit2(&zs, -15) != Z_OK) {
			for (size_t i; (i = next_block++) < m_blocks.size(); )
				m_block_status[i] = 3;
			return;
		}
		unsigned char dummy;
		for (size_t i; (i = next_block++) < m_blocks.size(); ) {
			const Block &b = m_blocks[i];
			const unsigned char *block = &m_cbuf[b.offset];
			size_t xlen = le16(block + 10);
			unsigned char *dst = b.length ? (unsigned char *)&m_data[b.out_offset] : &dummy;

			inflateReset(&zs);
			zs.next_in = (Bytef *)(block + 12 + xlen);
			zs.avail_in = (uInt)(b.size - 12 - xlen - 8);
			zs.next_out = dst;
			zs.avail_out = b.length;
			int ret = inflate(&zs, Z_FINISH);
			if (ret != Z_STREAM_END || zs.avail_out || zs.avail_in)
				m_block_status[i] = 1;
			else if (crc32(crc32(0L, Z_NULL, 0), dst, b.length) != le32(block + b.size - 8))
				m_block_status[i] = 2;
		}
		inflateEnd(&zs);
	};

	unsigned nthreads = (unsigned)min<size_t>(m_nthreads, m_blocks.size());
	vector<thread> threads;
	try {
		threads.reserve(nthreads ? nthreads - 1 : 0);
		for (unsigned i = 1; i < nthreads; ++i)
			threads.emplace_back(work);
	} catch (const system_error &) {
		// a thread that cannot be started leaves its share to the others
	} catch (const bad_alloc &) {
	}
	work();
	for (auto &t : threads)
		t.join();

	for (size_t i = 0; i < m_blocks.size(); ++i) {
		if (m_block_status[i]) {
			char msg[200];
			snprintf(msg, sizeof(msg), "%s of the BGZF block at offset %llu",
			         m_block_status[i] == 1 ? "failed to inflate the data" :
			         m_block_status[i] == 2 ? "CRC32 mismatch in the data" : "out of memory to inflate the data",
			         (unsigned long long)(m_cbuf_file_offset + m_blocks[i].offset));
			if (m_block_status[i] == 3)
				throw bad_alloc();
			truncated(msg);
		}
	}

	// keep the bytes of an incomplete block for the next batch
	memmove(m_cbuf.data(), m_cbuf.data() + at, m_cbuf_len - at);
	m_cbuf_len -= at;
	m_cbuf_file_offset += at;
	m_data_len = out;
	return true;
}

bool BamReader::next(Record &rec)
{
	if (!fill(4)) {
		if (m_data_len > m_data_pos)
			truncated("its last record is incomplete");
		return false;
	}

	const unsigned char *p = (const unsigned char *)&m_data[m_data_pos];
	int32_t block_size = le32s(p);
	uint64_t irec = m_num_records + 1;
	if (block_size < 32)
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: record %llu is invalid (length %d)", m_filename.c_str(), (unsigned long long)irec, block_size);
	if (!fill(4 + (size_t)block_size))
		truncated("its last record is incomplete");
	p = (const unsigned char *)&m_data[m_data_pos + 4];

	rec.refid = le32s(p);
	rec.pos = le32s(p + 4);
	uint32_t l_read_name = p[8];
	rec.mapq = p[9];
	uint32_t n_cigar = le16(p + 12);
	rec.flag = le16(p + 14);
	int32_t l_seq = le32s(p + 16);
	rec.next_refid = le32s(p + 20);
	rec.next_pos = le32s(p + 24);
	rec.tlen = le32s(p + 28);

	// the checks of htslib's bam_read1 and sam_read1
	if (l_seq < 0 || l_read_name < 1 ||
	    32 + (uint64_t)l_read_name + 4 * (uint64_t)n_cigar + ((uint64_t)l_seq + 1) / 2 + (uint64_t)l_seq > (uint64_t)block_size)
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: record %llu is invalid (its field lengths exceed its length)",
		                    m_filename.c_str(), (unsigned long long)irec);
	if (rec.refid < -1 || rec.refid >= (int32_t)m_ref_names.size() || rec.next_refid < -1 || rec.next_refid >= (int32_t)m_ref_names.size())
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: record %llu refers to a reference that its header does not list",
		                    m_filename.c_str(), (unsigned long long)irec);

	// Reference and query lengths of the CIGAR. A long CIGAR kept in a CG tag (a placeholder
	// "<l_seq>S<span>N" in the record) needs no reading: its span is the placeholder's.
	const int64_t MAX_SPAN = (int64_t)1 << 40; // as the SAM parser's cap
	const unsigned char *cigar = p + 32 + l_read_name;
	int64_t span = 0, qlen = 0;
	bool valid = true;
	for (uint32_t i = 0; i < n_cigar; ++i) {
		uint32_t v = le32(cigar + 4 * i);
		uint32_t op = v & 0xf, len = v >> 4;
		switch (op) {
		case 0: case 7: case 8: // M = X
			span += len;
			qlen += len;
			break;
		case 2: case 3:         // D N
			span += len;
			break;
		case 1: case 4:         // I S
			qlen += len;
			break;
		case 5: case 6:         // H P
			break;
		default:                // B and the codes SAM has no letter for
			valid = false;
		}
	}
	rec.ref_span = valid && span < MAX_SPAN ? span : -1;

	if (n_cigar && l_seq > 0 && !(rec.flag & 0x4) && qlen != l_seq)
		TGLError<BamReader>(INVALID_FORMAT, "BAM file %s: CIGAR and query sequence lengths differ for %.*s (record %llu)",
		                    m_filename.c_str(), (int)strnlen((const char *)p + 32, l_read_name), (const char *)p + 32, (unsigned long long)irec);

	m_data_pos += 4 + (size_t)block_size;
	++m_num_records;
	return true;
}

void BamReader::append_sam_line(const Record &rec, string &out) const
{
	out.push_back('\t');                         // QNAME
	append_uint(out, rec.flag);                  // FLAG
	out.push_back('\t');
	if (rec.refid >= 0)                          // RNAME
		out.append(m_ref_names[rec.refid]);
	else
		out.push_back('*');
	out.push_back('\t');
	append_int(out, (int64_t)rec.pos + 1);       // POS
	out.push_back('\t');
	append_uint(out, rec.mapq);                  // MAPQ
	out.push_back('\t');
	if (rec.ref_span > 0) {                      // CIGAR
		append_int(out, rec.ref_span);
		out.push_back('M');
	} else
		out.push_back('*');
	out.append("\t\t");                          // RNEXT
	append_int(out, (int64_t)rec.next_pos + 1);  // PNEXT
	out.push_back('\t');
	append_int(out, rec.tlen);                   // TLEN
	out.append("\t*\n");                         // SEQ
}
