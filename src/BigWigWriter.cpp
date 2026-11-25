#include <errno.h>
#include <sys/types.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <new>
#include <thread>
#include <zlib.h>

#include "BigWigWriter.h"
#include "TGLException.h"
#include "rdbutils.h"

// See BigWigWriter.h. All integers are little-endian, as in BigWig.cpp, the reader.
//
// Ctrl-C is checked (rdb::check_interrupt) in the loops over records and data blocks, on the
// calling thread only: the threads of parallel_for compress blocks and never check it.

namespace {

const uint32_t BIGWIG_MAGIC = 0x888FFC26;
const uint32_t CHROM_TREE_MAGIC = 0x78CA8C91;
const uint32_t RTREE_MAGIC = 0x2468ACE0;
const uint64_t STATE_MAGIC = 0x3157575342574d50ULL;   // the bytes "PMWBSWW1"

// header (64 bytes), MAX_ZOOM_LEVELS zoom headers (24 bytes each), total summary (40 bytes), then
// the chromosome tree
const uint64_t SUMMARY_OFFSET = 64 + 24 * BigWigWriter::MAX_ZOOM_LEVELS;
const uint64_t CHROM_TREE_OFFSET = SUMMARY_OFFSET + 40;

// A data or zoom block is at most BUF_SIZE bytes before compression: a 24-byte section header and
// as many items as fit (libBigWig's layout), or 1023 zoom records of 32 bytes
const uint32_t BUF_SIZE = 32768;
const uint32_t BEDGRAPH_ITEMS = (BUF_SIZE - 24) / 12;   // 2728
const uint32_t FIXED_ITEMS = (BUF_SIZE - 24) / 4;       // 8185
const uint32_t ZOOM_ITEMS = BUF_SIZE / 32;              // 1024
const uint32_t ZOOM_BLOCK_ITEMS = ZOOM_ITEMS - 1;

const uint32_t BPT_BLOCK_SIZE = 256;     // keys per chromosome tree node, as in UCSC's writers
const uint32_t RTREE_BLOCK_SIZE = 256;   // children per R-tree node, as in UCSC's writers
const uint32_t ITEMS_PER_SLOT = 1024;    // an R-tree header field that readers do not use; UCSC's value
const int ZLIB_LEVEL = 6;                // zlib's default, as UCSC's and libBigWig's writers

// Blocks compressed together, on several threads, before they are written
const size_t BLOCKS_PER_FLUSH = 256;

// Records of one chromosome gathered before they are laid out as blocks
const size_t ITEMS_PER_FLUSH = 1 << 20;

// Records appended between two checks for Ctrl-C
const uint64_t ITEMS_PER_INTERRUPT_CHECK = 1 << 16;

enum { BEDGRAPH = 1, FIXEDSTEP = 3 };

template <class T> void put(vector<char> &buf, T v)
{
	char b[sizeof(T)];
	memcpy(b, &v, sizeof(T));
	buf.insert(buf.end(), b, b + sizeof(T));
}

template <class T> T get(const char *p)
{
	T v;
	memcpy(&v, p, sizeof(T));
	return v;
}

// An R-tree leaf item: the block's bounds (startChromIx, startBase, endChromIx, endBase), offset and
// size, as written to the file
struct Leaf {
	uint32_t start_chrom;
	uint32_t start;
	uint32_t end_chrom;
	uint32_t end;
	uint64_t offset;
	uint64_t size;
};
static_assert(sizeof(Leaf) == 32, "an R-tree leaf item is 32 bytes");

struct StateHeader {
	uint64_t magic;
	uint64_t num_chroms;
	uint64_t data_offset;    // of the data block count; the blocks follow it
	uint64_t file_end;       // the size of the file once the step that wrote it ended
	uint64_t num_items;      // records written
	uint64_t valid_bases;    // bases they cover
	double   min, max, sum, sumsq;   // of their values, each base counted once
	int64_t  last_chromid;   // of the last record written, -1 if none
	int64_t  last_end;
	uint64_t num_leaves;     // data blocks
};

// What create() and append() hand to the next step: StateHeader, then the chromosome sizes (uint32),
// their names (uint32 length, then the bytes), and the leaf items of the data blocks
struct State {
	StateHeader      h;
	vector<uint32_t> sizes;
	vector<string>   names;
	vector<Leaf>     leaves;

	void parse(const vector<char> &buf)
	{
		const char *p = buf.data();
		const char *end = p + buf.size();

		auto need = [&](uint64_t n) {
			if ((uint64_t)(end - p) < n)
				TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "Invalid bigWig writer state");
		};

		need(sizeof(h));
		memcpy(&h, p, sizeof(h));
		p += sizeof(h);
		if (h.magic != STATE_MAGIC || h.num_chroms > (uint64_t)(end - p) / 8)
			TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "Invalid bigWig writer state");

		sizes.resize(h.num_chroms);
		memcpy(sizes.data(), p, h.num_chroms * 4);
		p += h.num_chroms * 4;

		names.resize(h.num_chroms);
		for (uint64_t i = 0; i < h.num_chroms; ++i) {
			need(4);
			uint32_t len = get<uint32_t>(p);
			p += 4;
			need(len);
			names[i].assign(p, len);
			p += len;
		}

		if ((uint64_t)(end - p) != h.num_leaves * sizeof(Leaf))
			TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "Invalid bigWig writer state");
		leaves.resize(h.num_leaves);
		memcpy(leaves.data(), p, h.num_leaves * sizeof(Leaf));
	}

	void serialize(vector<char> &buf)
	{
		h.num_chroms = sizes.size();
		h.num_leaves = leaves.size();
		buf.clear();
		buf.reserve(sizeof(h) + sizes.size() * 4 + leaves.size() * sizeof(Leaf));
		buf.insert(buf.end(), (const char *)&h, (const char *)&h + sizeof(h));
		buf.insert(buf.end(), (const char *)sizes.data(), (const char *)(sizes.data() + sizes.size()));
		for (const string &name : names) {
			put<uint32_t>(buf, (uint32_t)name.size());
			buf.insert(buf.end(), name.begin(), name.end());
		}
		buf.insert(buf.end(), (const char *)leaves.data(), (const char *)(leaves.data() + leaves.size()));
	}
};

// A file written (and read) at known offsets; closed by the destructor if an error leaves it open
class File {
public:
	~File() { if (m_fp) fclose(m_fp); }

	void open(const string &name, const char *mode)
	{
		m_name = name;
		m_fp = fopen(name.c_str(), mode);
		if (!m_fp)
			TGLError<BigWigWriter>(BigWigWriter::FILE_ERROR, "Failed to open bigWig file %s: %s", name.c_str(), strerror(errno));
		m_pos = 0;
	}

	uint64_t size()
	{
		if (fseeko(m_fp, 0, SEEK_END))
			error("seek in");
		off_t size = ftello(m_fp);
		if (size < 0)
			error("seek in");
		m_pos = size;
		return size;
	}

	void seek(uint64_t offset)
	{
		if (fseeko(m_fp, (off_t)offset, SEEK_SET))
			error("seek in");
		m_pos = offset;
	}

	uint64_t tell() const { return m_pos; }

	void write(const void *buf, size_t size)
	{
		if (size && fwrite(buf, 1, size, m_fp) != size)
			error("write to");
		m_pos += size;
	}

	void write(const vector<char> &buf) { write(buf.data(), buf.size()); }

	void read(void *buf, size_t size)
	{
		if (size && fread(buf, 1, size, m_fp) != size)
			error("read");
		m_pos += size;
	}

	void close()
	{
		int rc = fclose(m_fp);
		m_fp = NULL;
		if (rc)
			error("write to");
	}

private:
	FILE    *m_fp{NULL};
	string   m_name;
	uint64_t m_pos{0};

	[[noreturn]] void error(const char *what)
	{
		TGLError<BigWigWriter>(BigWigWriter::FILE_ERROR, "Failed to %s bigWig file %s: %s", what, m_name.c_str(),
		                       errno ? strerror(errno) : "unknown error");
	}
};

// Opens a file that create() made and the previous step wrote up to state.h.file_end, and puts the
// position there
void reopen(File &file, const string &filename, const State &state)
{
	file.open(filename, "r+b");
	if (file.size() != state.h.file_end)
		TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "bigWig file %s was changed while it was being written", filename.c_str());
}

// Runs work(i) for i = 0..n-1 on up to num_threads threads (the calling thread included). An
// exception thrown by work is raised here once every thread is joined. Threads that cannot be
// started leave their share to the calling thread.
template <class F> void parallel_for(size_t n, unsigned num_threads, F work)
{
	num_threads = (unsigned)max<size_t>(1, min<size_t>(num_threads, n));

	vector<exception_ptr> errors(num_threads);
	auto run = [&](unsigned t) {
		try {
			for (size_t i = t; i < n; i += num_threads)
				work(i);
		} catch (...) {
			errors[t] = current_exception();
		}
	};

	vector<thread> threads;
	try {
		threads.reserve(num_threads - 1);
		for (unsigned t = 1; t < num_threads; ++t)
			threads.emplace_back(run, t);
	} catch (...) {
	}
	run(0);
	for (unsigned t = threads.size() + 1; t < num_threads; ++t)
		run(t);
	for (thread &th : threads)
		th.join();
	for (exception_ptr &e : errors) {
		if (e)
			rethrow_exception(e);
	}
}

// Blocks, uncompressed, with their leaf items (without offset and size); compressed on several
// threads and written to a file at its position, or kept in memory (blob) with offsets in it
class BlockSink {
public:
	vector<Leaf> leaves;
	vector<char> blob;

	BlockSink(File *file, unsigned num_threads) : m_file(file), m_num_threads(num_threads) {}

	void add(vector<char> &raw, const Leaf &leaf)
	{
		m_raw.push_back(vector<char>());
		m_raw.back().swap(raw);
		m_leaves.push_back(leaf);
		if (m_raw.size() >= BLOCKS_PER_FLUSH)
			flush();
	}

	void flush()
	{
		vector<vector<char> > compressed(m_raw.size());

		parallel_for(m_raw.size(), m_num_threads, [&](size_t i) {
			uLongf size = compressBound(m_raw[i].size());
			compressed[i].resize(size);
			int rc = compress2((Bytef *)compressed[i].data(), &size, (const Bytef *)m_raw[i].data(), m_raw[i].size(), ZLIB_LEVEL);
			if (rc != Z_OK)
				throw bad_alloc();   // compress2 fails only for lack of memory with a large enough buffer
			compressed[i].resize(size);
		});

		for (size_t i = 0; i < m_raw.size(); ++i) {
			Leaf leaf = m_leaves[i];
			leaf.size = compressed[i].size();
			if (m_file) {
				leaf.offset = m_file->tell();
				m_file->write(compressed[i]);
			} else {
				leaf.offset = blob.size();
				blob.insert(blob.end(), compressed[i].begin(), compressed[i].end());
			}
			leaves.push_back(leaf);
		}
		m_raw.clear();
		m_leaves.clear();
	}

private:
	File                   *m_file;
	unsigned                m_num_threads;
	vector<vector<char> >   m_raw;
	vector<Leaf>            m_leaves;
};

// The records of one chromosome, laid out as data blocks: runs of at least MIN_FIXED_RUN records of
// equal width at a constant step as fixedStep sections, the others as bedGraph sections
class DataLayout {
public:
	DataLayout(BlockSink &sink) : m_sink(sink) {}

	void add(uint32_t chromid, uint32_t start, uint32_t end, float value)
	{
		if (chromid != m_chromid)
			flush(true);
		m_chromid = chromid;
		m_starts.push_back(start);
		m_ends.push_back(end);
		m_values.push_back(value);
		if (m_starts.size() >= ITEMS_PER_FLUSH)
			flush(false);
	}

	// Lays out the records gathered. Unless final, the records that the next ones may continue (a run,
	// and the bedGraph records before it that do not fill a block) are kept for the next call.
	void flush(bool final)
	{
		size_t n = m_starts.size();
		size_t i = 0;
		size_t bg = 0;   // the first record not laid out, before i: bedGraph records not filling a block
		size_t keep = n;

		while (i < n) {
			uint32_t span = m_ends[i] - m_starts[i];
			size_t j = i + 1;

			if (j < n && m_ends[j] - m_starts[j] == span) {
				uint32_t step = m_starts[j] - m_starts[i];

				while (j < n && m_ends[j] - m_starts[j] == span && m_starts[j] - m_starts[j - 1] == step)
					++j;
			}

			size_t len = j - i;

			if (j == n && !final) {
				if (len >= BigWigWriter::MIN_FIXED_RUN && len >= FIXED_ITEMS) {
					size_t whole = len - len % FIXED_ITEMS;

					bedgraph(bg, i);
					fixedstep(i, i + whole);
					keep = i + whole;
				} else
					keep = bg;
				break;
			}

			if (len >= BigWigWriter::MIN_FIXED_RUN) {
				bedgraph(bg, i);
				fixedstep(i, j);
				bg = j;
			} else {
				while (j - bg >= BEDGRAPH_ITEMS) {
					bedgraph(bg, bg + BEDGRAPH_ITEMS);
					bg += BEDGRAPH_ITEMS;
				}
			}
			i = j;
		}

		if (final || i >= n) {
			bedgraph(bg, n);
			keep = n;
		}

		m_starts.erase(m_starts.begin(), m_starts.begin() + keep);
		m_ends.erase(m_ends.begin(), m_ends.begin() + keep);
		m_values.erase(m_values.begin(), m_values.begin() + keep);
	}

private:
	BlockSink        &m_sink;
	uint32_t          m_chromid{0};
	vector<uint32_t>  m_starts;
	vector<uint32_t>  m_ends;
	vector<float>     m_values;
	vector<char>      m_raw;

	void section_header(uint32_t start, uint32_t end, uint32_t step, uint32_t span, uint8_t type, uint16_t count)
	{
		m_raw.clear();
		put<uint32_t>(m_raw, m_chromid);
		put<uint32_t>(m_raw, start);
		put<uint32_t>(m_raw, end);
		put<uint32_t>(m_raw, step);
		put<uint32_t>(m_raw, span);
		put<uint8_t>(m_raw, type);
		put<uint8_t>(m_raw, 0);
		put<uint16_t>(m_raw, count);
	}

	void emit(uint32_t start, uint32_t end)
	{
		Leaf leaf{m_chromid, start, m_chromid, end, 0, 0};
		m_sink.add(m_raw, leaf);
	}

	// records [from, to) as bedGraph sections
	void bedgraph(size_t from, size_t to)
	{
		for (size_t b = from; b < to; b += BEDGRAPH_ITEMS) {
			size_t e = min<size_t>(to, b + BEDGRAPH_ITEMS);

			section_header(m_starts[b], m_ends[e - 1], 0, 0, BEDGRAPH, (uint16_t)(e - b));
			m_raw.reserve(24 + 12 * (e - b));
			for (size_t k = b; k < e; ++k) {
				put<uint32_t>(m_raw, m_starts[k]);
				put<uint32_t>(m_raw, m_ends[k]);
				put<float>(m_raw, m_values[k]);
			}
			emit(m_starts[b], m_ends[e - 1]);
		}
	}

	// records [from, to), of equal width at a constant step, as fixedStep sections
	void fixedstep(size_t from, size_t to)
	{
		uint32_t span = m_ends[from] - m_starts[from];
		uint32_t step = to - from > 1 ? m_starts[from + 1] - m_starts[from] : span;

		for (size_t b = from; b < to; b += FIXED_ITEMS) {
			size_t e = min<size_t>(to, b + FIXED_ITEMS);

			section_header(m_starts[b], m_ends[e - 1], step, span, FIXEDSTEP, (uint16_t)(e - b));
			m_raw.reserve(24 + 4 * (e - b));
			for (size_t k = b; k < e; ++k)
				put<float>(m_raw, m_values[k]);
			emit(m_starts[b], m_ends[e - 1]);
		}
	}
};

// Writes an R-tree index over leaves at the file's position (index_offset), as UCSC's cirTree: a
// 48-byte header, then the nodes level by level from the root down, every leaf at the same depth.
// endFileOffset is index_offset, the end of the blocks it indexes.
void write_rtree(File &file, const vector<Leaf> &leaves, uint32_t items_per_slot)
{
	uint64_t index_offset = file.tell();
	uint64_t n = leaves.size();
	vector<char> buf;

	put<uint32_t>(buf, RTREE_MAGIC);
	put<uint32_t>(buf, RTREE_BLOCK_SIZE);
	put<uint64_t>(buf, n);
	put<uint32_t>(buf, n ? leaves.front().start_chrom : 0);
	put<uint32_t>(buf, n ? leaves.front().start : 0);
	put<uint32_t>(buf, n ? leaves.back().end_chrom : 0);
	put<uint32_t>(buf, n ? leaves.back().end : 0);
	put<uint64_t>(buf, index_offset);
	put<uint32_t>(buf, items_per_slot);
	put<uint32_t>(buf, 0);

	// nodes[h]: the number of nodes at level h (0: the leaf nodes); span[h]: the leaf items under
	// each node of level h
	vector<uint64_t> nodes(1, max<uint64_t>(1, (n + RTREE_BLOCK_SIZE - 1) / RTREE_BLOCK_SIZE));
	vector<uint64_t> span(1, RTREE_BLOCK_SIZE);

	while (nodes.back() > 1) {
		nodes.push_back((nodes.back() + RTREE_BLOCK_SIZE - 1) / RTREE_BLOCK_SIZE);
		span.push_back(span.back() * RTREE_BLOCK_SIZE);
	}

	int top = (int)nodes.size() - 1;
	vector<vector<uint64_t> > offsets(nodes.size());
	uint64_t pos = index_offset + buf.size();

	for (int h = top; h >= 0; --h) {
		offsets[h].resize(nodes[h]);
		for (uint64_t j = 0; j < nodes[h]; ++j) {
			uint64_t first = j * span[h];
			uint64_t count = h ? (min(n, first + span[h]) - first + span[h - 1] - 1) / span[h - 1] : min(n, first + span[h]) - first;

			offsets[h][j] = pos;
			pos += 4 + (h ? 24 : 32) * count;
		}
	}

	for (int h = top; h >= 0; --h) {
		for (uint64_t j = 0; j < nodes[h]; ++j) {
			uint64_t first = j * span[h];
			uint64_t last = min(n, first + span[h]);

			if (!h) {
				put<uint8_t>(buf, 1);
				put<uint8_t>(buf, 0);
				put<uint16_t>(buf, (uint16_t)(last - first));
				buf.insert(buf.end(), (const char *)(leaves.data() + first), (const char *)(leaves.data() + last));
			} else {
				uint64_t count = (last - first + span[h - 1] - 1) / span[h - 1];

				put<uint8_t>(buf, 0);
				put<uint8_t>(buf, 0);
				put<uint16_t>(buf, (uint16_t)count);
				for (uint64_t c = 0; c < count; ++c) {
					uint64_t cfirst = first + c * span[h - 1];
					uint64_t clast = min(n, cfirst + span[h - 1]) - 1;

					put<uint32_t>(buf, leaves[cfirst].start_chrom);
					put<uint32_t>(buf, leaves[cfirst].start);
					put<uint32_t>(buf, leaves[clast].end_chrom);
					put<uint32_t>(buf, leaves[clast].end);
					put<uint64_t>(buf, offsets[h - 1][j * RTREE_BLOCK_SIZE + c]);
				}
			}
		}
		file.write(buf);
		buf.clear();
	}
}

// Writes the chromosome B+ tree at the file's position, as UCSC's bptFileBulkIndexToOpenFile: a
// 32-byte header, then the nodes level by level from the root down, each padded to the block size.
void write_chrom_tree(File &file, const vector<BigWigWriter::Chrom> &chroms)
{
	uint64_t n = chroms.size();
	uint32_t block_size = (uint32_t)min<uint64_t>(n, BPT_BLOCK_SIZE);
	uint32_t key_size = 0;
	vector<char> buf;

	for (const BigWigWriter::Chrom &chrom : chroms)
		key_size = max(key_size, (uint32_t)chrom.name.size());

	put<uint32_t>(buf, CHROM_TREE_MAGIC);
	put<uint32_t>(buf, block_size);
	put<uint32_t>(buf, key_size);
	put<uint32_t>(buf, 8);   // valSize: chromId, chromSize
	put<uint64_t>(buf, n);
	put<uint64_t>(buf, 0);

	// nodes[h]: the number of nodes at level h (0: the leaves); span[h]: the keys under each node
	vector<uint64_t> nodes(1, (n + block_size - 1) / block_size);
	vector<uint64_t> span(1, block_size);

	while (nodes.back() > 1) {
		nodes.push_back((nodes.back() + block_size - 1) / block_size);
		span.push_back(span.back() * block_size);
	}

	int top = (int)nodes.size() - 1;
	uint64_t node_size = 4 + (uint64_t)block_size * (key_size + 8);
	vector<uint64_t> level_offset(nodes.size());
	uint64_t pos = file.tell() + buf.size();

	for (int h = top; h >= 0; --h) {
		level_offset[h] = pos;
		pos += nodes[h] * node_size;
	}

	auto put_key = [&](const string &name) {
		buf.insert(buf.end(), name.begin(), name.end());
		buf.insert(buf.end(), key_size - name.size(), '\0');
	};

	for (int h = top; h >= 0; --h) {
		for (uint64_t j = 0; j < nodes[h]; ++j) {
			uint64_t first = j * span[h];
			uint64_t last = min(n, first + span[h]);
			uint64_t count = h ? (last - first + span[h - 1] - 1) / span[h - 1] : last - first;

			put<uint8_t>(buf, h ? 0 : 1);
			put<uint8_t>(buf, 0);
			put<uint16_t>(buf, (uint16_t)count);
			for (uint64_t c = 0; c < count; ++c) {
				if (h) {
					put_key(chroms[first + c * span[h - 1]].name);
					put<uint64_t>(buf, level_offset[h - 1] + (j * block_size + c) * node_size);
				} else {
					put_key(chroms[first + c].name);
					put<uint32_t>(buf, (uint32_t)(first + c));
					put<uint32_t>(buf, (uint32_t)chroms[first + c].size);
				}
			}
			buf.insert(buf.end(), (block_size - count) * (key_size + 8), '\0');
		}
		file.write(buf);
		buf.clear();
	}
}

// A zoom record being built: the bases of the records in one bin of the level's size
struct ZoomRecord {
	uint32_t chromid;
	uint64_t bin;
	uint32_t start;
	uint32_t end;
	uint64_t count;
	double   min, max, sum, sumsq;
};

// One zoom level: records are added in file order, as the records of the data (add_item) or as
// the records of the level below (add_record), and are laid out as zoom blocks of ZOOM_BLOCK_ITEMS
// records. A record that is complete is also handed to the next level.
class ZoomLevel {
public:
	uint32_t   reduction;
	uint64_t   num_records{0};
	BlockSink  sink;
	ZoomLevel *next{NULL};

	ZoomLevel(uint32_t _reduction, File *file, unsigned num_threads) : reduction(_reduction), sink(file, num_threads) {}

	void add_item(uint32_t chromid, uint32_t start, uint32_t end, float value)
	{
		double v = value;

		for (uint64_t bin = start / reduction; (uint64_t)bin * reduction < end; ++bin) {
			uint32_t s = (uint32_t)max<uint64_t>(start, bin * reduction);
			uint32_t e = (uint32_t)min<uint64_t>(end, (bin + 1) * reduction);
			double w = e - s;

			add(ZoomRecord{chromid, bin, s, e, (uint64_t)(e - s), v, v, v * w, v * v * w});
		}
	}

	void add_record(const ZoomRecord &rec)
	{
		ZoomRecord r = rec;

		r.bin = rec.start / reduction;
		add(r);
	}

	// Completes the open record and lays out the remaining records, here and in the levels above
	void finish()
	{
		if (m_open) {
			complete(m_rec);
			m_open = false;
		}
		if (!m_records.empty())
			emit_block();
		sink.flush();
		if (next)
			next->finish();
	}

private:
	bool               m_open{false};
	ZoomRecord         m_rec;
	vector<ZoomRecord> m_records;
	vector<char>       m_raw;

	void add(const ZoomRecord &r)
	{
		if (m_open && (r.chromid != m_rec.chromid || r.bin != m_rec.bin)) {
			complete(m_rec);
			m_open = false;
		}
		if (!m_open) {
			m_rec = r;
			m_open = true;
		} else {
			m_rec.end = r.end;
			m_rec.count += r.count;
			m_rec.min = min(m_rec.min, r.min);
			m_rec.max = max(m_rec.max, r.max);
			m_rec.sum += r.sum;
			m_rec.sumsq += r.sumsq;
		}
	}

	void complete(const ZoomRecord &r)
	{
		++num_records;
		m_records.push_back(r);
		if (m_records.size() >= ZOOM_BLOCK_ITEMS)
			emit_block();
		if (next)
			next->add_record(r);
	}

	void emit_block()
	{
		m_raw.clear();
		m_raw.reserve(32 * m_records.size());
		for (const ZoomRecord &r : m_records) {
			put<uint32_t>(m_raw, r.chromid);
			put<uint32_t>(m_raw, r.start);
			put<uint32_t>(m_raw, r.end);
			put<uint32_t>(m_raw, (uint32_t)r.count);
			put<float>(m_raw, (float)r.min);
			put<float>(m_raw, (float)r.max);
			put<float>(m_raw, (float)r.sum);
			put<float>(m_raw, (float)r.sumsq);
		}

		Leaf leaf{m_records.front().chromid, m_records.front().start, m_records.back().chromid, m_records.back().end, 0, 0};

		sink.add(m_raw, leaf);
		m_records.clear();
	}
};

// libBigWig's choice of zoom level sizes (bwWrite.c, makeZoomLevels): 4 times 4 times the mean
// record width, at least 10, capped at the longest chromosome; then 4 times the previous one while
// within the longest chromosome
vector<uint32_t> zoom_reductions(uint64_t valid_bases, uint64_t num_items, uint64_t max_chrom_size)
{
	vector<uint32_t> res;
	uint64_t zoom = max<uint64_t>(10, 16 * (valid_bases / num_items));

	zoom = min(zoom, max_chrom_size);
	while (res.size() < BigWigWriter::MAX_ZOOM_LEVELS && zoom <= max_chrom_size && zoom <= 0xFFFFFFFFULL) {
		res.push_back((uint32_t)zoom);
		zoom *= 4;
	}
	return res;
}

// Hands the records of data block leaf, read back from the file, to level
void read_block(File &file, const Leaf &leaf, ZoomLevel &level, vector<char> &cbuf, vector<char> &ubuf)
{
	cbuf.resize(leaf.size);
	file.seek(leaf.offset);
	file.read(cbuf.data(), leaf.size);

	uLongf size = BUF_SIZE;
	ubuf.resize(BUF_SIZE);
	if (uncompress((Bytef *)ubuf.data(), &size, (const Bytef *)cbuf.data(), leaf.size) != Z_OK || size < 24)
		TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "Failed to read back a data block of the bigWig file");

	const char *p = ubuf.data();
	uint32_t chromid = get<uint32_t>(p);
	uint32_t start = get<uint32_t>(p + 4);
	uint32_t step = get<uint32_t>(p + 12);
	uint32_t span = get<uint32_t>(p + 16);
	uint8_t type = p[20];
	uint16_t count = get<uint16_t>(p + 22);

	if (size != 24 + (uint64_t)count * (type == BEDGRAPH ? 12 : 4))
		TGLError<BigWigWriter>(BigWigWriter::BAD_STATE, "Failed to read back a data block of the bigWig file");

	p += 24;
	for (uint16_t i = 0; i < count; ++i) {
		if (type == BEDGRAPH) {
			level.add_item(chromid, get<uint32_t>(p), get<uint32_t>(p + 4), get<float>(p + 8));
			p += 12;
		} else {
			uint32_t s = start + i * step;

			level.add_item(chromid, s, s + span, get<float>(p));
			p += 4;
		}
	}
}

}

void BigWigWriter::create(const string &filename, const vector<Chrom> &chroms, vector<char> &state)
{
	if (chroms.empty())
		TGLError<BigWigWriter>(BAD_ARG, "A bigWig file needs at least one chromosome");

	for (uint64_t i = 0; i < chroms.size(); ++i) {
		const Chrom &chrom = chroms[i];

		if (chrom.name.empty() || chrom.name.find('\0') != string::npos)
			TGLError<BigWigWriter>(BAD_ARG, "Invalid chromosome name \"%s\" for a bigWig file", chrom.name.c_str());
		if (chrom.size < 1 || chrom.size > 0xFFFFFFFFULL)
			TGLError<BigWigWriter>(BAD_ARG, "Chromosome %s is %llu bases long; a bigWig file holds chromosomes of 1 to 4294967295 bases",
			                       chrom.name.c_str(), (unsigned long long)chrom.size);
		if (i && !(chroms[i - 1].name < chrom.name))
			TGLError<BigWigWriter>(BAD_ARG, "The chromosomes of a bigWig file must be sorted by name and unique: %s comes after %s",
			                       chrom.name.c_str(), chroms[i - 1].name.c_str());
	}

	File file;
	vector<char> zeros(CHROM_TREE_OFFSET, 0);

	file.open(filename, "wb");
	file.write(zeros);
	write_chrom_tree(file, chroms);

	State s;

	memset(&s.h, 0, sizeof(s.h));
	s.h.magic = STATE_MAGIC;
	s.h.data_offset = file.tell();
	s.h.min = INFINITY;
	s.h.max = -INFINITY;
	s.h.last_chromid = -1;
	for (const Chrom &chrom : chroms) {
		s.sizes.push_back((uint32_t)chrom.size);
		s.names.push_back(chrom.name);
	}

	uint64_t num_blocks = 0;

	file.write(&num_blocks, sizeof(num_blocks));   // set by finish()
	s.h.file_end = file.tell();
	file.close();
	s.serialize(state);
}

void BigWigWriter::append(const string &filename, vector<char> &state, uint64_t n, const int64_t *chromids,
                          const int64_t *starts, const int64_t *ends, const double *values, unsigned num_threads)
{
	State s;
	File file;

	s.parse(state);
	reopen(file, filename, s);

	BlockSink sink(&file, num_threads);
	DataLayout layout(sink);
	StateHeader &h = s.h;
	int64_t num_chroms = s.sizes.size();

	for (uint64_t i = 0; i < n; ++i) {
		if (!(i % ITEMS_PER_INTERRUPT_CHECK))
			rdb::check_interrupt();

		int64_t chromid = chromids[i];
		int64_t start = starts[i];
		int64_t end = ends[i];

		if (chromid < 0 || chromid >= num_chroms)
			TGLError<BigWigWriter>(BAD_RECORD, "Invalid chromosome id %lld of a bigWig record", (long long)chromid);

		const char *chrom = s.names[chromid].c_str();

		if (start < 0 || end < start || end > (int64_t)s.sizes[chromid])
			TGLError<BigWigWriter>(BAD_RECORD, "Invalid bigWig record (%s, %lld, %lld): chromosome %s is %llu bases long",
			                       chrom, (long long)start, (long long)end, chrom, (unsigned long long)s.sizes[chromid]);

		if (std::isnan(values[i]) || start == end)
			continue;

		if (chromid < h.last_chromid)
			TGLError<BigWigWriter>(BAD_RECORD, "The records of a bigWig file must be sorted by chromosome name: %s comes after %s",
			                       chrom, s.names[h.last_chromid].c_str());
		if (chromid == h.last_chromid && start < h.last_end)
			TGLError<BigWigWriter>(BAD_RECORD, "The records of a bigWig file must be sorted and must not overlap: (%s, %lld, %lld) starts before %lld, the end of the previous record",
			                       chrom, (long long)start, (long long)end, (long long)h.last_end);

		float value = (float)values[i];
		double v = value;
		double w = end - start;

		h.last_chromid = chromid;
		h.last_end = end;
		++h.num_items;
		h.valid_bases += end - start;
		h.min = min(h.min, v);
		h.max = max(h.max, v);
		h.sum += v * w;
		h.sumsq += v * v * w;
		layout.add((uint32_t)chromid, (uint32_t)start, (uint32_t)end, value);
	}

	layout.flush(true);
	sink.flush();
	s.leaves.insert(s.leaves.end(), sink.leaves.begin(), sink.leaves.end());
	h.file_end = file.tell();
	file.close();
	s.serialize(state);
}

uint64_t BigWigWriter::finish(const string &filename, const vector<char> &state, unsigned num_threads)
{
	State s;
	File file;

	s.parse(state);
	reopen(file, filename, s);

	const StateHeader &h = s.h;
	uint64_t index_offset = file.tell();

	write_rtree(file, s.leaves, ITEMS_PER_SLOT);

	// Zoom levels, from the data blocks read back: the first is written as it is made, the others
	// are kept in memory (each holds at most a quarter of the records of the one below) and written
	// after it
	vector<uint32_t> reductions;
	uint64_t max_chrom_size = *max_element(s.sizes.begin(), s.sizes.end());

	if (h.num_items)
		reductions = zoom_reductions(h.valid_bases, h.num_items, max_chrom_size);

	struct ZoomHeader {
		uint32_t reduction;
		uint64_t data_offset;
		uint64_t index_offset;
	};
	vector<ZoomHeader> zoom_headers;

	if (!reductions.empty()) {
		vector<ZoomLevel *> levels;
		File reader;
		vector<char> cbuf, ubuf;
		uint64_t first_data_offset = file.tell();
		uint32_t count = 0;

		file.write(&count, sizeof(count));   // the number of records, set below

		try {
			for (uint64_t i = 0; i < reductions.size(); ++i) {
				levels.push_back(new ZoomLevel(reductions[i], i ? NULL : &file, num_threads));
				if (i)
					levels[i - 1]->next = levels[i];
			}

			reader.open(filename, "rb");
			for (const Leaf &leaf : s.leaves) {
				rdb::check_interrupt();
				read_block(reader, leaf, *levels[0], cbuf, ubuf);
			}
			reader.close();
			levels[0]->finish();

			// the first level, already written
			uint64_t end = file.tell();

			count = (uint32_t)min<uint64_t>(levels[0]->num_records, 0xFFFFFFFFULL);
			file.seek(first_data_offset);
			file.write(&count, sizeof(count));
			file.seek(end);
			zoom_headers.push_back(ZoomHeader{levels[0]->reduction, first_data_offset, end});
			write_rtree(file, levels[0]->sink.leaves, ZOOM_ITEMS);

			uint64_t fewest = levels[0]->num_records;

			for (uint64_t i = 1; i < levels.size(); ++i) {
				ZoomLevel &level = *levels[i];

				if (level.num_records >= fewest)
					continue;
				fewest = level.num_records;

				uint64_t data_offset = file.tell();
				uint64_t blob_offset = data_offset + sizeof(count);

				count = (uint32_t)min<uint64_t>(level.num_records, 0xFFFFFFFFULL);
				file.write(&count, sizeof(count));
				file.write(level.sink.blob);
				for (Leaf &leaf : level.sink.leaves)
					leaf.offset += blob_offset;
				zoom_headers.push_back(ZoomHeader{level.reduction, data_offset, file.tell()});
				write_rtree(file, level.sink.leaves, ZOOM_ITEMS);
			}
		} catch (...) {
			for (ZoomLevel *level : levels)
				delete level;
			throw;
		}
		for (ZoomLevel *level : levels)
			delete level;
	}

	uint32_t magic = BIGWIG_MAGIC;

	file.write(&magic, sizeof(magic));

	// header
	vector<char> buf;

	put<uint32_t>(buf, BIGWIG_MAGIC);
	put<uint16_t>(buf, 4);                        // version
	put<uint16_t>(buf, (uint16_t)zoom_headers.size());
	put<uint64_t>(buf, CHROM_TREE_OFFSET);
	put<uint64_t>(buf, h.data_offset);
	put<uint64_t>(buf, index_offset);
	put<uint16_t>(buf, 0);                        // fieldCount
	put<uint16_t>(buf, 0);                        // definedFieldCount
	put<uint64_t>(buf, 0);                        // autoSqlOffset
	put<uint64_t>(buf, SUMMARY_OFFSET);
	put<uint32_t>(buf, BUF_SIZE);                 // uncompressBufSize
	put<uint64_t>(buf, 0);                        // extensionOffset
	for (const ZoomHeader &zh : zoom_headers) {
		put<uint32_t>(buf, zh.reduction);
		put<uint32_t>(buf, 0);
		put<uint64_t>(buf, zh.data_offset);
		put<uint64_t>(buf, zh.index_offset);
	}
	file.seek(0);
	file.write(buf);

	// total summary: validCount, minVal, maxVal, sumData, sumSquares
	buf.clear();
	put<uint64_t>(buf, h.valid_bases);
	put<double>(buf, h.num_items ? h.min : 0);
	put<double>(buf, h.num_items ? h.max : 0);
	put<double>(buf, h.sum);
	put<double>(buf, h.sumsq);
	file.seek(SUMMARY_OFFSET);
	file.write(buf);

	uint64_t num_blocks = s.leaves.size();

	file.seek(h.data_offset);
	file.write(&num_blocks, sizeof(num_blocks));
	file.close();
	return h.num_items;
}
