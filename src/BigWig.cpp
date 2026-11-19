#include <errno.h>
#include <strings.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <zlib.h>

#include "BigWig.h"

// bigWig and bigBed layout: https://genome.ucsc.edu/goldenPath/help/bigWig.html, bigBed.html and
// Kent et al. 2010, Bioinformatics 26:2204, supplement. All integers are little-endian.
static const uint32_t BIGWIG_MAGIC = 0x888FFC26;
static const uint32_t BIGBED_MAGIC = 0x8789F2EB;
static const uint32_t CHROM_TREE_MAGIC = 0x78CA8C91;
static const uint32_t RTREE_MAGIC = 0x2468ACE0;
static const int MAX_TREE_DEPTH = 64;   // a corrupt child offset must not recurse forever
static const uint64_t MAX_BLOCK_SIZE = 24 + 65535 * 12;   // section header + itemCount (16 bits) x largest item

template <class T> static T get(const char *p)
{
	T v;
	memcpy(&v, p, sizeof(T));
	return v;
}

// whether [start, end) overlaps one of the intervals of query (sorted, non-overlapping)
static bool overlaps(const GIntervals &query, int64_t start, int64_t end)
{
	GIntervals::const_iterator iq = lower_bound(query.begin(), query.end(), start,
	                                            [](const GInterval &q, int64_t s) { return q.end <= s; });
	return iq != query.end() && iq->start < end;
}

// autoSql text split into words, quoted strings and single punctuation characters; # starts a comment
static bool autosql_tokens(const string &text, vector<string> &tokens)
{
	for (size_t i = 0; i < text.size(); ) {
		unsigned char c = text[i];

		if (isspace(c))
			++i;
		else if (c == '#') {
			while (i < text.size() && text[i] != '\n')
				++i;
		} else if (c == '"') {
			size_t j = text.find('"', i + 1);

			if (j == string::npos)
				return false;
			tokens.push_back(text.substr(i, j + 1 - i));
			i = j + 1;
		} else if (isalnum(c) || c == '_') {
			size_t j = i;

			while (j < text.size() && (isalnum((unsigned char)text[j]) || text[j] == '_'))
				++j;
			tokens.push_back(text.substr(i, j - i));
			i = j;
		} else {
			tokens.push_back(string(1, c));
			++i;
		}
	}
	return true;
}

// The columns of the first object an autoSql text defines (kent src/lib/asParse.c):
//   table|simple|object NAME "comment" ( COLUMN... )
//   COLUMN: TYPE [OBJECT-NAME] [ [SIZE] | (VALUES) ] NAME [primary|unique|index [SIZE]] [auto] ; "comment"
// A column is numeric if its type is a number and it is not an array.
static bool parse_autosql(const string &text, vector<BigWig::Field> &fields)
{
	static const char *NUMERIC_TYPES[] = { "int", "uint", "short", "ushort", "byte", "ubyte", "float", "double", "bigint" };
	vector<string> tok;

	if (!autosql_tokens(text, tok))
		return false;

	auto is_word = [&](size_t k) { return k < tok.size() && (isalnum((unsigned char)tok[k][0]) || tok[k][0] == '_'); };
	size_t i = 2;

	if (!is_word(0) || !is_word(1))
		return false;
	if (i < tok.size() && tok[i][0] == '"')
		++i;
	if (i >= tok.size() || tok[i] != "(")
		return false;

	for (++i; i < tok.size() && tok[i] != ")"; ) {
		if (!is_word(i))
			return false;

		const string &type = tok[i++];
		bool numeric = find(begin(NUMERIC_TYPES), end(NUMERIC_TYPES), type) != end(NUMERIC_TYPES);

		if (type == "object" || type == "simple" || type == "table")
			++i;   // the name of the embedded object

		if (i < tok.size() && (tok[i] == "[" || tok[i] == "(")) {   // array size, or enum / set values
			const char *close = tok[i] == "[" ? "]" : ")";

			numeric = false;
			while (i < tok.size() && tok[i] != close)
				++i;
			++i;
		}

		if (!is_word(i))
			return false;
		fields.push_back(BigWig::Field{tok[i++], numeric});

		while (i < tok.size() && tok[i] != ";")   // index, auto
			++i;
		if (++i < tok.size() && tok[i][0] == '"')
			++i;
	}
	return i < tok.size();
}

bool BigWig::is_bigwig(const char *filename)
{
	size_t len = strlen(filename);

	if ((len > 3 && !strcasecmp(filename + len - 3, ".bw")) || (len > 7 && !strcasecmp(filename + len - 7, ".bigwig")))
		return true;

	FILE *fp = fopen(filename, "rb");
	uint32_t magic = 0;
	bool res = fp && fread(&magic, sizeof(magic), 1, fp) == 1 &&
		(magic == BIGWIG_MAGIC || magic == __builtin_bswap32(BIGWIG_MAGIC));

	if (fp)
		fclose(fp);
	return res;
}

void BigWig::read_at(uint64_t offset, void *buf, uint64_t size)
{
	if (m_bfile.seek(offset, SEEK_SET) || m_bfile.read(buf, size) != size)
		TGLError<BigWig>(FILE_ERROR, "Failed to read %s file %s (%llu bytes at offset %llu): %s",
		                 type_name(), m_bfile.file_name().c_str(), (unsigned long long)size, (unsigned long long)offset,
		                 m_bfile.error() ? strerror(errno) : "the file is truncated or corrupt");
}

void BigWig::init(const GenomeChromKey &chromkey, const string &filename, bool bigbed)
{
	m_chromkey = &chromkey;
	m_bigbed = bigbed;
	m_blocks.clear();
	m_bw_names.clear();
	m_fields.clear();
	m_unknown_chroms.clear();
	m_last_chromid = -1;
	m_num_clipped = m_num_dropped = 0;
	m_first_clipped.clear();
	m_first_clipped_offset = 0;
	m_first_clipped_item = 0;

	if (m_bfile.open(filename.c_str(), "rb"))
		TGLError<BigWig>(FILE_ERROR, "Failed to open %s file %s: %s", type_name(), filename.c_str(), strerror(errno));

	const char *fname = m_bfile.file_name().c_str();
	char hdr[64];

	read_at(0, hdr, sizeof(hdr));

	// ponytail: big-endian files are refused; add byte swapping if one ever turns up
	uint32_t magic = bigbed ? BIGBED_MAGIC : BIGWIG_MAGIC;

	if (get<uint32_t>(hdr) == __builtin_bswap32(magic))
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s is big-endian, which is not supported. Convert it with UCSC's %s and import that.",
		                 type_name(), fname, bigbed ? "bigBedToBed" : "bigWigToBedGraph");
	if (get<uint32_t>(hdr) != magic)
		TGLError<BigWig>(INVALID_FORMAT, "%s is not a %s file", fname, type_name());

	uint64_t chrom_tree_offset = get<uint64_t>(hdr + 8);
	uint64_t index_offset = get<uint64_t>(hdr + 24);

	m_data_offset = get<uint64_t>(hdr + 16);
	m_uncompress_bufsize = get<uint32_t>(hdr + 52);

	// chromosome B+ tree: magic, blockSize, keySize, valSize, itemCount, reserved
	char ct[32];
	read_at(chrom_tree_offset, ct, sizeof(ct));
	if (get<uint32_t>(ct) != CHROM_TREE_MAGIC)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid chromosome tree", type_name(), fname);

	uint32_t keysize = get<uint32_t>(ct + 8);
	if (!keysize || keysize > 1024)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid chromosome tree", type_name(), fname);

	// Each visited tree node is charged its size on disk. The nodes of a valid file are disjoint, so
	// the total stays within the file size; a crafted tree whose nodes share children would otherwise
	// take exponential time and memory.
	m_tree_budget = m_bfile.file_size();

	// Names are sized by the chromosome ids actually found, never by a header count: ids are 0..n-1.
	vector<pair<uint32_t, string> > chroms;

	read_chrom_tree(chrom_tree_offset + sizeof(ct), keysize, 0, chroms);

	uint64_t num_chroms = chroms.size();

	vector<bool> seen(num_chroms, false);

	m_bw_names.resize(num_chroms);
	m_blocks.resize(num_chroms);
	for (vector<pair<uint32_t, string> >::const_iterator ichrom = chroms.begin(); ichrom != chroms.end(); ++ichrom) {
		if (ichrom->first >= num_chroms || seen[ichrom->first])
			TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid chromosome tree", type_name(), fname);
		seen[ichrom->first] = true;
		m_bw_names[ichrom->first] = ichrom->second;
	}

	// R-tree index: magic, blockSize, itemCount, 4 x bounds, endFileOffset, itemsPerSlot, reserved
	char rt[48];
	read_at(index_offset, rt, sizeof(rt));
	if (get<uint32_t>(rt) != RTREE_MAGIC)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid data index", type_name(), fname);

	read_rtree(index_offset + sizeof(rt), 0);

	m_db2bw.assign(chromkey.get_num_chroms(), -1);
	m_bw2db.resize(num_chroms);
	for (uint32_t bwid = 0; bwid < num_chroms; ++bwid)
		m_bw2db[bwid] = chromkey.find_chromid(m_bw_names[bwid]);

	// A bigBed block can hold records of several chromosomes: they are matched as the records are
	// read, and the header's record count is checked once they all are.
	if (bigbed) {
		read_fields(get<uint64_t>(hdr + 36), get<uint16_t>(hdr + 32), get<uint16_t>(hdr + 34));
		return;
	}

	// every data block the header counts must be reachable from the index
	uint64_t num_blocks = 0;
	char data_count[8];

	for (vector<vector<Block> >::const_iterator iblocks = m_blocks.begin(); iblocks != m_blocks.end(); ++iblocks)
		num_blocks += iblocks->size();
	read_at(m_data_offset, data_count, sizeof(data_count));
	if (get<uint64_t>(data_count) != num_blocks)
		TGLError<BigWig>(INCOMPLETE, "bigWig file %s: its index lists %llu data blocks, but its header counts %llu",
		                 fname, (unsigned long long)num_blocks, (unsigned long long)get<uint64_t>(data_count));

	// A chromosome the tree lists without data (libBigWig writers list every chromosome they
	// were given) has nothing to import, so it is neither matched nor reported as skipped.
	for (uint32_t bwid = 0; bwid < num_chroms; ++bwid) {
		if (m_blocks[bwid].empty())
			continue;

		int chromid = m_bw2db[bwid];

		if (chromid < 0)
			m_unknown_chroms.record(m_bw_names[bwid]);
		else if (m_db2bw[chromid] >= 0)
			TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s has data for both %s and %s, which are the same database chromosome %s",
			                 fname, m_bw_names[m_db2bw[chromid]].c_str(), m_bw_names[bwid].c_str(), chromkey.id2chrom(chromid).c_str());
		else
			m_db2bw[chromid] = bwid;
	}
}

void BigWig::read_chrom_tree(uint64_t offset, uint32_t keysize, int depth, vector<pair<uint32_t, string> > &chroms)
{
	char node[4];   // isLeaf, reserved, count

	if (depth > MAX_TREE_DEPTH)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid chromosome tree", type_name(), m_bfile.file_name().c_str());

	read_at(offset, node, sizeof(node));

	uint16_t count = get<uint16_t>(node + 2);
	uint64_t item_size = keysize + 8;   // leaf: key, chromId, chromSize; non-leaf: key, child offset

	if (m_tree_budget < sizeof(node) + count * item_size)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid chromosome tree", type_name(), m_bfile.file_name().c_str());
	m_tree_budget -= sizeof(node) + count * item_size;

	vector<uint64_t> children;
	{
		vector<char> items(count * item_size);

		read_at(offset + sizeof(node), items.data(), items.size());
		for (uint16_t i = 0; i < count; ++i) {
			const char *item = &items[i * item_size];

			if (node[0])
				chroms.push_back(make_pair(get<uint32_t>(item + keysize), string(item, strnlen(item, keysize))));
			else
				children.push_back(get<uint64_t>(item + keysize));
		}
	}   // items is freed before recursing, so a crafted deep tree cannot hold one buffer per level

	for (vector<uint64_t>::const_iterator ichild = children.begin(); ichild != children.end(); ++ichild)
		read_chrom_tree(*ichild, keysize, depth + 1, chroms);
}

void BigWig::read_rtree(uint64_t offset, int depth)
{
	char node[4];   // isLeaf, reserved, count

	if (depth > MAX_TREE_DEPTH)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid data index", type_name(), m_bfile.file_name().c_str());

	read_at(offset, node, sizeof(node));

	uint16_t count = get<uint16_t>(node + 2);
	uint64_t item_size = node[0] ? 32 : 24;   // 4 x bounds, then leaf: data offset, data size; non-leaf: child offset

	if (m_tree_budget < sizeof(node) + count * item_size)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid data index", type_name(), m_bfile.file_name().c_str());
	m_tree_budget -= sizeof(node) + count * item_size;

	vector<uint64_t> children;
	{
		vector<char> items(count * item_size);

		read_at(offset + sizeof(node), items.data(), items.size());
		for (uint16_t i = 0; i < count; ++i) {
			const char *item = &items[i * item_size];

			if (node[0]) {
				// bounds: startChromIx, startBase, endChromIx, endBase. A bigWig data block holds a single
				// chromosome; a bigBed block is listed under its first one.
				uint32_t bwid = get<uint32_t>(item);
				Block block{get<uint64_t>(item + 16), get<uint64_t>(item + 24), get<uint32_t>(item + 4),
				            get<uint32_t>(item + 8) == bwid ? (int64_t)get<uint32_t>(item + 12) : INT64_MAX};
				uint64_t file_size = m_bfile.file_size();

				if (bwid >= m_blocks.size() || block.size > file_size || block.offset > file_size - block.size)
					TGLError<BigWig>(INVALID_FORMAT, "%s file %s: invalid data index", type_name(), m_bfile.file_name().c_str());
				m_blocks[bwid].push_back(block);
			} else
				children.push_back(get<uint64_t>(item + 16));
		}
	}   // items is freed before recursing, as in read_chrom_tree

	for (vector<uint64_t>::const_iterator ichild = children.begin(); ichild != children.end(); ++ichild)
		read_rtree(*ichild, depth + 1);
}

// the BED standard field names (kent src/lib/basicBed.c, bedAsDef) and whether each is numeric
static const BigWig::Field BED_FIELDS[] = {
	{"chrom", false}, {"chromStart", true}, {"chromEnd", true}, {"name", false}, {"score", true}, {"strand", false},
	{"thickStart", true}, {"thickEnd", true}, {"reserved", true}, {"blockCount", true}, {"blockSizes", false},
	{"chromStarts", false}, {"expCount", true}, {"expIds", false}, {"expScores", false}
};

void BigWig::read_fields(uint64_t autosql_offset, uint16_t field_count, uint16_t defined_field_count)
{
	const char *fname = m_bfile.file_name().c_str();

	if (field_count < 3 || defined_field_count < 3 || defined_field_count > field_count)
		TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: invalid field counts in its header", fname);

	if (autosql_offset) {
		// a NUL-terminated string
		string text;
		uint64_t file_size = m_bfile.file_size();

		for (uint64_t pos = autosql_offset; text.find('\0') == string::npos; ) {
			char buf[4096];
			uint64_t size = min((uint64_t)sizeof(buf), pos < file_size ? file_size - pos : 0);

			if (!size)
				TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: its autoSql is not terminated", fname);
			read_at(pos, buf, size);
			text.append(buf, size);
			pos += size;
		}
		text.resize(text.find('\0'));

		if (!parse_autosql(text, m_fields))
			TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: cannot parse its autoSql:\n%s", fname, text.c_str());
		if (m_fields.size() != field_count)
			TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: its autoSql declares %llu fields, but its header %u",
			                 fname, (unsigned long long)m_fields.size(), (unsigned)field_count);
	} else {
		if (field_count > defined_field_count)
			TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s has %u fields past its %u BED standard ones, and no autoSql to name them",
			                 fname, (unsigned)(field_count - defined_field_count), (unsigned)defined_field_count);
		if (defined_field_count > sizeof(BED_FIELDS) / sizeof(BED_FIELDS[0]))
			TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: invalid field counts in its header", fname);
		m_fields.assign(BED_FIELDS, BED_FIELDS + defined_field_count);
	}

	// The BED standard's itemRgb is "R,G,B" text, although autoSql declares it a uint
	if (defined_field_count >= 9)
		m_fields[8].numeric = false;

	m_fields.erase(m_fields.begin(), m_fields.begin() + 3);   // chrom, start and end are in the record header
}

const char *BigWig::read_block(const Block &block, uint64_t &size)
{
	m_cbuf.resize(block.size);
	read_at(block.offset, m_cbuf.data(), block.size);

	if (!m_uncompress_bufsize) {
		size = block.size;
		return m_cbuf.data();
	}

	// deflate cannot expand data more than 1032-fold, so a bigBed header's buffer size is not trusted further
	uLongf usize = min((uint64_t)m_uncompress_bufsize, m_bigbed ? block.size * 1032 + 64 : MAX_BLOCK_SIZE);
	int rc;

	m_ubuf.resize(usize);
	rc = uncompress((Bytef *)m_ubuf.data(), &usize, (const Bytef *)m_cbuf.data(), block.size);
	if (rc != Z_OK)
		TGLError<BigWig>(INVALID_FORMAT, "%s file %s: data block at offset %llu failed to decompress (zlib error %d)",
		                 type_name(), m_bfile.file_name().c_str(), (unsigned long long)block.offset, rc);
	size = usize;
	return m_ubuf.data();
}

bool BigWig::clip(uint32_t bwid, int chromid, const Block &block, uint64_t item, int64_t start, int64_t &end)
{
	int64_t chromsize = m_chromkey->get_chrom_size(chromid);

	if (end <= chromsize)
		return true;

	// As in the other importers: a record that runs past the end of the chromosome is clipped to it,
	// and dropped if it starts at or past the end. The caller reports them.
	// The chromosomes are read in database order: the first is the one of the lowest offset.
	if (m_first_clipped.empty() || block.offset < m_first_clipped_offset ||
	    (block.offset == m_first_clipped_offset && item < m_first_clipped_item)) {
		m_first_clipped_offset = block.offset;
		m_first_clipped_item = item;
		char buf[1000];
		snprintf(buf, sizeof(buf), "(%s, %lld, %lld), chromosome size %lld", m_bw_names[bwid].c_str(), (long long)start,
		         (long long)end, (long long)chromsize);
		m_first_clipped = buf;
	}
	if (start >= chromsize) {
		++m_num_dropped;
		return false;
	}
	++m_num_clipped;
	end = chromsize;
	return true;
}

bool BigWig::get_data(int chromid, uint64_t iblock, GIntervals &intervals, const GIntervals *query)
{
	intervals.clear();
	if (m_db2bw[chromid] < 0 || iblock >= m_blocks[m_db2bw[chromid]].size())
		return false;

	const char *fname = m_bfile.file_name().c_str();
	uint32_t bwid = m_db2bw[chromid];
	const Block &block = m_blocks[bwid][iblock];

	if (!iblock || chromid != m_last_chromid)
		m_last_end = 0;
	m_last_chromid = chromid;

	uint64_t size;
	const char *data = read_block(block, size);

	// section header: chromId, start, end, itemStep, itemSpan, type, reserved, itemCount
	if (size < 24)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: malformed data block at offset %llu", fname, (unsigned long long)block.offset);

	int64_t start = get<uint32_t>(data + 4);
	int64_t step = get<uint32_t>(data + 12);
	int64_t span = get<uint32_t>(data + 16);
	uint8_t type = data[20];
	uint16_t num_items = get<uint16_t>(data + 22);
	uint64_t item_size = type == 1 ? 12 : type == 2 ? 8 : type == 3 ? 4 : 0;   // bedGraph, variableStep, fixedStep

	if (get<uint32_t>(data) != bwid || !item_size || size != 24 + num_items * item_size)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: malformed data block at offset %llu", fname, (unsigned long long)block.offset);

	const char *item = data + 24;

	for (uint16_t i = 0; i < num_items; ++i, item += item_size) {
		int64_t s, e;
		float v;

		if (type == 1) {
			s = get<uint32_t>(item);
			e = get<uint32_t>(item + 4);
			v = get<float>(item + 8);
		} else if (type == 2) {
			s = get<uint32_t>(item);
			e = s + span;
			v = get<float>(item + 4);
		} else {
			s = start + i * step;
			e = s + span;
			v = get<float>(item);
		}

		if (s >= e)   // bedGraphToBigWig stores zero-length records; bigWigToWig skipped them
			continue;

		if (query && !overlaps(*query, s, e))
			continue;

		if (!clip(bwid, chromid, block, i, s, e))
			continue;

		if (s < m_last_end)
			TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: records of %s are overlapping or unsorted at %lld",
			                 fname, m_bw_names[bwid].c_str(), (long long)s);

		m_last_end = e;
		intervals.push_back(GInterval(chromid, s, e, 0, GInterval::cast2udata(v)));
	}
	return true;
}

void BigWig::get_data(int chromid, const GIntervals &query, GIntervals &intervals)
{
	GIntervals data;

	intervals.clear();
	if (m_db2bw[chromid] < 0)
		return;

	const vector<Block> &blocks = m_blocks[m_db2bw[chromid]];

	for (uint64_t iblock = 0; iblock < blocks.size(); ++iblock) {
		if (overlaps(query, blocks[iblock].start, blocks[iblock].end)) {
			get_data(chromid, iblock, data, &query);
			intervals.insert(intervals.end(), data.begin(), data.end());
		}
	}
}

void BigWig::read_bed(BedRecords &records)
{
	const char *fname = m_bfile.file_name().c_str();
	uint64_t num_fields = m_fields.size();
	uint64_t num_records = 0;
	char data_count[8];

	records.intervals.clear();
	records.values.assign(num_fields, vector<double>());
	records.texts.assign(num_fields, vector<string>());

	for (uint32_t ibwid = 0; ibwid < m_blocks.size(); ++ibwid) {
		for (vector<Block>::const_iterator iblock = m_blocks[ibwid].begin(); iblock != m_blocks[ibwid].end(); ++iblock) {
			uint64_t size;
			const char *data = read_block(*iblock, size);
			const char *data_end = data + size;

			// a record: chromId, start, end, then the other fields, tab separated, up to a NUL
			for (uint64_t item = 0; data < data_end; ++item, ++num_records) {
				const char *rest = data + 12;
				const char *rest_end = rest < data_end ? (const char *)memchr(rest, '\0', data_end - rest) : NULL;

				if (!rest_end)
					TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: malformed data block at offset %llu", fname, (unsigned long long)iblock->offset);

				uint32_t bwid = get<uint32_t>(data);
				int64_t s = get<uint32_t>(data + 4);
				int64_t e = get<uint32_t>(data + 8);

				data = rest_end + 1;
				if (bwid >= m_bw2db.size() || s > e)
					TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: malformed data block at offset %llu", fname, (unsigned long long)iblock->offset);

				int chromid = m_bw2db[bwid];

				if (chromid < 0) {
					m_unknown_chroms.record(m_bw_names[bwid]);
					continue;
				}

				// fields 4+, tab separated
				if (num_fields ? (uint64_t)count(rest, rest_end, '\t') != num_fields - 1 : rest != rest_end)
					TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: the record (%s, %lld, %lld) does not have the %llu fields the file declares",
					                 fname, m_bw_names[bwid].c_str(), (long long)s, (long long)e, (unsigned long long)num_fields + 3);

				int64_t file_end = e;

				m_db2bw[chromid] = bwid;
				if (!clip(bwid, chromid, *iblock, item, s, e))
					continue;

				const char *field = rest;

				for (uint64_t ifield = 0; ifield < num_fields; ++ifield) {
					const char *field_end = find(field, rest_end, '\t');

					if (m_fields[ifield].numeric) {
						string str(field, field_end);
						char *endp;
						double v = strtod(str.c_str(), &endp);

						if (str.empty() || *endp)
							TGLError<BigWig>(INVALID_FORMAT, "bigBed file %s: field %s of the record (%s, %lld, %lld) is not a number: \"%s\"",
							                 fname, m_fields[ifield].name.c_str(), m_bw_names[bwid].c_str(), (long long)s,
							                 (long long)file_end, str.c_str());
						records.values[ifield].push_back(v);
					} else
						records.texts[ifield].push_back(string(field, field_end));
					field = field_end + 1;
				}

				records.intervals.push_back(GInterval(chromid, s, e, 0));
			}
		}
	}

	// every record the header counts must have been read
	read_at(m_data_offset, data_count, sizeof(data_count));
	if (get<uint64_t>(data_count) != num_records)
		TGLError<BigWig>(INCOMPLETE, "bigBed file %s: %llu records were read, but its header counts %llu",
		                 fname, (unsigned long long)num_records, (unsigned long long)get<uint64_t>(data_count));
}

uint64_t BigWig::get_num_matched_chroms() const
{
	uint64_t num_matched = 0;

	for (vector<int>::const_iterator ibwid = m_db2bw.begin(); ibwid != m_db2bw.end(); ++ibwid) {
		if (*ibwid >= 0)
			++num_matched;
	}
	return num_matched;
}
