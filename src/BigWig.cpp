#include <errno.h>
#include <strings.h>
#include <cstring>
#include <zlib.h>

#include "BigWig.h"

// bigWig layout: https://genome.ucsc.edu/goldenPath/help/bigWig.html and Kent et al. 2010,
// Bioinformatics 26:2204, supplement. All integers are little-endian.
static const uint32_t BIGWIG_MAGIC = 0x888FFC26;
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
		TGLError<BigWig>(FILE_ERROR, "Failed to read bigWig file %s (%llu bytes at offset %llu): %s",
		                 m_bfile.file_name().c_str(), (unsigned long long)size, (unsigned long long)offset,
		                 m_bfile.error() ? strerror(errno) : "the file is truncated or corrupt");
}

void BigWig::init(const GenomeChromKey &chromkey, const string &filename)
{
	m_chromkey = &chromkey;
	m_blocks.clear();
	m_bw_names.clear();
	m_unknown_chroms.clear();
	m_num_clipped = m_num_dropped = 0;
	m_first_clipped.clear();
	m_first_clipped_offset = 0;
	m_first_clipped_item = 0;

	if (m_bfile.open(filename.c_str(), "rb"))
		TGLError<BigWig>(FILE_ERROR, "Failed to open bigWig file %s: %s", filename.c_str(), strerror(errno));

	const char *fname = m_bfile.file_name().c_str();
	char hdr[64];

	read_at(0, hdr, sizeof(hdr));

	// ponytail: big-endian bigWigs are refused; add byte swapping if one ever turns up
	if (get<uint32_t>(hdr) != BIGWIG_MAGIC)
		TGLError<BigWig>(INVALID_FORMAT, get<uint32_t>(hdr) == __builtin_bswap32(BIGWIG_MAGIC) ?
		                 "bigWig file %s is big-endian, which is not supported. Convert it with UCSC's bigWigToBedGraph and import that." :
		                 "%s is not a bigWig file", fname);

	uint64_t chrom_tree_offset = get<uint64_t>(hdr + 8);
	uint64_t data_offset = get<uint64_t>(hdr + 16);
	uint64_t index_offset = get<uint64_t>(hdr + 24);

	m_uncompress_bufsize = get<uint32_t>(hdr + 52);

	// chromosome B+ tree: magic, blockSize, keySize, valSize, itemCount, reserved
	char ct[32];
	read_at(chrom_tree_offset, ct, sizeof(ct));
	if (get<uint32_t>(ct) != CHROM_TREE_MAGIC)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid chromosome tree", fname);

	uint32_t keysize = get<uint32_t>(ct + 8);
	if (!keysize || keysize > 1024)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid chromosome tree", fname);

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
			TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid chromosome tree", fname);
		seen[ichrom->first] = true;
		m_bw_names[ichrom->first] = ichrom->second;
	}

	// R-tree index: magic, blockSize, itemCount, 4 x bounds, endFileOffset, itemsPerSlot, reserved
	char rt[48];
	read_at(index_offset, rt, sizeof(rt));
	if (get<uint32_t>(rt) != RTREE_MAGIC)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid data index", fname);

	read_rtree(index_offset + sizeof(rt), 0);

	// every data block the header counts must be reachable from the index
	uint64_t num_blocks = 0;
	char data_count[8];

	for (vector<vector<Block> >::const_iterator iblocks = m_blocks.begin(); iblocks != m_blocks.end(); ++iblocks)
		num_blocks += iblocks->size();
	read_at(data_offset, data_count, sizeof(data_count));
	if (get<uint64_t>(data_count) != num_blocks)
		TGLError<BigWig>(INCOMPLETE, "bigWig file %s: its index lists %llu data blocks, but its header counts %llu",
		                 fname, (unsigned long long)num_blocks, (unsigned long long)get<uint64_t>(data_count));

	// A chromosome the tree lists without data (libBigWig writers list every chromosome they
	// were given) has nothing to import, so it is neither matched nor reported as skipped.
	m_db2bw.assign(chromkey.get_num_chroms(), -1);
	for (uint32_t bwid = 0; bwid < num_chroms; ++bwid) {
		if (m_blocks[bwid].empty())
			continue;

		int chromid = chromkey.find_chromid(m_bw_names[bwid]);

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
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid chromosome tree", m_bfile.file_name().c_str());

	read_at(offset, node, sizeof(node));

	uint16_t count = get<uint16_t>(node + 2);
	uint64_t item_size = keysize + 8;   // leaf: key, chromId, chromSize; non-leaf: key, child offset

	if (m_tree_budget < sizeof(node) + count * item_size)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid chromosome tree", m_bfile.file_name().c_str());
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
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid data index", m_bfile.file_name().c_str());

	read_at(offset, node, sizeof(node));

	uint16_t count = get<uint16_t>(node + 2);
	uint64_t item_size = node[0] ? 32 : 24;   // 4 x bounds, then leaf: data offset, data size; non-leaf: child offset

	if (m_tree_budget < sizeof(node) + count * item_size)
		TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid data index", m_bfile.file_name().c_str());
	m_tree_budget -= sizeof(node) + count * item_size;

	vector<uint64_t> children;
	{
		vector<char> items(count * item_size);

		read_at(offset + sizeof(node), items.data(), items.size());
		for (uint16_t i = 0; i < count; ++i) {
			const char *item = &items[i * item_size];

			if (node[0]) {
				uint32_t bwid = get<uint32_t>(item);   // a bigWig data block holds a single chromosome
				Block block{get<uint64_t>(item + 16), get<uint64_t>(item + 24)};
				uint64_t file_size = m_bfile.file_size();

				if (bwid >= m_blocks.size() || block.size > file_size || block.offset > file_size - block.size)
					TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: invalid data index", m_bfile.file_name().c_str());
				m_blocks[bwid].push_back(block);
			} else
				children.push_back(get<uint64_t>(item + 16));
		}
	}   // items is freed before recursing, as in read_chrom_tree

	for (vector<uint64_t>::const_iterator ichild = children.begin(); ichild != children.end(); ++ichild)
		read_rtree(*ichild, depth + 1);
}

bool BigWig::get_data(int chromid, uint64_t iblock, GIntervals &intervals)
{
	intervals.clear();
	if (m_db2bw[chromid] < 0 || iblock >= m_blocks[m_db2bw[chromid]].size())
		return false;

	const char *fname = m_bfile.file_name().c_str();
	uint32_t bwid = m_db2bw[chromid];
	int64_t chromsize = m_chromkey->get_chrom_size(chromid);
	const Block &block = m_blocks[bwid][iblock];

	if (!iblock)
		m_last_end = 0;

	const char *data;
	uint64_t size;

	m_cbuf.resize(block.size);
	read_at(block.offset, m_cbuf.data(), block.size);

	if (m_uncompress_bufsize) {
		uLongf usize = min((uint64_t)m_uncompress_bufsize, MAX_BLOCK_SIZE);
		int rc;

		m_ubuf.resize(usize);
		rc = uncompress((Bytef *)m_ubuf.data(), &usize, (const Bytef *)m_cbuf.data(), block.size);
		if (rc != Z_OK)
			TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: data block at offset %llu failed to decompress (zlib error %d)",
			                 fname, (unsigned long long)block.offset, rc);
		data = m_ubuf.data();
		size = usize;
	} else {
		data = m_cbuf.data();
		size = block.size;
	}

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

		// As in the other importers: a record that runs past the end of the chromosome is clipped to it,
		// and dropped if it starts at or past the end. The caller reports them.
		if (e > chromsize) {
			// the chromosomes are read in database order: the first is the one of the lowest offset
			if (m_first_clipped.empty() || block.offset < m_first_clipped_offset ||
			    (block.offset == m_first_clipped_offset && i < m_first_clipped_item)) {
				m_first_clipped_offset = block.offset;
				m_first_clipped_item = i;
				char buf[1000];
				snprintf(buf, sizeof(buf), "(%s, %lld, %lld), chromosome size %lld", m_bw_names[bwid].c_str(), (long long)s,
				         (long long)e, (long long)chromsize);
				m_first_clipped = buf;
			}
			if (s >= chromsize) {
				++m_num_dropped;
				continue;
			}
			++m_num_clipped;
		}

		if (s < m_last_end)
			TGLError<BigWig>(INVALID_FORMAT, "bigWig file %s: records of %s are overlapping or unsorted at %lld",
			                 fname, m_bw_names[bwid].c_str(), (long long)s);

		m_last_end = min(e, chromsize);
		intervals.push_back(GInterval(chromid, s, m_last_end, 0, GInterval::cast2udata(v)));
	}
	return true;
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
