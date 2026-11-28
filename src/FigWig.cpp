/*
 * FigWig.cpp
 *
 * A reader of figwig files, version 1. Every check below is one the specification
 * (docs/format.rst of https://github.com/jmschrei/figwig) asks a reader to make before it
 * reads or allocates by a value of the file; a file that fails one is refused as corrupt.
 */

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <set>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <zlib.h>

#include "FigWig.h"
#include "TGLException.h"

using namespace std;

namespace {

const uint32_t MAGIC = 0x57474946;  // "FIGW"
const uint64_t HEADER_SIZE = 128;
const uint64_t MAX_TILES = (uint64_t)1 << 26;
const int MAX_JSON_DEPTH = 512;

// element codes: u1 u2 u4 i1 i2 i4 f4 one f8
const int ELEM_WIDTH[9] = { 1, 2, 4, 1, 2, 4, 4, 0, 8 };
const int ELEM_ONE = 7;

inline uint16_t get_u2(const unsigned char *p) { uint16_t v; memcpy(&v, p, 2); return v; }
inline uint32_t get_u4(const unsigned char *p) { uint32_t v; memcpy(&v, p, 4); return v; }
inline uint64_t get_u8(const unsigned char *p) { uint64_t v; memcpy(&v, p, 8); return v; }

inline uint64_t pad4(uint64_t x) { return (x + 3) & ~(uint64_t)3; }
inline uint64_t pad8(uint64_t x) { return (x + 7) & ~(uint64_t)7; }

inline uint32_t crc32_of(const unsigned char *p, uint64_t n)
{
	uLong crc = crc32(0L, Z_NULL, 0);
	while (n) {
		uInt step = (uInt)min<uint64_t>(n, 1u << 30);
		crc = crc32(crc, p, step);
		p += step;
		n -= step;
	}
	return (uint32_t)crc;
}

// The value of element i of an array of element type 'code' (not the one code) at p
inline double elem(const unsigned char *p, int code, uint64_t i)
{
	switch (code) {
	case 0: return p[i];
	case 1: return get_u2(p + 2 * i);
	case 2: return get_u4(p + 4 * i);
	case 3: return (int8_t)p[i];
	case 4: { int16_t v; memcpy(&v, p + 2 * i, 2); return v; }
	case 5: { int32_t v; memcpy(&v, p + 4 * i, 4); return v; }
	case 6: { float v; memcpy(&v, p + 4 * i, 4); return v; }
	case 8: { double v; memcpy(&v, p + 8 * i, 8); return v; }
	default: return 1;
	}
}

// Inflates the zlib stream src[0, n) into exactly 'expected' bytes at dst; false if it is not a
// valid stream of that length with nothing after its end
bool inflate_exact(const unsigned char *src, uint64_t n, unsigned char *dst, uint64_t expected)
{
	if (n > numeric_limits<uInt>::max() || expected > numeric_limits<uInt>::max())
		return false;

	unsigned char dummy;
	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	if (inflateInit(&zs) != Z_OK)
		return false;
	zs.next_in = (Bytef *)src;
	zs.avail_in = (uInt)n;
	zs.next_out = expected ? dst : &dummy;
	zs.avail_out = expected ? (uInt)expected : 0;
	int ret = inflate(&zs, Z_FINISH);
	if (ret == Z_BUF_ERROR && !expected && zs.avail_in) {
		// an empty stream may need room it does not fill
		zs.next_out = &dummy;
		zs.avail_out = 1;
		ret = inflate(&zs, Z_FINISH);
	}
	bool ok = ret == Z_STREAM_END && zs.avail_in == 0 && zs.total_out == expected;
	inflateEnd(&zs);
	return ok;
}

//------------------------------------------------- JSON -------------------------------------------------

struct JVal {
	enum Type { NUL, BOOL, NUM, STR, ARR, OBJ };

	Type                              type{NUL};
	string                            s;   // STR: decoded UTF-8; NUM: the token; BOOL: true / false
	vector<JVal>                      arr;
	vector<pair<string, JVal>>        obj;

	const JVal *get(const char *key) const {
		for (const auto &kv : obj) {
			if (kv.first == key)
				return &kv.second;
		}
		return NULL;
	}
};

class JParser {
public:
	JParser(const string &text) : m_p(text.data()), m_end(text.data() + text.size()) {}

	bool parse(JVal &v) {
		ws();
		if (!value(v, 0))
			return false;
		ws();
		return m_p == m_end;
	}

private:
	const char *m_p;
	const char *m_end;

	void ws() { while (m_p < m_end && (*m_p == ' ' || *m_p == '\t' || *m_p == '\n' || *m_p == '\r')) ++m_p; }

	bool lit(const char *word) {
		size_t n = strlen(word);
		if ((size_t)(m_end - m_p) < n || strncmp(m_p, word, n))
			return false;
		m_p += n;
		return true;
	}

	static void put_utf8(string &out, uint32_t cp) {
		if (cp < 0x80)
			out += (char)cp;
		else if (cp < 0x800) {
			out += (char)(0xc0 | (cp >> 6));
			out += (char)(0x80 | (cp & 0x3f));
		} else if (cp < 0x10000) {
			out += (char)(0xe0 | (cp >> 12));
			out += (char)(0x80 | ((cp >> 6) & 0x3f));
			out += (char)(0x80 | (cp & 0x3f));
		} else {
			out += (char)(0xf0 | (cp >> 18));
			out += (char)(0x80 | ((cp >> 12) & 0x3f));
			out += (char)(0x80 | ((cp >> 6) & 0x3f));
			out += (char)(0x80 | (cp & 0x3f));
		}
	}

	bool hex4(uint32_t &cp) {
		if (m_end - m_p < 4)
			return false;
		cp = 0;
		for (int i = 0; i < 4; ++i) {
			char c = *m_p++;
			cp <<= 4;
			if (c >= '0' && c <= '9') cp |= c - '0';
			else if (c >= 'a' && c <= 'f') cp |= c - 'a' + 10;
			else if (c >= 'A' && c <= 'F') cp |= c - 'A' + 10;
			else return false;
		}
		return true;
	}

	bool str(string &out) {
		if (m_p >= m_end || *m_p != '"')
			return false;
		++m_p;
		while (m_p < m_end && *m_p != '"') {
			unsigned char c = *m_p;
			if (c < 0x20)
				return false;
			if (c != '\\') {
				out += (char)c;
				++m_p;
				continue;
			}
			if (++m_p >= m_end)
				return false;
			char e = *m_p++;
			switch (e) {
			case '"': out += '"'; break;
			case '\\': out += '\\'; break;
			case '/': out += '/'; break;
			case 'b': out += '\b'; break;
			case 'f': out += '\f'; break;
			case 'n': out += '\n'; break;
			case 'r': out += '\r'; break;
			case 't': out += '\t'; break;
			case 'u': {
				uint32_t cp;
				if (!hex4(cp))
					return false;
				if (cp >= 0xd800 && cp < 0xdc00 && m_end - m_p >= 6 && m_p[0] == '\\' && m_p[1] == 'u') {
					const char *save = m_p;
					uint32_t lo;
					m_p += 2;
					if (hex4(lo) && lo >= 0xdc00 && lo < 0xe000)
						cp = 0x10000 + ((cp - 0xd800) << 10) + (lo - 0xdc00);
					else
						m_p = save;
				}
				put_utf8(out, cp);
				break;
			}
			default:
				return false;
			}
		}
		if (m_p >= m_end)
			return false;
		++m_p;
		return true;
	}

	bool number(string &out) {
		const char *start = m_p;
		auto digits = [&]() { const char *d = m_p; while (m_p < m_end && *m_p >= '0' && *m_p <= '9') ++m_p; return m_p > d; };

		if (m_p < m_end && *m_p == '-')
			++m_p;
		if (m_p < m_end && *m_p == '0')
			++m_p;
		else if (!digits())
			return false;
		if (m_p < m_end && *m_p == '.') {
			++m_p;
			if (!digits())
				return false;
		}
		if (m_p < m_end && (*m_p == 'e' || *m_p == 'E')) {
			++m_p;
			if (m_p < m_end && (*m_p == '+' || *m_p == '-'))
				++m_p;
			if (!digits())
				return false;
		}
		out.assign(start, m_p - start);
		return true;
	}

	bool value(JVal &v, int depth) {
		if (depth > MAX_JSON_DEPTH || m_p >= m_end)
			return false;
		char c = *m_p;
		if (c == '{') {
			v.type = JVal::OBJ;
			++m_p;
			ws();
			if (m_p < m_end && *m_p == '}') {
				++m_p;
				return true;
			}
			while (1) {
				string key;
				JVal item;
				ws();
				if (!str(key))
					return false;
				ws();
				if (m_p >= m_end || *m_p != ':')
					return false;
				++m_p;
				ws();
				if (!value(item, depth + 1))
					return false;
				// a repeated key keeps its first place and its last value, as Python's json does
				bool found = false;
				for (auto &kv : v.obj) {
					if (kv.first == key) {
						kv.second = std::move(item);
						found = true;
						break;
					}
				}
				if (!found)
					v.obj.emplace_back(std::move(key), std::move(item));
				ws();
				if (m_p < m_end && *m_p == ',') {
					++m_p;
					continue;
				}
				if (m_p < m_end && *m_p == '}') {
					++m_p;
					return true;
				}
				return false;
			}
		}
		if (c == '[') {
			v.type = JVal::ARR;
			++m_p;
			ws();
			if (m_p < m_end && *m_p == ']') {
				++m_p;
				return true;
			}
			while (1) {
				v.arr.emplace_back();
				ws();
				if (!value(v.arr.back(), depth + 1))
					return false;
				ws();
				if (m_p < m_end && *m_p == ',') {
					++m_p;
					continue;
				}
				if (m_p < m_end && *m_p == ']') {
					++m_p;
					return true;
				}
				return false;
			}
		}
		if (c == '"') {
			v.type = JVal::STR;
			return str(v.s);
		}
		if (lit("true")) {
			v.type = JVal::BOOL;
			v.s = "true";
			return true;
		}
		if (lit("false")) {
			v.type = JVal::BOOL;
			v.s = "false";
			return true;
		}
		if (lit("null")) {
			v.type = JVal::NUL;
			return true;
		}
		v.type = JVal::NUM;
		return number(v.s);
	}
};

// A string as Python's json.dumps writes it: ASCII only, with \uXXXX escapes
void dump_str(const string &s, string &out)
{
	static const char *hex = "0123456789abcdef";
	auto u = [&](uint32_t cp) {
		out += "\\u";
		for (int shift = 12; shift >= 0; shift -= 4)
			out += hex[(cp >> shift) & 0xf];
	};

	out += '"';
	for (size_t i = 0; i < s.size();) {
		unsigned char c = s[i];
		uint32_t cp = c;
		size_t len = 1;
		if (c >= 0xf0 && i + 3 < s.size() + 0) { cp = ((c & 0x07) << 18) | ((s[i + 1] & 0x3f) << 12) | ((s[i + 2] & 0x3f) << 6) | (s[i + 3] & 0x3f); len = 4; }
		else if (c >= 0xe0 && i + 2 < s.size()) { cp = ((c & 0x0f) << 12) | ((s[i + 1] & 0x3f) << 6) | (s[i + 2] & 0x3f); len = 3; }
		else if (c >= 0xc0 && i + 1 < s.size()) { cp = ((c & 0x1f) << 6) | (s[i + 1] & 0x3f); len = 2; }
		i += len;

		if (cp == '"') out += "\\\"";
		else if (cp == '\\') out += "\\\\";
		else if (cp == '\n') out += "\\n";
		else if (cp == '\r') out += "\\r";
		else if (cp == '\t') out += "\\t";
		else if (cp == '\b') out += "\\b";
		else if (cp == '\f') out += "\\f";
		else if (cp < 0x20 || cp > 0x7e) {
			if (cp >= 0x10000) {
				cp -= 0x10000;
				u(0xd800 + (cp >> 10));
				u(0xdc00 + (cp & 0x3ff));
			} else
				u(cp);
		} else
			out += (char)cp;
	}
	out += '"';
}

// A value as Python's json.dumps writes it (numbers as their token in the file)
void dump(const JVal &v, string &out)
{
	switch (v.type) {
	case JVal::NUL: out += "null"; break;
	case JVal::BOOL:
	case JVal::NUM: out += v.s; break;
	case JVal::STR: dump_str(v.s, out); break;
	case JVal::ARR:
		out += '[';
		for (size_t i = 0; i < v.arr.size(); ++i) {
			if (i)
				out += ", ";
			dump(v.arr[i], out);
		}
		out += ']';
		break;
	case JVal::OBJ:
		out += '{';
		for (size_t i = 0; i < v.obj.size(); ++i) {
			if (i)
				out += ", ";
			dump_str(v.obj[i].first, out);
			out += ": ";
			dump(v.obj[i].second, out);
		}
		out += '}';
		break;
	}
}

// varints of at most three bytes: offsets within a tile are below 2^16
inline bool read_varint(const unsigned char *&p, const unsigned char *end, uint32_t &v)
{
	v = 0;
	for (int shift = 0; shift < 21; shift += 7) {
		if (p >= end)
			return false;
		unsigned char b = *p++;
		v |= (uint32_t)(b & 0x7f) << shift;
		if (!(b & 0x80))
			return true;
	}
	return false;
}

} // namespace

//------------------------------------------------- FigWig -------------------------------------------------

FigWig::~FigWig()
{
	if (m_fd >= 0)
		::close(m_fd);
}

void FigWig::corrupt(const char *fmt, ...) const
{
	char buf[1000];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);
	TGLError<FigWig>("figwig file %s is corrupt or unsupported: %s", m_path.c_str(), buf);
}

void FigWig::read_at(uint64_t offset, void *buf, uint64_t size, const char *what)
{
	if (offset > m_length || size > m_length - offset)
		corrupt("%s lies past the end of the file", what);

	char *p = (char *)buf;
	while (size) {
		ssize_t n = pread(m_fd, p, min<uint64_t>(size, (uint64_t)1 << 30), (off_t)offset);
		if (n < 0) {
			if (errno == EINTR)
				continue;
			TGLError<FigWig>("Reading figwig file %s: %s", m_path.c_str(), strerror(errno));
		}
		if (n == 0)
			corrupt("%s lies past the end of the file", what);
		p += n;
		offset += n;
		size -= n;
	}
}

void FigWig::open(const char *path)
{
	{
		uint16_t probe = 1;
		unsigned char first;
		memcpy(&first, &probe, 1);
		if (first != 1)
			TGLError<FigWig>("Reading figwig files needs a little-endian machine");
	}

	m_path = path;
	m_fd = ::open(path, O_RDONLY);
	if (m_fd < 0)
		TGLError<FigWig>("Opening figwig file %s: %s", path, strerror(errno));

	struct stat st;
	if (fstat(m_fd, &st))
		TGLError<FigWig>("Reading figwig file %s: %s", path, strerror(errno));
	uint64_t file_size = st.st_size;
	if (file_size < HEADER_SIZE)
		TGLError<FigWig>("%s is not a figwig file", path);

	m_length = file_size;
	unsigned char h[HEADER_SIZE];
	read_at(0, h, HEADER_SIZE, "the header");

	if (get_u4(h) != MAGIC)
		TGLError<FigWig>("%s is not a figwig file", path);

	{
		unsigned char copy[HEADER_SIZE];
		memcpy(copy, h, HEADER_SIZE);
		memset(copy + 104, 0, 4);
		if (crc32_of(copy, HEADER_SIZE) != get_u4(h + 104))
			corrupt("the header's CRC-32 does not match");
	}

	uint16_t version = get_u2(h + 4);
	if (version != 1)
		TGLError<FigWig>("figwig file %s is of format version %u, and only version 1 can be read: it was written by a newer figwig",
		                 path, version);
	if (get_u4(h + 96))
		corrupt("the header names required features (0x%x) this reader does not know", get_u4(h + 96));

	uint16_t flags = get_u2(h + 6);
	if (flags & ~1u)
		corrupt("unknown header flags 0x%x", flags);
	m_track_order = flags & 1;

	m_tile = get_u4(h + 8);
	if (m_tile < 16 || m_tile > 65536 || (m_tile & (m_tile - 1)))
		corrupt("invalid tile size %u", m_tile);

	m_codec = get_u2(h + 12);
	if (m_codec > 1)
		corrupt("unknown codec %u", m_codec);

	uint32_t n_chroms = get_u4(h + 16);
	uint32_t n_tracks = get_u4(h + 20);
	m_num_groups = get_u4(h + 24);
	uint32_t tables_crc = get_u4(h + 28);
	uint64_t chroms_off = get_u8(h + 32);
	uint64_t tracks_off = get_u8(h + 40);
	uint64_t index_off = get_u8(h + 48);
	uint64_t metadata_off = get_u8(h + 64);
	uint64_t metadata_len = get_u8(h + 72);
	uint64_t length = get_u8(h + 88);

	if (length < HEADER_SIZE || length > file_size)
		corrupt("invalid length %llu (the file has %llu bytes)", (unsigned long long)length, (unsigned long long)file_size);
	m_length = length;
	if (!(HEADER_SIZE <= chroms_off && chroms_off <= tracks_off && tracks_off <= metadata_off && metadata_off <= index_off && index_off <= length))
		corrupt("the table offsets are out of order");
	if (!n_chroms || !n_tracks || !m_num_groups)
		corrupt("no chromosomes, tracks or groups");
	if ((uint64_t)m_num_groups * 8 > length - index_off)
		corrupt("the index directory lies past the end of the file");
	if (metadata_len > index_off - metadata_off)
		corrupt("the metadata overlaps the index");
	if (m_track_order && m_num_groups != n_tracks)
		corrupt("a track-order file must have one group per track");

	m_chroms_offset = chroms_off;
	m_index_offset = index_off;

	// the four tables, checked together
	uint64_t tables_end = index_off + 8 * (uint64_t)m_num_groups;
	vector<unsigned char> tables(tables_end - chroms_off);
	read_at(chroms_off, tables.data(), tables.size(), "the tables");
	if (crc32_of(tables.data(), tables.size()) != tables_crc)
		corrupt("the tables' CRC-32 does not match");

	// the chromosome table
	{
		const unsigned char *p = tables.data();
		const unsigned char *end = tables.data() + (tracks_off - chroms_off);
		set<string> names;
		uint64_t first_tile = 0;

		for (uint32_t i = 0; i < n_chroms; ++i) {
			if (end - p < 2)
				corrupt("the chromosome table is truncated");
			uint16_t len = get_u2(p);
			p += 2;
			if (!len || (uint64_t)(end - p) < (uint64_t)len + 8)
				corrupt("the chromosome table is truncated or holds an empty name");
			Chrom chrom;
			chrom.name.assign((const char *)p, len);
			if (chrom.name.find('\0') != string::npos)
				corrupt("a chromosome name holds a NUL byte");
			p += len;
			chrom.size = get_u8(p);
			p += 8;
			if (!chrom.size || chrom.size > 0xffffffffULL)
				corrupt("chromosome %s has an invalid length", chrom.name.c_str());
			if (!names.insert(chrom.name).second)
				corrupt("chromosome %s appears twice", chrom.name.c_str());
			chrom.first_tile = first_tile;
			first_tile += (chrom.size + m_tile - 1) / m_tile;
			m_chroms.push_back(chrom);
		}
		if (p != end)
			corrupt("the chromosome table does not end where the statistics table starts");
		if (first_tile > MAX_TILES)
			corrupt("the genome has more than 2^26 tiles");
		m_num_tiles = first_tile;
	}

	if ((uint64_t)n_tracks * 32 > metadata_off - tracks_off)
		corrupt("the statistics table is truncated");

	parse_metadata(string((const char *)tables.data() + (metadata_off - chroms_off), metadata_len));
	if (m_tracks.size() != n_tracks)
		corrupt("the metadata describes %llu tracks where the header has %u", (unsigned long long)m_tracks.size(), n_tracks);

	m_directory.resize(m_num_groups);
	for (uint32_t g = 0; g < m_num_groups; ++g) {
		m_directory[g] = get_u8(tables.data() + (index_off - chroms_off) + 8 * g);
		if (m_directory[g] < tables_end || m_directory[g] >= length)
			corrupt("the record of group %u lies outside the index", g);
	}

	if (m_track_order)
		m_blocks.resize(n_tracks);
	else
		m_live.resize(m_num_groups);
}

void FigWig::parse_metadata(const string &json)
{
	JVal doc;
	if (!JParser(json).parse(doc) || doc.type != JVal::OBJ)
		corrupt("the metadata is not a JSON object");

	const JVal *attributes = doc.get("attributes");
	if (attributes && attributes->type != JVal::OBJ)
		corrupt("the metadata's attributes are not an object");

	const JVal *genome = doc.get("genome");
	if (genome && genome->type == JVal::STR) {
		if (genome->s.find('\0') != string::npos)
			corrupt("the genome name holds a NUL byte");
		m_has_genome = true;
		m_genome = genome->s;
	}

	const JVal *tracks = doc.get("tracks");
	if (!tracks || tracks->type != JVal::ARR)
		corrupt("the metadata has no list of tracks");

	static const char *dtypes[] = { "uint8", "uint16", "uint32", "int8", "int16", "int32", "float32", "float64" };
	set<string> names;

	m_group_tracks.assign(m_num_groups, vector<int>());
	for (const JVal &rec : tracks->arr) {
		if (rec.type != JVal::OBJ)
			corrupt("a track record is not an object");
		const JVal *name = rec.get("name");
		const JVal *dtype = rec.get("dtype");
		const JVal *group = rec.get("group");
		if (!name || name->type != JVal::STR || name->s.empty() || name->s.find('\0') != string::npos)
			corrupt("a track has no name, or one holding a NUL byte");
		if (!names.insert(name->s).second)
			corrupt("track %s appears twice", name->s.c_str());
		if (!dtype || dtype->type != JVal::STR || find(begin(dtypes), end(dtypes), dtype->s) == end(dtypes))
			corrupt("track %s has an invalid dtype", name->s.c_str());
		if (!group || group->type != JVal::NUM || group->s.find_first_of(".eE-") != string::npos || group->s.size() > 10 ||
		    strtoull(group->s.c_str(), NULL, 10) >= m_num_groups)
			corrupt("track %s has an invalid group", name->s.c_str());

		Track track;
		track.name = name->s;
		track.dtype = dtype->s;
		track.group = (int)strtoul(group->s.c_str(), NULL, 10);
		track.slot = m_group_tracks[track.group].size();
		m_group_tracks[track.group].push_back(m_tracks.size());

		for (const auto &kv : rec.obj) {
			if (kv.first == "name" || kv.first == "dtype" || kv.first == "group")
				continue;
			Attr attr;
			attr.key = kv.first;
			attr.null = kv.second.type == JVal::NUL;
			if (kv.second.type == JVal::STR)
				attr.text = kv.second.s;
			else
				dump(kv.second, attr.text);
			track.attrs.push_back(attr);
		}
		m_tracks.push_back(track);
	}

	for (uint32_t g = 0; g < m_num_groups; ++g) {
		if (m_group_tracks[g].empty())
			corrupt("group %u holds no track", g);
		if (m_track_order && m_group_tracks[g].size() != 1)
			corrupt("a group of a track-order file holds several tracks");
	}
}

void FigWig::decode(uint64_t offset, uint32_t size, uint32_t raw, uint32_t crc, bool stored, const char *what)
{
	if (offset < HEADER_SIZE || offset > m_chroms_offset || size > m_chroms_offset - offset)
		corrupt("a %s lies outside the data region", what);

	m_decoded.resize(raw);
	if (!m_codec || stored) {
		if (size != raw)
			corrupt("a stored %s has %u bytes where it should have %u", what, size, raw);
		read_at(offset, m_decoded.data(), size, what);
		if (m_codec && crc32_of(m_decoded.data(), size) != crc)
			corrupt("a stored %s does not match its CRC-32", what);
		return;
	}
	m_filebuf.resize(size);
	read_at(offset, m_filebuf.data(), size, what);
	if (!inflate_exact(m_filebuf.data(), size, m_decoded.data(), raw))
		corrupt("a %s is not a valid zlib stream of %u bytes", what, raw);
}

void FigWig::load_live(int group)
{
	LiveIndex &li = m_live[group];
	if (li.loaded)
		return;

	uint64_t rec = m_directory[group];
	unsigned char head[8];
	read_at(rec, head, 8, "a group's record");
	uint32_t n_runs = get_u4(head);
	uint32_t crc = get_u4(head + 4);

	if ((uint64_t)n_runs * 32 + 8 > m_length - rec - 8)
		corrupt("the runs of group %d lie past the end of the file", group);
	vector<unsigned char> buf((uint64_t)n_runs * 32 + 8);
	read_at(rec + 8, buf.data(), buf.size(), "a group's record");
	uint32_t n_live = get_u4(buf.data() + 32 * (uint64_t)n_runs);
	uint32_t n_bytes = get_u4(buf.data() + 32 * (uint64_t)n_runs + 4);
	if (n_live > m_num_tiles)
		corrupt("group %d lists more live tiles than the genome has", group);
	buf.resize(buf.size() + n_bytes);
	read_at(rec + 8 + (uint64_t)n_runs * 32 + 8, buf.data() + (uint64_t)n_runs * 32 + 8, n_bytes, "a group's record");
	if (crc32_of(buf.data(), buf.size()) != crc)
		corrupt("the record of group %d does not match its CRC-32", group);

	vector<unsigned char> live(17 * (uint64_t)n_live);
	if (!inflate_exact(buf.data() + (uint64_t)n_runs * 32 + 8, n_bytes, live.data(), live.size()))
		corrupt("the live tiles of group %d are not a valid zlib stream", group);

	uint64_t tracks = m_group_tracks[group].size();
	uint64_t max_raw = 8 * tracks * m_tile + 16 * tracks + 64;
	const unsigned char *p = live.data();

	li.tiles.resize(n_live);
	li.kinds.resize(n_live);
	li.sizes.resize(n_live);
	li.raws.resize(n_live);
	li.crcs.resize(n_live);
	li.offsets.resize(n_live);
	for (uint32_t i = 0; i < n_live; ++i) {
		li.tiles[i] = get_u4(p + 4 * (uint64_t)i);
		li.kinds[i] = p[4 * (uint64_t)n_live + i];
		li.sizes[i] = get_u4(p + 5 * (uint64_t)n_live + 4 * (uint64_t)i);
		li.raws[i] = get_u4(p + 9 * (uint64_t)n_live + 4 * (uint64_t)i);
		li.crcs[i] = get_u4(p + 13 * (uint64_t)n_live + 4 * (uint64_t)i);

		if (li.tiles[i] >= m_num_tiles || (i && li.tiles[i] <= li.tiles[i - 1]))
			corrupt("the live tiles of group %d are out of order or out of range", group);
		uint8_t kind = li.kinds[i] & 7;
		uint8_t code = (li.kinds[i] >> 3) & 15;
		bool stored = li.kinds[i] & 0x80;
		if (kind > 4 || code > 8 || ((kind == 1 || kind == 4) && code == ELEM_ONE))
			corrupt("a chunk of group %d has an invalid kind byte 0x%x", group, li.kinds[i]);
		if (stored && !m_codec)
			corrupt("a chunk of group %d is marked stored in a file of codec 0", group);
		if (!li.sizes[i] || !li.raws[i] || li.raws[i] % 8 || li.raws[i] > max_raw || li.sizes[i] > max_raw + 64)
			corrupt("a chunk of group %d has an invalid size", group);
		if ((!m_codec || stored) && li.sizes[i] != li.raws[i])
			corrupt("a stored chunk of group %d has a size other than its decoded size", group);
	}

	// the runs: each live tile in exactly one, its chunks back to back
	struct RunRec { uint64_t first, last, offset, bytes; };
	vector<RunRec> runs(n_runs);
	for (uint32_t r = 0; r < n_runs; ++r) {
		const unsigned char *q = buf.data() + 32 * (uint64_t)r;
		runs[r] = RunRec{ get_u8(q), get_u8(q + 8), get_u8(q + 16), get_u8(q + 24) };
		if (runs[r].first > runs[r].last)
			corrupt("a run of group %d is inverted", group);
	}
	sort(runs.begin(), runs.end(), [](const RunRec &a, const RunRec &b) { return a.first < b.first; });

	uint32_t k = 0;
	for (uint32_t r = 0; r < n_runs; ++r) {
		if (r && runs[r].first <= runs[r - 1].last)
			corrupt("the runs of group %d overlap", group);
		if (k < n_live && li.tiles[k] < runs[r].first)
			corrupt("a live tile of group %d lies in no run", group);
		uint64_t offset = runs[r].offset;
		uint64_t bytes = 0;
		for (; k < n_live && li.tiles[k] <= runs[r].last; ++k) {
			li.offsets[k] = offset + bytes;
			bytes += li.sizes[k];
		}
		if (bytes != runs[r].bytes)
			corrupt("the chunks of a run of group %d do not add up to its bytes", group);
		if (bytes && (offset < HEADER_SIZE || offset > m_chroms_offset || bytes > m_chroms_offset - offset))
			corrupt("a run of group %d lies outside the data region", group);
	}
	if (k != n_live)
		corrupt("a live tile of group %d lies in no run", group);

	li.loaded = true;
}

bool FigWig::read_tile(int group, uint64_t tile, const vector<int> &slots, double *out)
{
	load_live(group);
	LiveIndex &li = m_live[group];

	auto it = lower_bound(li.tiles.begin(), li.tiles.end(), (uint32_t)tile);
	if (it == li.tiles.end() || *it != tile)
		return false;
	size_t k = it - li.tiles.begin();
	uint8_t kind = li.kinds[k] & 7;
	int code = (li.kinds[k] >> 3) & 15;
	int w = ELEM_WIDTH[code];
	uint64_t tracks = m_group_tracks[group].size();
	uint64_t T = m_tile;

	if (kind == 0)
		return false;

	decode(li.offsets[k], li.sizes[k], li.raws[k], li.crcs[k], li.kinds[k] & 0x80, "chunk");
	const unsigned char *d = m_decoded.data();
	uint64_t raw = li.raws[k];

	for (size_t i = 0; i < slots.size(); ++i)
		fill(out + i * T, out + (i + 1) * T, 0.);

	if (kind == 1) {  // dense
		if (tracks * T * w > raw)
			corrupt("a dense chunk is shorter than its tracks");
		for (size_t i = 0; i < slots.size(); ++i) {
			double *o = out + i * T;
			for (uint64_t j = 0; j < T; ++j)
				o[j] = elem(d, code, (uint64_t)slots[i] * T + j);
		}
		return true;
	}

	if (kind == 2 || kind == 4) {  // sparse or runs
		if (4 + 4 * tracks > raw)
			corrupt("a chunk is shorter than its header");
		uint32_t total = get_u4(d);
		uint64_t sum_lengths = 0;
		for (uint64_t s = 0; s < tracks; ++s)
			sum_lengths += get_u4(d + 4 + 4 * s);
		uint64_t values_start = pad4(4 + 4 * tracks + sum_lengths);
		if (values_start > raw || (uint64_t)total * w > raw - values_start)
			corrupt("a chunk is shorter than its entries");
		const unsigned char *vals = d + values_start;

		// the entries of each slot, and where its offsets and values start
		vector<uint64_t> first_entry(tracks + 1, 0), first_byte(tracks + 1, 0);
		const unsigned char *offs = d + 4 + 4 * tracks;
		for (uint64_t s = 0; s < tracks; ++s) {
			uint32_t len = get_u4(d + 4 + 4 * s);
			uint64_t n = 0;
			for (uint32_t b = 0; b < len; ++b)
				n += !(offs[first_byte[s] + b] & 0x80);
			first_byte[s + 1] = first_byte[s] + len;
			first_entry[s + 1] = first_entry[s] + n;
		}
		if (first_entry[tracks] != total)
			corrupt("a chunk's entries do not add up to its total");

		for (size_t i = 0; i < slots.size(); ++i) {
			int s = slots[i];
			const unsigned char *p = offs + first_byte[s];
			const unsigned char *pend = offs + first_byte[s + 1];
			uint64_t e = first_entry[s];
			uint64_t pos = 0;
			double *o = out + i * T;
			double cur = 0;
			uint64_t from = 0;
			bool first = true;

			while (p < pend) {
				uint32_t v;
				if (!read_varint(p, pend, v))
					corrupt("a chunk holds an invalid varint");
				pos = first ? v : pos + v;
				if (pos >= T || (!first && !v))
					corrupt("a chunk's offsets are out of range or out of order");
				double value = code == ELEM_ONE ? 1. : elem(vals, code, e);
				if (kind == 2)
					o[pos] = value;
				else {
					for (uint64_t j = from; j < pos; ++j)
						o[j] = cur;
					cur = value;
					from = pos;
				}
				first = false;
				++e;
			}
			if (kind == 4) {
				for (uint64_t j = from; j < T; ++j)
					o[j] = cur;
			}
		}
		return true;
	}

	// matrix
	if (raw < 4)
		corrupt("a chunk is shorter than its header");
	uint32_t n_present = get_u4(d);
	if (n_present > tracks || 8 + 8 * (uint64_t)n_present > raw)
		corrupt("a matrix chunk is shorter than its header");
	const unsigned char *present = d + 4;
	const unsigned char *ptr = d + 4 + 4 * (uint64_t)n_present;
	uint32_t total = get_u4(ptr + 4 * (uint64_t)n_present);
	uint64_t offs_start = 8 + 8 * (uint64_t)n_present;
	uint64_t values_start = pad4(offs_start + 2 * (uint64_t)total);
	if (values_start > raw || (uint64_t)total * w > raw - values_start || get_u4(ptr))
		corrupt("a matrix chunk is shorter than its entries");
	const unsigned char *offs = d + offs_start;
	const unsigned char *vals = d + values_start;
	vector<int> slot_index(tracks, -1);

	for (size_t i = 0; i < slots.size(); ++i)
		slot_index[slots[i]] = i;

	for (uint32_t j = 0; j < n_present; ++j) {
		uint32_t s = get_u4(present + 4 * (uint64_t)j);
		uint32_t a = get_u4(ptr + 4 * (uint64_t)j);
		uint32_t b = get_u4(ptr + 4 * (uint64_t)j + 4);
		if (s >= tracks || (j && s <= get_u4(present + 4 * (uint64_t)j - 4)) || b < a || b > total)
			corrupt("a matrix chunk's slots or pointers are invalid");
		if (slot_index[s] < 0)
			continue;
		double *o = out + (uint64_t)slot_index[s] * T;
		for (uint32_t e = a; e < b; ++e) {
			uint16_t pos = get_u2(offs + 2 * (uint64_t)e);
			if (pos >= T || (e > a && pos <= get_u2(offs + 2 * (uint64_t)e - 2)))
				corrupt("a matrix chunk's offsets are out of range or out of order");
			o[pos] = code == ELEM_ONE ? 1. : elem(vals, code, e);
		}
	}
	return true;
}

void FigWig::load_blocks(int track)
{
	BlockIndex &bi = m_blocks[track];
	if (bi.loaded)
		return;

	int group = m_tracks[track].group;
	uint64_t rec = m_directory[group];
	unsigned char head[8];
	read_at(rec, head, 8, "a track's record");
	uint32_t n_blocks = get_u4(head);
	uint32_t crc = get_u4(head + 4);
	if ((uint64_t)n_blocks * 40 > m_length - rec - 8)
		corrupt("the blocks of track %s lie past the end of the file", m_tracks[track].name.c_str());
	vector<unsigned char> buf((uint64_t)n_blocks * 40);
	read_at(rec + 8, buf.data(), buf.size(), "a track's record");
	if (crc32_of(buf.data(), buf.size()) != crc)
		corrupt("the record of track %s does not match its CRC-32", m_tracks[track].name.c_str());

	bi.blocks.resize(n_blocks);
	bi.chrom_start.assign(m_chroms.size() + 1, 0);
	for (uint32_t i = 0; i < n_blocks; ++i) {
		const unsigned char *q = buf.data() + 40 * (uint64_t)i;
		Block &b = bi.blocks[i];
		b = Block{ get_u4(q), get_u4(q + 4), get_u4(q + 8), get_u4(q + 12), get_u8(q + 16), get_u4(q + 24),
		           get_u4(q + 28), get_u4(q + 32), get_u4(q + 36) };
		int code = (b.kind >> 3) & 15;
		if (b.chrom >= m_chroms.size() || b.first > b.last || b.last >= m_chroms[b.chrom].size || !b.n || b.n > 65536)
			corrupt("a block of track %s is out of range", m_tracks[track].name.c_str());
		if (i && (b.chrom < bi.blocks[i - 1].chrom || (b.chrom == bi.blocks[i - 1].chrom && b.first <= bi.blocks[i - 1].last)))
			corrupt("the blocks of track %s are out of order", m_tracks[track].name.c_str());
		if (b.kind > 255 || (b.kind & 7) != 1 || code > 8 || (b.kind & 0x80))
			corrupt("a block of track %s has an invalid kind byte", m_tracks[track].name.c_str());
		if (b.raw != pad8(pad8(4 * (uint64_t)b.n) + (uint64_t)b.n * ELEM_WIDTH[code]) || !b.size || (!m_codec && b.size != b.raw))
			corrupt("a block of track %s has an invalid size", m_tracks[track].name.c_str());
		bi.chrom_start[b.chrom + 1] = i + 1;
	}
	for (size_t c = 1; c <= m_chroms.size(); ++c)
		bi.chrom_start[c] = max(bi.chrom_start[c], bi.chrom_start[c - 1]);
	bi.loaded = true;
}

uint64_t FigWig::num_blocks(int track, int chromid)
{
	load_blocks(track);
	const BlockIndex &bi = m_blocks[track];
	return bi.chrom_start[chromid + 1] - bi.chrom_start[chromid];
}

void FigWig::read_block(int track, int chromid, uint64_t iblock, vector<uint32_t> &positions, vector<double> &values)
{
	load_blocks(track);
	const Block &b = m_blocks[track].blocks[m_blocks[track].chrom_start[chromid] + iblock];
	int code = (b.kind >> 3) & 15;

	decode(b.offset, b.size, b.raw, b.crc, false, "block");
	const unsigned char *d = m_decoded.data();
	const unsigned char *vals = d + pad8(4 * (uint64_t)b.n);

	positions.resize(b.n);
	values.resize(b.n);
	for (uint32_t i = 0; i < b.n; ++i) {
		positions[i] = get_u4(d + 4 * (uint64_t)i);
		if (i && positions[i] <= positions[i - 1])
			corrupt("the positions of a block of track %s do not ascend", m_tracks[track].name.c_str());
		values[i] = code == ELEM_ONE ? 1. : elem(vals, code, i);
	}
	if (positions.front() != b.first || positions.back() != b.last)
		corrupt("a block of track %s does not match its record", m_tracks[track].name.c_str());
}
