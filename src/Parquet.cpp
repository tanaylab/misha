/*
 * Parquet.cpp
 *
 * Thrift compact protocol: https://github.com/apache/thrift/blob/master/doc/specs/thrift-compact-protocol.md
 * Encodings (PLAIN, RLE / bit-packed hybrid, RLE_DICTIONARY, BYTE_STREAM_SPLIT): apache/parquet-format, Encodings.md
 * Snappy: https://github.com/google/snappy/blob/main/format_description.txt
 */

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <map>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <unordered_map>
#include <zlib.h>

#include "Parquet.h"
#include "TGLException.h"
#include "rdbutils.h"

using namespace parquet;

// PLAIN is little-endian; values are copied to the file as they are in memory
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "Parquet output assumes a little-endian host");

//------------------------------ Thrift compact protocol ------------------------------

namespace {

enum { CT_TRUE = 1, CT_FALSE = 2, CT_I8 = 3, CT_I16 = 4, CT_I32 = 5, CT_I64 = 6, CT_DOUBLE = 7, CT_BINARY = 8, CT_LIST = 9, CT_SET = 10, CT_MAP = 11, CT_STRUCT = 12, CT_UUID = 13 };

void put_varint(string &out, uint64_t v)
{
	for (; v >= 0x80; v >>= 7)
		out += (char)(v | 0x80);
	out += (char)v;
}

uint64_t zigzag(int64_t v) { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }

// 'last' is the id of the previous field written in the same struct
void put_field(string &out, int16_t &last, int16_t id, int type)
{
	if (id > last && id - last <= 15)
		out += (char)((id - last) << 4 | type);
	else {
		out += (char)type;
		put_varint(out, zigzag(id));
	}
	last = id;
}

void put_i32(string &out, int16_t &last, int16_t id, int64_t v)
{
	put_field(out, last, id, CT_I32);
	put_varint(out, zigzag(v));
}

void put_i64(string &out, int16_t &last, int16_t id, int64_t v)
{
	put_field(out, last, id, CT_I64);
	put_varint(out, zigzag(v));
}

void put_binary(string &out, const string &v)
{
	put_varint(out, v.size());
	out += v;
}

void put_binary(string &out, int16_t &last, int16_t id, const string &v)
{
	put_field(out, last, id, CT_BINARY);
	put_binary(out, v);
}

void put_list(string &out, int16_t &last, int16_t id, int elem_type, size_t size)
{
	put_field(out, last, id, CT_LIST);
	if (size < 15)
		out += (char)(size << 4 | elem_type);
	else {
		out += (char)(0xf0 | elem_type);
		put_varint(out, size);
	}
}

template <typename T> void put_struct(string &out, int16_t &last, int16_t id, const T &v)
{
	put_field(out, last, id, CT_STRUCT);
	v.encode(out);
}

template <typename T> void put_struct_list(string &out, int16_t &last, int16_t id, const vector<T> &v)
{
	put_list(out, last, id, CT_STRUCT, v.size());
	for (const T &x : v)
		x.encode(out);
}

void put_stop(string &out) { out += '\0'; }

}

void Statistics::encode(string &out) const
{
	int16_t last = 0;
	put_i64(out, last, 3, null_count);
	put_binary(out, last, 5, max_value);
	put_binary(out, last, 6, min_value);
	put_stop(out);
}

void SchemaElement::encode(string &out) const
{
	int16_t last = 0;
	if (type >= 0)
		put_i32(out, last, 1, type);
	if (repetition_type >= 0)
		put_i32(out, last, 3, repetition_type);
	put_binary(out, last, 4, name);
	if (num_children >= 0)
		put_i32(out, last, 5, num_children);
	if (converted_type >= 0)
		put_i32(out, last, 6, converted_type);
	if (logical_type == LogicalType::STRING) {
		// LogicalType is a union: field 1 (STRING) holding an empty StringType struct
		put_field(out, last, 10, CT_STRUCT);
		int16_t union_last = 0;
		put_field(out, union_last, 1, CT_STRUCT);
		put_stop(out);
		put_stop(out);
	}
	put_stop(out);
}

void DataPageHeader::encode(string &out) const
{
	int16_t last = 0;
	put_i32(out, last, 1, num_values);
	put_i32(out, last, 2, encoding);
	put_i32(out, last, 3, definition_level_encoding);
	put_i32(out, last, 4, repetition_level_encoding);
	put_stop(out);
}

void DictionaryPageHeader::encode(string &out) const
{
	int16_t last = 0;
	put_i32(out, last, 1, num_values);
	put_i32(out, last, 2, encoding);
	put_stop(out);
}

void PageHeader::encode(string &out) const
{
	int16_t last = 0;
	put_i32(out, last, 1, type);
	put_i32(out, last, 2, uncompressed_page_size);
	put_i32(out, last, 3, compressed_page_size);
	if (type == PageType::DATA_PAGE)
		put_struct(out, last, 5, data_page_header);
	else if (type == PageType::DICTIONARY_PAGE)
		put_struct(out, last, 7, dictionary_page_header);
	put_stop(out);
}

void ColumnMetaData::encode(string &out) const
{
	int16_t last = 0;
	put_i32(out, last, 1, type);
	put_list(out, last, 2, CT_I32, encodings.size());
	for (int e : encodings)
		put_varint(out, zigzag(e));
	put_list(out, last, 3, CT_BINARY, path_in_schema.size());
	for (const string &s : path_in_schema)
		put_binary(out, s);
	put_i32(out, last, 4, codec);
	put_i64(out, last, 5, num_values);
	put_i64(out, last, 6, total_uncompressed_size);
	put_i64(out, last, 7, total_compressed_size);
	put_i64(out, last, 9, data_page_offset);
	if (dictionary_page_offset >= 0)
		put_i64(out, last, 11, dictionary_page_offset);
	if (has_statistics)
		put_struct(out, last, 12, statistics);
	put_stop(out);
}

void ColumnChunk::encode(string &out) const
{
	int16_t last = 0;
	put_i64(out, last, 2, file_offset);
	put_struct(out, last, 3, meta_data);
	put_stop(out);
}

void RowGroup::encode(string &out) const
{
	int16_t last = 0;
	put_struct_list(out, last, 1, columns);
	put_i64(out, last, 2, total_byte_size);
	put_i64(out, last, 3, num_rows);
	put_i64(out, last, 5, file_offset);
	put_i64(out, last, 6, total_compressed_size);
	put_stop(out);
}

void FileMetaData::encode(string &out) const
{
	int16_t last = 0;
	put_i32(out, last, 1, version);
	put_struct_list(out, last, 2, schema);
	put_i64(out, last, 3, num_rows);
	put_struct_list(out, last, 4, row_groups);
	if (!key_value_metadata.empty()) {
		// list<KeyValue>, KeyValue { 1: key; 2: value }
		put_list(out, last, 5, CT_STRUCT, key_value_metadata.size());
		for (const auto &kv : key_value_metadata) {
			int16_t kv_last = 0;
			put_binary(out, kv_last, 1, kv.first);
			put_binary(out, kv_last, 2, kv.second);
			put_stop(out);
		}
	}
	put_binary(out, last, 6, created_by);
	// ColumnOrder is a union: field 1 (TYPE_ORDER) holding an empty TypeDefinedOrder struct
	put_list(out, last, 7, CT_STRUCT, schema.empty() ? 0 : schema.size() - 1);
	for (size_t i = 1; i < schema.size(); ++i) {
		int16_t union_last = 0;
		put_field(out, union_last, 1, CT_STRUCT);
		put_stop(out);
		put_stop(out);
	}
	put_stop(out);
}

//------------------------------ Thrift compact protocol: decoding ------------------------------

// Decoding errors are thrown without the file, row group and column; ParquetReader adds them.
[[noreturn]] static void corrupt(const char *what)
{
	TGLError<ParquetReader>(ParquetReader::FORMAT_ERROR, "corrupt data (%s)", what);
}

static uint64_t get_varint(const uint8_t *&p, const uint8_t *end)
{
	uint64_t v = 0;
	for (int shift = 0; shift < 64; shift += 7) {
		if (p == end)
			corrupt("truncated varint");
		uint8_t b = *p++;
		v |= (uint64_t)(b & 0x7f) << shift;
		if (!(b & 0x80))
			return v;
	}
	corrupt("varint too long");
}

struct parquet::ThriftIn {
	const uint8_t *p;
	const uint8_t *end;
	int            depth = 0;

	ThriftIn(const void *data, size_t size) : p((const uint8_t *)data), end((const uint8_t *)data + size) {}

	void need(uint64_t n) { if ((uint64_t)(end - p) < n) corrupt("truncated"); }
	uint8_t byte() { need(1); return *p++; }

	uint64_t varint() { return get_varint(p, end); }

	int64_t zigzag_varint() { uint64_t v = varint(); return (int64_t)(v >> 1) ^ -(int64_t)(v & 1); }

	void expect(int type, int expected) { if (type != expected) corrupt("unexpected field type"); }

	int64_t i32(int type) {
		expect(type, CT_I32);
		int64_t v = zigzag_varint();
		if (v < INT32_MIN || v > INT32_MAX)
			corrupt("i32 out of range");
		return v;
	}
	int64_t i64(int type) { expect(type, CT_I64); return zigzag_varint(); }
	bool boolean(int type) {
		if (type != CT_TRUE && type != CT_FALSE)
			corrupt("unexpected field type");
		return type == CT_TRUE;
	}
	string binary(int type) {
		expect(type, CT_BINARY);
		uint64_t n = varint();
		need(n);
		string s((const char *)p, n);
		p += n;
		return s;
	}

	// Calls field(id, type) for every field of a struct; a field it does not take (returns
	// false) is skipped
	template <typename F> void read_struct(F field) {
		if (++depth > 64)
			corrupt("nested too deeply");
		int16_t last = 0;
		for (uint8_t h; (h = byte()); ) {
			int type = h & 0xf;
			int16_t id = (h >> 4) ? (int16_t)(last + (h >> 4)) : (int16_t)zigzag_varint();
			last = id;
			if (!field(id, type))
				skip(type);
		}
		--depth;
	}

	// Calls elem(type) for every element of a list or set
	template <typename F> void read_list(int type, F elem) {
		if (type != CT_LIST && type != CT_SET)
			corrupt("unexpected field type");
		uint8_t h = byte();
		uint64_t n = h >> 4;
		if (n == 15)
			n = varint();
		for (uint64_t i = 0; i < n; ++i)
			elem(h & 0xf);
	}

	template <typename T> void read_struct_list(int type, vector<T> &v) {
		v.clear();
		read_list(type, [&](int elem_type) {
			expect(elem_type, CT_STRUCT);
			v.emplace_back();
			v.back().decode(*this);
		});
	}

	// Skips a value of a field (bool: none) or of a container element (bool: a byte)
	void skip(int type, bool element = false) {
		switch (type) {
		case CT_TRUE:
		case CT_FALSE:
			if (element)
				byte();
			return;
		case CT_I8: byte(); return;
		case CT_I16:
		case CT_I32:
		case CT_I64: varint(); return;
		case CT_DOUBLE: need(8); p += 8; return;
		case CT_UUID: need(16); p += 16; return;
		case CT_BINARY: { uint64_t n = varint(); need(n); p += n; return; }
		case CT_STRUCT: read_struct([](int16_t, int) { return false; }); return;
		case CT_LIST:
		case CT_SET:
		case CT_MAP:
			if (++depth > 64)
				corrupt("nested too deeply");
			if (type == CT_MAP) {
				uint64_t n = varint();
				uint8_t kv = n ? byte() : 0;
				for (uint64_t i = 0; i < n; ++i) {
					skip(kv >> 4, true);
					skip(kv & 0xf, true);
				}
			} else
				read_list(type, [&](int elem_type) { skip(elem_type, true); });
			--depth;
			return;
		default:
			corrupt("unknown field type");
		}
	}
};

void Statistics::decode(ThriftIn &in)
{
	*this = Statistics();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 3: null_count = in.i64(ct); return true;
		case 5: max_value = in.binary(ct); return true;
		case 6: min_value = in.binary(ct); return true;
		}
		return false;
	});
}

void SchemaElement::decode(ThriftIn &in)
{
	*this = SchemaElement();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: type = in.i32(ct); return true;
		case 3: repetition_type = in.i32(ct); return true;
		case 4: name = in.binary(ct); return true;
		case 5: num_children = in.i32(ct); return true;
		case 6: converted_type = in.i32(ct); return true;
		case 10:  // LogicalType, a union: which member, and INTEGER's isSigned
			in.expect(ct, CT_STRUCT);
			in.read_struct([&](int16_t lid, int lct) {
				if (lct != CT_STRUCT)
					return false;
				logical_type = lid;
				if (lid != LogicalType::INTEGER)
					return false;
				in.read_struct([&](int16_t iid, int ict) {
					if (iid != 2)
						return false;
					integer_signed = in.boolean(ict);
					return true;
				});
				return true;
			});
			return true;
		}
		return false;
	});
}

void DataPageHeader::decode(ThriftIn &in)
{
	*this = DataPageHeader();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: num_values = in.i32(ct); return true;
		case 2: encoding = in.i32(ct); return true;
		case 3: definition_level_encoding = in.i32(ct); return true;
		case 4: repetition_level_encoding = in.i32(ct); return true;
		}
		return false;
	});
}

void DictionaryPageHeader::decode(ThriftIn &in)
{
	*this = DictionaryPageHeader();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: num_values = in.i32(ct); return true;
		case 2: encoding = in.i32(ct); return true;
		}
		return false;
	});
}

void DataPageHeaderV2::decode(ThriftIn &in)
{
	*this = DataPageHeaderV2();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: num_values = in.i32(ct); return true;
		case 2: num_nulls = in.i32(ct); return true;
		case 3: num_rows = in.i32(ct); return true;
		case 4: encoding = in.i32(ct); return true;
		case 5: definition_levels_byte_length = in.i32(ct); return true;
		case 6: repetition_levels_byte_length = in.i32(ct); return true;
		case 7: is_compressed = in.boolean(ct); return true;
		}
		return false;
	});
}

void PageHeader::decode(ThriftIn &in)
{
	*this = PageHeader();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: type = in.i32(ct); return true;
		case 2: uncompressed_page_size = in.i32(ct); return true;
		case 3: compressed_page_size = in.i32(ct); return true;
		case 5: in.expect(ct, CT_STRUCT); data_page_header.decode(in); return true;
		case 7: in.expect(ct, CT_STRUCT); dictionary_page_header.decode(in); return true;
		case 8: in.expect(ct, CT_STRUCT); data_page_header_v2.decode(in); return true;
		}
		return false;
	});
}

void ColumnMetaData::decode(ThriftIn &in)
{
	*this = ColumnMetaData();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: type = in.i32(ct); return true;
		case 2:
			in.read_list(ct, [&](int elem_type) { encodings.push_back(in.i32(elem_type)); });
			return true;
		case 3:
			in.read_list(ct, [&](int elem_type) { path_in_schema.push_back(in.binary(elem_type)); });
			return true;
		case 4: codec = in.i32(ct); return true;
		case 5: num_values = in.i64(ct); return true;
		case 6: total_uncompressed_size = in.i64(ct); return true;
		case 7: total_compressed_size = in.i64(ct); return true;
		case 9: data_page_offset = in.i64(ct); return true;
		case 11: dictionary_page_offset = in.i64(ct); return true;
		case 12: in.expect(ct, CT_STRUCT); statistics.decode(in); has_statistics = true; return true;
		}
		return false;
	});
}

void ColumnChunk::decode(ThriftIn &in)
{
	*this = ColumnChunk();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: file_path = in.binary(ct); return true;
		case 2: file_offset = in.i64(ct); return true;
		case 3: in.expect(ct, CT_STRUCT); meta_data.decode(in); return true;
		}
		return false;
	});
}

void RowGroup::decode(ThriftIn &in)
{
	*this = RowGroup();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: in.read_struct_list(ct, columns); return true;
		case 2: total_byte_size = in.i64(ct); return true;
		case 3: num_rows = in.i64(ct); return true;
		case 5: file_offset = in.i64(ct); return true;
		case 6: total_compressed_size = in.i64(ct); return true;
		}
		return false;
	});
}

// column_orders (field 7) is not kept: encode() writes TYPE_ORDER for every leaf column
void FileMetaData::decode(ThriftIn &in)
{
	*this = FileMetaData();
	in.read_struct([&](int16_t id, int ct) {
		switch (id) {
		case 1: version = in.i32(ct); return true;
		case 2: in.read_struct_list(ct, schema); return true;
		case 3: num_rows = in.i64(ct); return true;
		case 4: in.read_struct_list(ct, row_groups); return true;
		case 5:
			key_value_metadata.clear();
			in.read_list(ct, [&](int elem_type) {
				in.expect(elem_type, CT_STRUCT);
				pair<string, string> kv;
				in.read_struct([&](int16_t kv_id, int kv_ct) {
					if (kv_id == 1)
						kv.first = in.binary(kv_ct);
					else if (kv_id == 2)
						kv.second = in.binary(kv_ct);
					else
						return false;
					return true;
				});
				key_value_metadata.push_back(kv);
			});
			return true;
		case 6: created_by = in.binary(ct); return true;
		}
		return false;
	});
}

//------------------------------ encodings ------------------------------

// The bit-packed run of the RLE / bit-packed hybrid: values packed LSB first, in groups of 8;
// the last group is padded with zeros.
static void put_bit_packed(string &out, const uint32_t *v, uint64_t n, int bit_width)
{
	uint64_t num_groups = (n + 7) / 8;
	uint64_t acc = 0;
	int nbits = 0;

	put_varint(out, num_groups << 1 | 1);
	for (uint64_t i = 0; i < num_groups * 8; ++i) {
		acc |= (uint64_t)(i < n ? v[i] : 0) << nbits;
		for (nbits += bit_width; nbits >= 8; nbits -= 8) {
			out += (char)acc;
			acc >>= 8;
		}
	}
}

// RLE / bit-packed hybrid without the length prefix (as in RLE_DICTIONARY data pages).
// Repeats of 8 or more become RLE runs; values in between are bit-packed, the bit-packed run
// taking as many values of the following repeat as it needs to fill its last group of 8.
static void put_rle_hybrid(string &out, const vector<uint32_t> &v, int bit_width)
{
	uint64_t n = v.size();
	uint64_t pending = 0;  // first value not written yet

	for (uint64_t i = 0; i < n; ) {
		uint64_t j = i + 1;
		while (j < n && v[j] == v[i])
			++j;
		if (j - i >= 8) {
			uint64_t fill = (8 - (i - pending) % 8) % 8;
			if (i + fill > pending)
				put_bit_packed(out, &v[pending], i + fill - pending, bit_width);
			put_varint(out, (j - i - fill) << 1);
			for (int b = 0; b < bit_width; b += 8)
				out += (char)(v[i] >> b);
			pending = j;
		}
		i = j;
	}
	if (pending < n)
		put_bit_packed(out, &v[pending], n - pending, bit_width);
}

template <typename T> static string plain(T v) { return string((const char *)&v, sizeof(T)); }

//------------------------------ Arrow schema ------------------------------

// pyarrow keeps an Arrow schema in a Parquet file's key-value metadata, under this key: an Arrow
// IPC message (a FlatBuffer: https://github.com/apache/arrow/blob/main/format/Schema.fbs and
// Message.fbs) prefixed by 0xFFFFFFFF and its length, base64-encoded. Its fields carry the field
// metadata pyarrow shows as schema.field(name).metadata.
static const char *ARROW_SCHEMA_KEY = "ARROW:schema";

// Arrow's union ids and enums used here
enum { ARROW_MESSAGE_SCHEMA = 1 };                                // MessageHeader
enum { ARROW_TYPE_INT = 2, ARROW_TYPE_FLOAT = 3, ARROW_TYPE_UTF8 = 5 };  // Type
enum { ARROW_METADATA_V4 = 3 };                                   // MetadataVersion: read by pyarrow 0.15 and later
enum { ARROW_PRECISION_DOUBLE = 2 };

static const char BASE64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static string base64_encode(const string &in)
{
	string out;
	size_t i = 0;
	for (; i + 2 < in.size(); i += 3) {
		uint32_t v = (uint8_t)in[i] << 16 | (uint8_t)in[i + 1] << 8 | (uint8_t)in[i + 2];
		for (int k = 3; k >= 0; --k)
			out += BASE64[(v >> (6 * k)) & 63];
	}
	if (i + 1 == in.size()) {
		uint32_t v = (uint8_t)in[i] << 16;
		out += BASE64[v >> 18];
		out += BASE64[(v >> 12) & 63];
		out += "==";
	} else if (i + 2 == in.size()) {
		uint32_t v = (uint8_t)in[i] << 16 | (uint8_t)in[i + 1] << 8;
		out += BASE64[v >> 18];
		out += BASE64[(v >> 12) & 63];
		out += BASE64[(v >> 6) & 63];
		out += '=';
	}
	return out;
}

// false if in is not base64
static bool base64_decode(const string &in, string &out)
{
	uint32_t v = 0;
	int bits = 0;

	out.clear();
	for (char c : in) {
		if (c == '=' || c == '\n' || c == '\r')
			continue;
		const char *p = strchr(BASE64, c);
		if (!c || !p)
			return false;
		v = v << 6 | (uint32_t)(p - BASE64);
		bits += 6;
		if (bits >= 8) {
			bits -= 8;
			out += (char)((v >> bits) & 0xff);
		}
	}
	return true;
}

// Builds a FlatBuffer front to back: an object is written before the objects it refers to, so every
// offset points forward, and the referring field is set once its target is written (link()).
class FlatBufferBuilder {
public:
	struct Field {
		int      slot;   // the field's index in the table's vtable
		int      size;   // 1, 2, 4 or 8; 0 for an offset to another object
		uint64_t value;  // of a scalar
	};

	string buf;

	FlatBufferBuilder() : buf(4, '\0') {}  // the offset to the root table

	void align(size_t a) { while (buf.size() % a) buf += '\0'; }

	template <typename T> void set(size_t at, T v) { memcpy(&buf[at], &v, sizeof(T)); }

	// Sets the offset at position at to the object at target, which follows it
	void link(size_t at, size_t target) { set<uint32_t>(at, (uint32_t)(target - at)); }

	// Writes a table (its vtable first) and returns its position; offset_pos[slot] is the position
	// of each offset field, to link() once its object is written
	size_t table(const vector<Field> &fields, vector<size_t> &offset_pos)
	{
		int num_slots = 0;
		size_t max_align = 4;
		for (const Field &f : fields) {
			num_slots = max(num_slots, f.slot + 1);
			max_align = max(max_align, (size_t)(f.size ? f.size : 4));
		}

		align(2);
		size_t vtable = buf.size();
		buf.append(4 + 2 * num_slots, '\0');
		align(max_align);
		size_t table = buf.size();
		buf.append(4, '\0');
		set<int32_t>(table, (int32_t)(table - vtable));

		// the largest fields first, which needs no padding between them
		vector<Field> sorted(fields);
		stable_sort(sorted.begin(), sorted.end(), [](const Field &a, const Field &b) { return (a.size ? a.size : 4) > (b.size ? b.size : 4); });
		offset_pos.assign(num_slots, 0);
		for (const Field &f : sorted) {
			int size = f.size ? f.size : 4;
			align(size);
			size_t pos = buf.size();
			buf.append((char *)&f.value, size);  // little-endian: the low bytes of value
			set<uint16_t>(vtable + 4 + 2 * f.slot, (uint16_t)(pos - table));
			if (!f.size)
				offset_pos[f.slot] = pos;
		}
		set<uint16_t>(vtable, (uint16_t)(4 + 2 * num_slots));
		set<uint16_t>(vtable + 2, (uint16_t)(buf.size() - table));
		return table;
	}

	size_t add_string(const string &s)
	{
		align(4);
		size_t pos = buf.size();
		buf.append(4, '\0');
		set<uint32_t>(pos, (uint32_t)s.size());
		buf += s;
		buf += '\0';
		return pos;
	}

	// A vector of n offsets; elem_pos[i] is the position of element i, to link()
	size_t offset_vector(size_t n, vector<size_t> &elem_pos)
	{
		align(4);
		size_t pos = buf.size();
		buf.append(4 + 4 * n, '\0');
		set<uint32_t>(pos, (uint32_t)n);
		elem_pos.resize(n);
		for (size_t i = 0; i < n; ++i)
			elem_pos[i] = pos + 4 + 4 * i;
		return pos;
	}
};

// The ARROW:schema value of a file of the columns of schema (the root first; chrom a string,
// start and end INT32 or INT64, the rest DOUBLE), all nullable, with the field metadata of the
// value columns
static string arrow_schema(const vector<SchemaElement> &schema, const vector<vector<pair<string, string>>> &value_metadata)
{
	typedef FlatBufferBuilder::Field F;
	FlatBufferBuilder b;
	vector<size_t> pos, elems;

	size_t message = b.table({ F{ 0, 2, ARROW_METADATA_V4 }, F{ 1, 1, ARROW_MESSAGE_SCHEMA }, F{ 2, 0, 0 } }, pos);
	b.link(0, message);
	size_t header_pos = pos[2];
	size_t schema_table = b.table({ F{ 1, 0, 0 } }, pos);  // endianness: Little, the default
	b.link(header_pos, schema_table);
	size_t fields_pos = pos[1];
	vector<size_t> field_elems;
	b.link(fields_pos, b.offset_vector(schema.size() - 1, field_elems));

	for (size_t i = 1; i < schema.size(); ++i) {
		const SchemaElement &col = schema[i];
		const vector<pair<string, string>> *md = i >= 4 && i - 4 < value_metadata.size() && !value_metadata[i - 4].empty() ? &value_metadata[i - 4] : NULL;
		int type = col.type == Type::BYTE_ARRAY ? ARROW_TYPE_UTF8 : col.type == Type::DOUBLE ? ARROW_TYPE_FLOAT : ARROW_TYPE_INT;
		vector<F> fields = { F{ 0, 0, 0 }, F{ 1, 1, 1 }, F{ 2, 1, (uint64_t)type }, F{ 3, 0, 0 }, F{ 5, 0, 0 } };
		if (md)
			fields.push_back(F{ 6, 0, 0 });

		size_t field = b.table(fields, pos);
		vector<size_t> field_pos = pos;
		b.link(field_elems[i - 1], field);
		b.link(field_pos[0], b.add_string(col.name));

		if (type == ARROW_TYPE_UTF8)
			b.link(field_pos[3], b.table({}, pos));
		else if (type == ARROW_TYPE_FLOAT)
			b.link(field_pos[3], b.table({ F{ 0, 2, ARROW_PRECISION_DOUBLE } }, pos));
		else
			b.link(field_pos[3], b.table({ F{ 0, 4, (uint64_t)(col.type == Type::INT64 ? 64 : 32) }, F{ 1, 1, 1 } }, pos));

		b.link(field_pos[5], b.offset_vector(0, elems));  // children: none, but present (pyarrow < 4 requires it)

		if (md) {
			vector<size_t> kv_elems;
			b.link(field_pos[6], b.offset_vector(md->size(), kv_elems));
			for (size_t k = 0; k < md->size(); ++k) {
				size_t kv = b.table({ F{ 0, 0, 0 }, F{ 1, 0, 0 } }, pos);
				vector<size_t> kv_pos = pos;
				b.link(kv_elems[k], kv);
				b.link(kv_pos[0], b.add_string((*md)[k].first));
				b.link(kv_pos[1], b.add_string((*md)[k].second));
			}
		}
	}
	b.align(8);

	string message_bytes = "\xff\xff\xff\xff";
	uint32_t len = b.buf.size();
	message_bytes.append((const char *)&len, 4);
	message_bytes += b.buf;
	return base64_encode(message_bytes);
}

// Reads a FlatBuffer; any offset out of the buffer throws
class FlatBufferReader {
public:
	FlatBufferReader(const string &buf) : m_buf(buf) {}

	template <typename T> T get(size_t at) const
	{
		if (at > m_buf.size() || m_buf.size() - at < sizeof(T))
			throw out_of_range("flatbuffer");
		T v;
		memcpy(&v, m_buf.data() + at, sizeof(T));
		return v;
	}

	size_t deref(size_t at) const { return at + get<uint32_t>(at); }

	// The position of field slot of the table at table, 0 if absent
	size_t field(size_t table, int slot) const
	{
		size_t vtable = table - get<int32_t>(table);
		uint16_t vtable_size = get<uint16_t>(vtable);
		if (4 + 2 * slot >= vtable_size)
			return 0;
		uint16_t off = get<uint16_t>(vtable + 4 + 2 * slot);
		return off ? table + off : 0;
	}

	string str(size_t at) const
	{
		uint32_t n = get<uint32_t>(at);
		if (m_buf.size() - at - 4 < n)
			throw out_of_range("flatbuffer");
		return m_buf.substr(at + 4, n);
	}

	// The tables of a vector of offsets at at
	vector<size_t> tables(size_t at) const
	{
		vector<size_t> res;
		uint32_t n = get<uint32_t>(at);
		for (uint32_t i = 0; i < n; ++i)
			res.push_back(deref(at + 4 + 4 * i));
		return res;
	}

private:
	const string &m_buf;
};

map<string, vector<pair<string, string>>> arrow_field_metadata(const FileMetaData &meta)
{
	map<string, vector<pair<string, string>>> res;

	for (const auto &kv : meta.key_value_metadata) {
		if (kv.first != ARROW_SCHEMA_KEY)
			continue;

		string msg;
		if (!base64_decode(kv.second, msg))
			return res;

		try {
			// 0xFFFFFFFF and the length, or only the length (before Arrow 0.15)
			size_t start = msg.size() >= 4 && !memcmp(msg.data(), "\xff\xff\xff\xff", 4) ? 8 : 4;
			if (msg.size() < start)
				return res;
			string fb = msg.substr(start);
			FlatBufferReader r(fb);

			size_t message = r.deref(0);
			size_t header_type = r.field(message, 1);
			size_t header = r.field(message, 2);
			if (!header_type || r.get<uint8_t>(header_type) != ARROW_MESSAGE_SCHEMA || !header)
				return res;
			size_t schema = r.deref(header);
			size_t fields = r.field(schema, 1);
			if (!fields)
				return res;
			for (size_t field : r.tables(r.deref(fields))) {
				size_t name = r.field(field, 0);
				size_t md = r.field(field, 6);
				if (!name || !md)
					continue;
				vector<pair<string, string>> &out = res[r.str(r.deref(name))];
				for (size_t kv_table : r.tables(r.deref(md))) {
					size_t key = r.field(kv_table, 0);
					size_t value = r.field(kv_table, 1);
					if (key)
						out.push_back(make_pair(r.str(r.deref(key)), value ? r.str(r.deref(value)) : string()));
				}
			}
		} catch (const out_of_range &) {
			res.clear();
		}
		return res;
	}
	return res;
}

//------------------------------ ParquetWriter ------------------------------

// Column chunks with more distinct values (by bit pattern) are written PLAIN: arrow's 1 MB
// dictionary page limit, counted in 4-byte values
static const uint32_t MAX_DICTIONARY_SIZE = 262144;

const uint64_t ParquetWriter::MAX_ROW_GROUP_ROWS;

ParquetWriter::~ParquetWriter()
{
	if (m_fd >= 0)
		::close(m_fd);
}

void ParquetWriter::open(const string &filename, const vector<string> &value_names, bool int64_coords, const string &created_by, const string &name,
                         const vector<vector<pair<string, string>>> &value_metadata)
{
	const char *coords[] = { "chrom", "start", "end" };

	for (size_t i = 0; i < value_names.size(); ++i) {
		for (const char *coord : coords) {
			if (value_names[i] == coord)
				TGLError<ParquetWriter>(FORMAT_ERROR, "Parquet output cannot have a value column named \"%s\": chrom, start and end are the interval columns", coord);
		}
		for (size_t j = 0; j < i; ++j) {
			if (value_names[i] == value_names[j])
				TGLError<ParquetWriter>(FORMAT_ERROR, "Parquet output cannot have two columns named \"%s\"", value_names[i].c_str());
		}
	}

	m_meta = FileMetaData();
	m_meta.created_by = created_by;
	m_meta.schema.resize(4 + value_names.size());
	m_meta.schema[0].name = "schema";
	m_meta.schema[0].num_children = 3 + value_names.size();
	for (size_t i = 1; i < m_meta.schema.size(); ++i) {
		SchemaElement &col = m_meta.schema[i];
		col.repetition_type = FieldRepetitionType::OPTIONAL;  // nullable, as arrow writes them
		col.type = i == 1 ? Type::BYTE_ARRAY : i <= 3 ? (int64_coords ? Type::INT64 : Type::INT32) : Type::DOUBLE;
		col.name = i <= 3 ? coords[i - 1] : value_names[i - 4];
	}
	m_meta.schema[1].converted_type = ConvertedType::UTF8;
	m_meta.schema[1].logical_type = LogicalType::STRING;

	for (const auto &md : value_metadata) {
		if (!md.empty()) {
			m_meta.key_value_metadata.push_back(make_pair(string(ARROW_SCHEMA_KEY), arrow_schema(m_meta.schema, value_metadata)));
			break;
		}
	}

	m_int64_coords = int64_coords;
	m_chrom.clear();
	m_starts.clear();
	m_ends.clear();
	m_values.assign(value_names.size(), vector<uint64_t>());
	create(filename, name);
}

// Creates the file and writes the header
void ParquetWriter::create(const string &filename, const string &name)
{
	m_filename = name;
	m_fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
	if (m_fd < 0) {
		size_t slash = filename.rfind('/');
		string dir = slash == string::npos ? "." : slash ? filename.substr(0, slash) : "/";
		TGLError<ParquetWriter>(FILE_ERROR, "Failed to write %s: cannot create a temporary file in %s: %s", name.c_str(), dir.c_str(), strerror(errno));
	}
	m_offset = 0;
	write("PAR1");
}

void ParquetWriter::merge(const string &target, const vector<string> &parts, const string &name)
{
	ParquetWriter out;  // the row groups are copied as they are; out adds the header and the footer
	string schema;

	out.create(target, name);
	for (size_t i = 0; i < parts.size(); ++i) {
		{  // the reader is closed before the part is unlinked (NFS silly-renames open files)
			ParquetReader part;
			part.open(parts[i]);
			const FileMetaData &meta = part.metadata();

			string part_schema;
			for (const SchemaElement &se : meta.schema)
				se.encode(part_schema);
			if (!i) {
				schema = part_schema;
				out.m_meta = meta;
				out.m_meta.num_rows = 0;
				out.m_meta.row_groups.clear();
			} else if (part_schema != schema)
				TGLError<ParquetWriter>(FORMAT_ERROR, "Cannot merge %s into %s: its schema differs from that of %s", parts[i].c_str(), name.c_str(), parts[0].c_str());

			// A part without rows holds one empty row group (see close())
			bool skip = !meta.num_rows && (i + 1 < parts.size() || !out.m_meta.row_groups.empty());

			if (!skip) {
				int64_t shift = (int64_t)out.m_offset - 4;
				out.append(part.m_fd, parts[i], 4, part.m_footer_offset - 4);
				for (RowGroup rg : meta.row_groups) {
					rg.file_offset += shift;
					for (ColumnChunk &cc : rg.columns) {
						cc.file_offset += shift;
						cc.meta_data.data_page_offset += shift;
						if (cc.meta_data.dictionary_page_offset >= 0)
							cc.meta_data.dictionary_page_offset += shift;
					}
					out.m_meta.row_groups.push_back(rg);
				}
				out.m_meta.num_rows += meta.num_rows;
			}
		}
		// Frees the part's disk now, so that the peak is not every part plus the merged file;
		// the caller's guard unlinks the rest on failure and ignores the ENOENT of this one
		unlink(parts[i].c_str());
	}
	out.close();
}

void ParquetWriter::add_row(const string &chrom, int64_t start, int64_t end, const double *values)
{
	if (!m_int64_coords && (start < INT32_MIN || end > INT32_MAX))
		TGLError<ParquetWriter>(FORMAT_ERROR, "Interval (%s, %lld, %lld) does not fit 32-bit Parquet coordinates", chrom.c_str(), (long long)start, (long long)end);

	if (chrom != m_chrom || m_starts.size() >= MAX_ROW_GROUP_ROWS) {
		if (!m_starts.empty())
			flush_row_group();
		m_chrom = chrom;
	}

	m_starts.push_back(start);
	m_ends.push_back(end);
	for (size_t i = 0; i < m_values.size(); ++i) {
		uint64_t bits;
		memcpy(&bits, &values[i], sizeof(bits));
		m_values[i].push_back(bits);
	}
}

void ParquetWriter::close()
{
	// A file without rows still gets one, empty, row group, as pyarrow writes an empty table:
	// pyarrow's ParquetFile.read() and iter_batches() fail on a file with no row groups.
	if (!m_starts.empty() || m_meta.row_groups.empty())
		flush_row_group();

	string footer;
	m_meta.encode(footer);
	uint32_t footer_len = footer.size();
	footer += plain(footer_len);
	footer += "PAR1";
	write(footer);

	int fd = m_fd;
	m_fd = -1;
	if (::close(fd))
		TGLError<ParquetWriter>(FILE_ERROR, "Failed to close file %s: %s", m_filename.c_str(), strerror(errno));
}

void ParquetWriter::flush_row_group()
{
	RowGroup rg;
	rg.num_rows = m_starts.size();
	rg.file_offset = m_offset;
	rg.columns.resize(m_meta.schema.size() - 1);

	// chrom: a one-entry dictionary
	ColumnChunk &chrom = rg.columns[0];
	chrom.file_offset = m_offset;
	chrom.meta_data.type = Type::BYTE_ARRAY;
	chrom.meta_data.path_in_schema = { "chrom" };
	chrom.meta_data.num_values = rg.num_rows;
	write_dictionary_chunk(chrom.meta_data, plain((uint32_t)m_chrom.size()) + m_chrom, 1, vector<uint32_t>(m_starts.size(), 0));
	chrom.meta_data.has_statistics = rg.num_rows > 0;
	chrom.meta_data.statistics.min_value = chrom.meta_data.statistics.max_value = m_chrom;

	write_coord_column(rg.columns[1], "start", m_starts);
	write_coord_column(rg.columns[2], "end", m_ends);
	for (size_t i = 0; i < m_values.size(); ++i)
		write_fixed_column(rg.columns[i + 3], Type::DOUBLE, m_meta.schema[i + 4].name, m_values[i]);

	for (const ColumnChunk &cc : rg.columns) {
		rg.total_byte_size += cc.meta_data.total_uncompressed_size;
		rg.total_compressed_size += cc.meta_data.total_compressed_size;
	}
	m_meta.num_rows += rg.num_rows;
	m_meta.row_groups.push_back(rg);

	m_starts.clear();
	m_ends.clear();
	for (vector<uint64_t> &v : m_values)
		v.clear();
}

void ParquetWriter::write_coord_column(ColumnChunk &cc, const string &name, const vector<int64_t> &coords)
{
	if (m_int64_coords)
		write_fixed_column(cc, Type::INT64, name, coords);
	else
		write_fixed_column(cc, Type::INT32, name, vector<int32_t>(coords.begin(), coords.end()));

	if (coords.empty())  // no min and max to write
		return;

	auto minmax = minmax_element(coords.begin(), coords.end());
	Statistics &stats = cc.meta_data.statistics;
	stats.min_value = m_int64_coords ? plain(*minmax.first) : plain((int32_t)*minmax.first);
	stats.max_value = m_int64_coords ? plain(*minmax.second) : plain((int32_t)*minmax.second);
	cc.meta_data.has_statistics = true;
}

// Dictionary-encoded unless the chunk has more than MAX_DICTIONARY_SIZE distinct values. The
// dictionary is keyed on the bit pattern, so NaN (NaN != NaN) still takes one entry.
template <typename T>
void ParquetWriter::write_fixed_column(ColumnChunk &cc, int type, const string &name, const vector<T> &vals)
{
	ColumnMetaData &md = cc.meta_data;
	cc.file_offset = m_offset;
	md.type = type;
	md.path_in_schema = { name };
	md.num_values = vals.size();

	unordered_map<T, uint32_t> ids;
	vector<uint32_t> indices(vals.size());
	string dict;
	size_t i;

	for (i = 0; i < vals.size(); ++i) {
		auto res = ids.emplace(vals[i], (uint32_t)ids.size());
		if (res.second) {
			if (ids.size() > MAX_DICTIONARY_SIZE)
				break;
			dict += plain(vals[i]);
		}
		indices[i] = res.first->second;
	}

	if (i == vals.size()) {
		write_dictionary_chunk(md, dict, ids.size(), indices);
		return;
	}

	PageHeader header;
	header.data_page_header.num_values = vals.size();
	header.data_page_header.encoding = Encoding::PLAIN;
	md.data_page_offset = m_offset;
	md.codec = CompressionCodec::GZIP;
	md.total_uncompressed_size = write_page(header, string((const char *)vals.data(), vals.size() * sizeof(T)));
	md.encodings = { Encoding::PLAIN, Encoding::RLE };
	md.total_compressed_size = m_offset - cc.file_offset;
}

// A PLAIN dictionary page followed by one RLE_DICTIONARY data page
void ParquetWriter::write_dictionary_chunk(ColumnMetaData &md, const string &dict, uint32_t dict_size, const vector<uint32_t> &indices)
{
	uint64_t start = m_offset;
	PageHeader header;

	header.type = PageType::DICTIONARY_PAGE;
	header.dictionary_page_header.num_values = dict_size;
	header.dictionary_page_header.encoding = Encoding::PLAIN;
	md.dictionary_page_offset = m_offset;
	md.codec = CompressionCodec::GZIP;
	md.total_uncompressed_size = write_page(header, dict);

	int bit_width = 1;
	while (((uint64_t)1 << bit_width) < dict_size)
		++bit_width;
	string data(1, (char)bit_width);
	put_rle_hybrid(data, indices, bit_width);

	header = PageHeader();
	header.data_page_header.num_values = indices.size();
	header.data_page_header.encoding = Encoding::RLE_DICTIONARY;
	md.data_page_offset = m_offset;
	md.total_uncompressed_size += write_page(header, data);

	md.encodings = { Encoding::PLAIN, Encoding::RLE, Encoding::RLE_DICTIONARY };
	md.total_compressed_size = m_offset - start;
}

// One gzip stream: the GZIP codec of Parquet (RFC 1952, not zlib or raw deflate)
static string gzip_compress(const string &in)
{
	z_stream z;
	memset(&z, 0, sizeof(z));
	if (deflateInit2(&z, 6, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		TGLError<ParquetWriter>(ParquetWriter::FORMAT_ERROR, "zlib initialization failed");
	string out(deflateBound(&z, in.size()), '\0');  // enough for deflate() to finish in one call
	z.next_in = (Bytef *)in.data();
	z.avail_in = in.size();
	z.next_out = (Bytef *)&out[0];
	z.avail_out = out.size();
	int res = deflate(&z, Z_FINISH);
	out.resize(z.total_out);
	deflateEnd(&z);
	if (res != Z_STREAM_END)
		TGLError<ParquetWriter>(ParquetWriter::FORMAT_ERROR, "zlib compression failed");
	return out;
}

// Writes the page GZIP-compressed and returns its uncompressed size, header included. A data
// page (v1) starts with the definition levels of its OPTIONAL column, all 1 as misha writes no
// nulls: the 4-byte length, then one RLE run (bit width 1) of num_values ones; no run when
// there are no values, as a run cannot be empty.
uint64_t ParquetWriter::write_page(PageHeader &header, const string &payload)
{
	string page;
	if (header.type == PageType::DATA_PAGE) {
		string levels;
		if (header.data_page_header.num_values) {
			put_varint(levels, (uint64_t)header.data_page_header.num_values << 1);
			levels += '\1';
		}
		page = plain((uint32_t)levels.size()) + levels;
	}
	page += payload;

	string compressed = gzip_compress(page);
	header.uncompressed_page_size = page.size();
	header.compressed_page_size = compressed.size();
	string bytes;
	header.encode(bytes);
	write(bytes);
	write(compressed);
	return bytes.size() + page.size();
}

void ParquetWriter::write(const string &bytes)
{
	const char *p = bytes.data();
	size_t left = bytes.size();

	while (left) {
		ssize_t n = ::write(m_fd, p, left);
		if (n < 0) {
			if (errno == EINTR)
				continue;
			TGLError<ParquetWriter>(FILE_ERROR, "Failed to write to file %s: %s", m_filename.c_str(), strerror(errno));
		}
		p += n;
		left -= n;
	}
	m_offset += bytes.size();
}

// Appends size bytes of fd (the file filename) from offset: with copy_file_range() where the
// kernel can copy between the two files, through a buffer where it cannot
void ParquetWriter::append(int fd, const string &filename, uint64_t offset, uint64_t size)
{
	uint64_t end = offset + size;
	bool kernel_copy = true;  // until copy_file_range() turns out not to work here
	string buf;

	while (offset < end) {
		rdb::check_interrupt();

		size_t chunk = min<uint64_t>(end - offset, 1 << 26);
		ssize_t n;

		if (kernel_copy) {
#if defined(__linux__) && defined(SYS_copy_file_range)
			// the system call: glibc < 2.27 has no copy_file_range()
			loff_t off = offset;
			n = syscall(SYS_copy_file_range, fd, &off, m_fd, NULL, chunk, 0);
#else
			n = -1;
			errno = ENOSYS;
#endif
			// not supported for these files (EINVAL is how some filesystems say so)
			if (n < 0 && (errno == ENOSYS || errno == EXDEV || errno == EOPNOTSUPP || errno == EINVAL)) {
				kernel_copy = false;
				continue;
			}
			if (n > 0)
				m_offset += n;
		} else {
			buf.resize(chunk);
			n = pread(fd, &buf[0], chunk, offset);
			if (n > 0) {
				buf.resize(n);
				write(buf);
			}
		}

		if (n < 0 && errno == EINTR)
			continue;
		if (n <= 0)
			TGLError<ParquetWriter>(FILE_ERROR, "Failed to copy %s to %s: %s", filename.c_str(), m_filename.c_str(), n ? strerror(errno) : "unexpected end of file");
		offset += n;
	}
}

//------------------------------ ParquetReader ------------------------------

static const char *type_name(int type)
{
	static const char *names[] = { "BOOLEAN", "INT32", "INT64", "INT96", "FLOAT", "DOUBLE", "BYTE_ARRAY", "FIXED_LEN_BYTE_ARRAY" };
	return type >= 0 && type < (int)(sizeof(names) / sizeof(names[0])) ? names[type] : "unknown";
}

static const char *encoding_name(int encoding)
{
	static const char *names[] = { "PLAIN", "GROUP_VAR_INT", "PLAIN_DICTIONARY", "RLE", "BIT_PACKED", "DELTA_BINARY_PACKED",
	                               "DELTA_LENGTH_BYTE_ARRAY", "DELTA_BYTE_ARRAY", "RLE_DICTIONARY", "BYTE_STREAM_SPLIT", "ALP" };
	return encoding >= 0 && encoding < (int)(sizeof(names) / sizeof(names[0])) ? names[encoding] : "unknown";
}

static const char *codec_name(int codec)
{
	static const char *names[] = { "UNCOMPRESSED", "SNAPPY", "GZIP", "LZO", "BROTLI", "LZ4", "ZSTD", "LZ4_RAW" };
	return codec >= 0 && codec < (int)(sizeof(names) / sizeof(names[0])) ? names[codec] : "unknown";
}

// Sets col.unsigned_int and col.other_type from the converted and logical types of se. Kept as
// numbers or strings: no annotation, INT_* / UINT_* / INTEGER, UTF8 / STRING, ENUM, JSON;
// anything else (DECIMAL, DATE, TIME, TIMESTAMP, ...) is other_type, refused rather than misread.
static void read_annotations(const SchemaElement &se, ParquetColumn &col)
{
	static const char *converted[] = { "UTF8", "MAP", "MAP_KEY_VALUE", "LIST", "ENUM", "DECIMAL", "DATE", "TIME_MILLIS", "TIME_MICROS",
	                                   "TIMESTAMP_MILLIS", "TIMESTAMP_MICROS", "UINT_8", "UINT_16", "UINT_32", "UINT_64", "INT_8", "INT_16",
	                                   "INT_32", "INT_64", "JSON", "BSON", "INTERVAL" };
	static const char *logical[] = { "", "STRING", "MAP", "LIST", "ENUM", "DECIMAL", "DATE", "TIME", "TIMESTAMP", "", "INTEGER", "UNKNOWN",
	                                 "JSON", "BSON", "UUID", "FLOAT16", "VARIANT", "GEOMETRY", "GEOGRAPHY", "FILE" };
	int ct = se.converted_type;
	int lt = se.logical_type;
	bool ct_int = ct >= ConvertedType::UINT_8 && ct <= ConvertedType::INT_64;
	bool ct_kept = ct < 0 || ct_int || ct == ConvertedType::UTF8 || ct == ConvertedType::ENUM || ct == ConvertedType::JSON;
	bool lt_kept = lt < 0 || lt == LogicalType::INTEGER || lt == LogicalType::STRING || lt == LogicalType::ENUM || lt == LogicalType::JSON;

	col.unsigned_int = (ct_int && ct <= ConvertedType::UINT_64) || (lt == LogicalType::INTEGER && !se.integer_signed);
	col.other_type.clear();
	if (!ct_kept)
		col.other_type = ct < (int)(sizeof(converted) / sizeof(converted[0])) ? converted[ct] : "an unknown converted type";
	else if (!lt_kept)
		col.other_type = lt < (int)(sizeof(logical) / sizeof(logical[0])) && *logical[lt] ? logical[lt] : "an unknown logical type";
}

[[noreturn]] static void unsupported_encoding(int encoding)
{
	TGLError<ParquetReader>(ParquetReader::UNSUPPORTED, "encoding %s (%d) is not supported", encoding_name(encoding), encoding);
}

static uint32_t get_le32(const uint8_t *p)
{
	uint32_t v;
	memcpy(&v, p, sizeof(v));
	return v;
}

// Snappy: a varint preamble with the uncompressed length, then literals and back-references.
// 'size' is the length the page header gives.
static void snappy_decompress(const uint8_t *src, size_t n, string &out, uint64_t size)
{
	const uint8_t *end = src + n;
	uint64_t len = 0;

	for (int shift = 0; ; shift += 7) {
		if (src == end || shift > 28)
			corrupt("snappy preamble");
		uint8_t b = *src++;
		len |= (uint64_t)(b & 0x7f) << shift;
		if (!(b & 0x80))
			break;
	}
	if (len != size)
		corrupt("snappy length differs from the page size");

	out.resize(len);
	char *dst = &out[0];
	uint64_t pos = 0;

	while (src < end) {
		uint8_t tag = *src++;
		uint64_t length, offset;

		if ((tag & 3) == 0) {  // literal
			length = tag >> 2;
			if (length >= 60) {
				int nbytes = length - 59;
				if (end - src < nbytes)
					corrupt("snappy literal");
				length = 0;
				for (int i = 0; i < nbytes; ++i)
					length |= (uint64_t)src[i] << (8 * i);
				src += nbytes;
			}
			++length;
			if ((uint64_t)(end - src) < length || len - pos < length)
				corrupt("snappy literal");
			memcpy(dst + pos, src, length);
			src += length;
			pos += length;
			continue;
		}

		if ((tag & 3) == 1) {
			if (src == end)
				corrupt("snappy copy");
			length = ((tag >> 2) & 7) + 4;
			offset = (uint64_t)(tag >> 5) << 8 | *src++;
		} else {
			int nbytes = (tag & 3) == 2 ? 2 : 4;
			if (end - src < nbytes)
				corrupt("snappy copy");
			length = (tag >> 2) + 1;
			offset = 0;
			for (int i = 0; i < nbytes; ++i)
				offset |= (uint64_t)src[i] << (8 * i);
			src += nbytes;
		}
		if (!offset || offset > pos || len - pos < length)
			corrupt("snappy copy");
		for (uint64_t i = 0; i < length; ++i, ++pos)  // byte by byte: the source may overlap the copy
			dst[pos] = dst[pos - offset];
	}

	if (pos != len)
		corrupt("snappy length");
}

static void gzip_decompress(const uint8_t *src, size_t n, string &out, size_t size)
{
	z_stream z;
	memset(&z, 0, sizeof(z));
	if (inflateInit2(&z, 16 + MAX_WBITS) != Z_OK)
		TGLError<ParquetReader>(ParquetReader::FORMAT_ERROR, "zlib initialization failed");

	out.resize(size);
	z.next_in = (Bytef *)src;
	z.avail_in = n;
	z.next_out = (Bytef *)&out[0];
	z.avail_out = size;
	int res = inflate(&z, Z_FINISH);
	uint64_t total = z.total_out;
	inflateEnd(&z);
	if (res != Z_STREAM_END || total != size)
		corrupt("gzip");
}

// Decompresses [page, page_end) into buf, which must come out at 'size' bytes, and points page
// and page_end at it
static void decompress(int codec, const uint8_t *&page, const uint8_t *&page_end, int64_t size, string &buf)
{
	if (codec == CompressionCodec::UNCOMPRESSED)
		return;
	if (size < 0)
		corrupt("page size");
	if (codec == CompressionCodec::SNAPPY)
		snappy_decompress(page, page_end - page, buf, size);
	else
		gzip_decompress(page, page_end - page, buf, size);
	page = (const uint8_t *)buf.data();
	page_end = page + buf.size();
}

// Appends n values of the RLE / bit-packed hybrid (without a length prefix) to out
static void get_rle_hybrid(const uint8_t *p, const uint8_t *end, int bit_width, uint64_t n, vector<uint32_t> &out)
{
	if (bit_width < 0 || bit_width > 32)
		corrupt("bit width");

	uint64_t target = out.size() + n;
	uint64_t mask = ((uint64_t)1 << bit_width) - 1;
	int value_bytes = (bit_width + 7) / 8;

	while (out.size() < target) {
		uint64_t header = get_varint(p, end);
		uint64_t count = header >> 1;

		if (header & 1) {  // bit-packed: count groups of 8 values
			uint64_t avail = end - p;
			uint64_t take = min(count, (target - out.size() + 7) / 8) * 8;
			take = min(take, target - out.size());
			if ((take * bit_width + 7) / 8 > avail)
				corrupt("bit-packed run");
			uint64_t acc = 0;
			int nbits = 0;
			const uint8_t *q = p;
			for (uint64_t i = 0; i < take; ++i) {
				for (; nbits < bit_width; nbits += 8)
					acc |= (uint64_t)*q++ << nbits;
				out.push_back(acc & mask);
				acc >>= bit_width;
				nbits -= bit_width;
			}
			p += bit_width && count > avail ? avail : min(count * bit_width, avail);  // the last run may be cut short
		} else {  // RLE: count repeats of one value
			if (end - p < value_bytes)
				corrupt("RLE run");
			uint32_t v = 0;
			for (int i = 0; i < value_bytes; ++i)
				v |= (uint32_t)p[i] << (8 * i);
			if (v > mask)  // a definition level of 2 would count as two values
				corrupt("RLE run value exceeds the bit width");
			p += value_bytes;
			out.insert(out.end(), min(count, target - out.size()), v);
		}
	}
}

template <typename T> static void get_plain(const uint8_t *p, const uint8_t *end, uint64_t n, vector<T> &out)
{
	if ((uint64_t)(end - p) / sizeof(T) < n)
		corrupt("PLAIN values");
	if (!n)  // an all-null page: &out[old] would be past the end
		return;
	size_t old = out.size();
	out.resize(old + n);
	memcpy(&out[old], p, n * sizeof(T));
}

// Longer strings are refused: a dictionary entry is copied into every row that refers to it, so a
// small file of one long entry could otherwise ask for terabytes. Chromosome names are far shorter.
static const uint32_t MAX_STRING_BYTES = 1024;

template <> void get_plain<string>(const uint8_t *p, const uint8_t *end, uint64_t n, vector<string> &out)
{
	for (uint64_t i = 0; i < n; ++i) {
		if (end - p < 4)
			corrupt("PLAIN values");
		uint32_t len = get_le32(p);
		p += 4;
		if ((uint64_t)(end - p) < len)
			corrupt("PLAIN values");
		if (len > MAX_STRING_BYTES)
			TGLError<ParquetReader>(ParquetReader::UNSUPPORTED, "a string of %u bytes; strings longer than %u bytes are not supported", len, MAX_STRING_BYTES);
		out.emplace_back((const char *)p, len);
		p += len;
	}
}

// Byte k of value i is at p[k * n + i]
template <typename T> static void get_byte_stream_split(const uint8_t *p, const uint8_t *end, uint64_t n, vector<T> &out)
{
	if ((uint64_t)(end - p) / sizeof(T) < n)
		corrupt("BYTE_STREAM_SPLIT values");
	for (uint64_t i = 0; i < n; ++i) {
		uint8_t bytes[sizeof(T)];
		for (size_t k = 0; k < sizeof(T); ++k)
			bytes[k] = p[k * n + i];
		T v;
		memcpy(&v, bytes, sizeof(T));
		out.push_back(v);
	}
}

template <> void get_byte_stream_split<string>(const uint8_t *, const uint8_t *, uint64_t, vector<string> &)
{
	unsupported_encoding(Encoding::BYTE_STREAM_SPLIT);
}

// Appends the n values of a data page to out
template <typename T>
static void get_values(int encoding, const uint8_t *p, const uint8_t *end, uint64_t n, const vector<T> *dict, vector<T> &out, vector<uint32_t> &scratch)
{
	switch (encoding) {
	case Encoding::PLAIN:
		get_plain(p, end, n, out);
		return;
	case Encoding::PLAIN_DICTIONARY:
	case Encoding::RLE_DICTIONARY:
		if (!dict)
			corrupt("dictionary-encoded page without a dictionary");
		if (!n)
			return;
		if (p == end)
			corrupt("dictionary indices");
		scratch.clear();
		get_rle_hybrid(p + 1, end, *p, n, scratch);
		for (uint32_t idx : scratch) {
			if (idx >= dict->size())
				corrupt("dictionary index out of range");
			out.push_back((*dict)[idx]);
		}
		return;
	case Encoding::BYTE_STREAM_SPLIT:
		get_byte_stream_split(p, end, n, out);
		return;
	default:
		unsupported_encoding(encoding);
	}
}

ParquetReader::~ParquetReader()
{
	if (m_fd >= 0)
		::close(m_fd);
}

void ParquetReader::pread_exact(void *buf, uint64_t size, uint64_t offset)
{
	char *p = (char *)buf;

	while (size) {
		ssize_t n = ::pread(m_fd, p, size, offset);
		if (n < 0) {
			if (errno == EINTR)
				continue;
			TGLError<ParquetReader>(FILE_ERROR, "Failed to read file %s: %s", m_filename.c_str(), strerror(errno));
		}
		if (!n)
			TGLError<ParquetReader>(FILE_ERROR, "Failed to read file %s: unexpected end of file", m_filename.c_str());
		p += n;
		size -= n;
		offset += n;
	}
}

void ParquetReader::open(const string &filename, const string &name)
{
	if (m_fd >= 0)
		::close(m_fd);
	m_fd = -1;
	m_filename = name.empty() ? filename : name;
	m_columns.clear();

	m_fd = ::open(filename.c_str(), O_RDONLY | O_CLOEXEC);
	if (m_fd < 0)
		TGLError<ParquetReader>(FILE_ERROR, "Failed to open file %s: %s", m_filename.c_str(), strerror(errno));

	struct stat st;
	if (fstat(m_fd, &st))
		TGLError<ParquetReader>(FILE_ERROR, "Failed to stat file %s: %s", m_filename.c_str(), strerror(errno));
	m_size = st.st_size;

	char head[4], tail[8];
	if (m_size >= 12) {
		pread_exact(head, 4, 0);
		pread_exact(tail, 8, m_size - 8);
	}
	if (m_size >= 12 && !memcmp(tail + 4, "PARE", 4))
		TGLError<ParquetReader>(UNSUPPORTED, "File %s: encrypted Parquet files are not supported", m_filename.c_str());
	if (m_size < 12 || memcmp(head, "PAR1", 4) || memcmp(tail + 4, "PAR1", 4))
		TGLError<ParquetReader>(FORMAT_ERROR, "File %s is not a Parquet file, or is truncated", m_filename.c_str());

	uint32_t footer_len = get_le32((const uint8_t *)tail);
	if (footer_len > m_size - 12)
		TGLError<ParquetReader>(FORMAT_ERROR, "File %s: invalid footer length %u", m_filename.c_str(), footer_len);

	m_footer_offset = m_size - 8 - footer_len;
	string footer(footer_len, '\0');
	pread_exact(&footer[0], footer_len, m_footer_offset);

	try {
		ThriftIn in(footer.data(), footer.size());
		m_meta.decode(in);
	} catch (TGLException &e) {
		TGLError<ParquetReader>(FORMAT_ERROR, "File %s: %s", m_filename.c_str(), e.msg());
	}

	const vector<SchemaElement> &schema = m_meta.schema;
	if (schema.empty() || schema[0].num_children < 0)
		TGLError<ParquetReader>(FORMAT_ERROR, "File %s: invalid schema", m_filename.c_str());
	if ((size_t)schema[0].num_children != schema.size() - 1)
		TGLError<ParquetReader>(UNSUPPORTED, "File %s: only flat schemas (a root group of plain columns) are supported", m_filename.c_str());

	for (size_t i = 1; i < schema.size(); ++i) {
		if (schema[i].num_children > 0)
			TGLError<ParquetReader>(UNSUPPORTED, "File %s: column %s is nested; only flat schemas (a root group of plain columns) are supported",
			                        m_filename.c_str(), schema[i].name.c_str());
		m_columns.push_back(ParquetColumn{ schema[i].name, schema[i].type, schema[i].repetition_type == FieldRepetitionType::OPTIONAL, false, "" });
		read_annotations(schema[i], m_columns.back());
	}

	for (size_t i = 0; i < m_meta.row_groups.size(); ++i) {
		const RowGroup &rg = m_meta.row_groups[i];
		if (rg.columns.size() != m_columns.size() || rg.num_rows < 0)
			TGLError<ParquetReader>(FORMAT_ERROR, "File %s, row group %llu: invalid metadata", m_filename.c_str(), (unsigned long long)i);
	}
}

int ParquetReader::column_index(const string &name) const
{
	for (size_t i = 0; i < m_columns.size(); ++i) {
		if (m_columns[i].name == name)
			return (int)i;
	}
	return -1;
}

template <typename T>
void ParquetReader::read_chunk(int row_group, int col, vector<T> &vals, vector<uint8_t> &defs)
{
	vals.clear();
	defs.clear();

	try {
		const RowGroup &rg = m_meta.row_groups.at(row_group);
		const ColumnChunk &cc = rg.columns.at(col);
		const ColumnMetaData &md = cc.meta_data;
		bool optional = m_columns[col].optional;

		if (m_meta.schema[col + 1].repetition_type == FieldRepetitionType::REPEATED)
			TGLError<ParquetReader>(UNSUPPORTED, "repeated columns are not supported");
		if (!cc.file_path.empty())
			TGLError<ParquetReader>(UNSUPPORTED, "the column chunk is in another file (%s)", cc.file_path.c_str());
		if (md.path_in_schema.size() != 1 || md.path_in_schema[0] != m_columns[col].name || md.type != m_columns[col].physical_type)
			corrupt("column chunk metadata does not match the schema");
		if (md.codec != CompressionCodec::UNCOMPRESSED && md.codec != CompressionCodec::SNAPPY && md.codec != CompressionCodec::GZIP)
			TGLError<ParquetReader>(UNSUPPORTED, "compression codec %s (%d) is not supported; supported are UNCOMPRESSED, SNAPPY and GZIP (rewrite the file with compression=\"snappy\")",
			                        codec_name(md.codec), md.codec);

		// nothing to read in a chunk of no values (pyarrow writes data_page_offset 0 for one)
		if (md.num_values == 0 && rg.num_rows == 0)
			return;

		int64_t start = md.dictionary_page_offset > 0 && md.dictionary_page_offset < md.data_page_offset ? md.dictionary_page_offset : md.data_page_offset;
		if (start < 4 || md.total_compressed_size < 0 || (uint64_t)start + md.total_compressed_size > m_size - 8)
			corrupt("column chunk offsets");

		string chunk(md.total_compressed_size, '\0');
		pread_exact(&chunk[0], chunk.size(), start);

		const uint8_t *p = (const uint8_t *)chunk.data();
		const uint8_t *chunk_end = p + chunk.size();
		vector<T> dict;
		bool has_dict = false;
		string buf;
		vector<uint32_t> levels;
		vector<uint32_t> scratch;
		int64_t num_values = 0;

		if (md.num_values != rg.num_rows)  // before anything is allocated from a page's counts
			corrupt("number of values differs from the number of rows");

		while (num_values < md.num_values) {
			if (p == chunk_end)
				corrupt("column chunk ends before its values");
			PageHeader ph;
			ThriftIn in(p, chunk_end - p);
			ph.decode(in);
			p = in.p;
			if (ph.compressed_page_size < 0 || ph.uncompressed_page_size < 0 || chunk_end - p < ph.compressed_page_size)
				corrupt("page size");
			const uint8_t *page = p;
			const uint8_t *page_end = p + ph.compressed_page_size;
			p = page_end;

			if (ph.type == PageType::DICTIONARY_PAGE) {
				int encoding = ph.dictionary_page_header.encoding;
				if (encoding != Encoding::PLAIN && encoding != Encoding::PLAIN_DICTIONARY)
					unsupported_encoding(encoding);
				decompress(md.codec, page, page_end, ph.uncompressed_page_size, buf);
				if (ph.dictionary_page_header.num_values < 0)
					corrupt("dictionary size");
				dict.clear();
				get_plain(page, page_end, ph.dictionary_page_header.num_values, dict);
				has_dict = true;
				continue;
			}

			if (ph.type != PageType::DATA_PAGE && ph.type != PageType::DATA_PAGE_V2)
				continue;  // an index page, or a page type from the future

			bool v2 = ph.type == PageType::DATA_PAGE_V2;
			int64_t n = v2 ? ph.data_page_header_v2.num_values : ph.data_page_header.num_values;
			int encoding = v2 ? ph.data_page_header_v2.encoding : ph.data_page_header.encoding;
			const uint8_t *levels_begin, *levels_end;  // the definition levels

			if (n < 0 || n > md.num_values - num_values)
				corrupt("number of values");

			if (v2) {
				// levels are never compressed; the values are, if is_compressed
				const DataPageHeaderV2 &h = ph.data_page_header_v2;
				if (h.repetition_levels_byte_length < 0 || h.definition_levels_byte_length < 0 ||
				    (int64_t)h.repetition_levels_byte_length + h.definition_levels_byte_length > page_end - page)
					corrupt("level lengths");
				levels_begin = page + h.repetition_levels_byte_length;
				levels_end = page = levels_begin + h.definition_levels_byte_length;
				if (h.is_compressed)
					decompress(md.codec, page, page_end, (int64_t)ph.uncompressed_page_size - h.repetition_levels_byte_length - h.definition_levels_byte_length, buf);
			} else {
				decompress(md.codec, page, page_end, ph.uncompressed_page_size, buf);
				// a flat column has no repetition levels; definition levels, if OPTIONAL, are
				// RLE with a 4-byte length prefix
				levels_begin = levels_end = page;
				if (optional) {
					if (ph.data_page_header.definition_level_encoding != Encoding::RLE)
						unsupported_encoding(ph.data_page_header.definition_level_encoding);
					if (page_end - page < 4 || get_le32(page) > (uint64_t)(page_end - page - 4))
						corrupt("definition levels");
					levels_begin = page + 4;
					levels_end = page = levels_begin + get_le32(page);
				}
			}

			uint64_t num_present = n;
			if (optional) {
				levels.clear();
				get_rle_hybrid(levels_begin, levels_end, 1, n, levels);
				num_present = 0;
				for (uint32_t l : levels) {
					defs.push_back(l);
					num_present += l;
				}
			}

			get_values(encoding, page, page_end, num_present, has_dict ? &dict : NULL, vals, scratch);
			num_values += n;
		}
	} catch (TGLException &e) {
		TGLError<ParquetReader>(e.code(), "File %s, row group %d, column %s: %s", m_filename.c_str(), row_group, m_columns[col].name.c_str(), e.msg());
	} catch (const out_of_range &) {
		TGLError<ParquetReader>(FORMAT_ERROR, "File %s: no row group %d or column %d", m_filename.c_str(), row_group, col);
	}
}

template <typename T, typename O>
void ParquetReader::read_as(int row_group, int col, vector<O> &out, const O *null_value)
{
	vector<T> vals;
	vector<uint8_t> defs;

	read_chunk(row_group, col, vals, defs);
	out.clear();
	if (defs.empty()) {
		out.assign(vals.begin(), vals.end());
		return;
	}

	out.reserve(defs.size());
	size_t j = 0;
	for (uint8_t d : defs) {
		if (d)
			out.push_back((O)vals[j++]);
		else if (null_value)
			out.push_back(*null_value);
		else
			TGLError<ParquetReader>(FORMAT_ERROR, "File %s, row group %d, column %s: null values are not allowed in this column",
			                        m_filename.c_str(), row_group, m_columns[col].name.c_str());
	}
}

[[noreturn]] static void wrong_type(const string &filename, int row_group, const ParquetColumn &col, const char *expected)
{
	TGLError<ParquetReader>(ParquetReader::UNSUPPORTED, "File %s, row group %d, column %s: type %s%s%s is not supported here; expected %s",
	                        filename.c_str(), row_group, col.name.c_str(), type_name(col.physical_type),
	                        col.other_type.empty() ? "" : " annotated ", col.other_type.c_str(), expected);
}

void ParquetReader::read_doubles(int row_group, int col, vector<double> &out)
{
	const double nan = numeric_limits<double>::quiet_NaN();
	const ParquetColumn &c = m_columns.at(col);

	if (c.other_type.empty()) {
		switch (c.physical_type) {
		case Type::INT32:
			c.unsigned_int ? read_as<uint32_t>(row_group, col, out, &nan) : read_as<int32_t>(row_group, col, out, &nan);
			return;
		case Type::INT64:
			c.unsigned_int ? read_as<uint64_t>(row_group, col, out, &nan) : read_as<int64_t>(row_group, col, out, &nan);
			return;
		case Type::FLOAT: read_as<float>(row_group, col, out, &nan); return;
		case Type::DOUBLE: read_as<double>(row_group, col, out, &nan); return;
		}
	}
	wrong_type(m_filename, row_group, c, "a number (INT32, INT64, FLOAT or DOUBLE)");
}

void ParquetReader::read_int64s(int row_group, int col, vector<int64_t> &out)
{
	const ParquetColumn &c = m_columns.at(col);

	if (c.other_type.empty()) {
		switch (c.physical_type) {
		case Type::INT32:
			c.unsigned_int ? read_as<uint32_t>(row_group, col, out, (const int64_t *)NULL) : read_as<int32_t>(row_group, col, out, (const int64_t *)NULL);
			return;
		case Type::INT64:
			if (!c.unsigned_int) {
				read_as<int64_t>(row_group, col, out, (const int64_t *)NULL);
				return;
			}
			read_as<uint64_t>(row_group, col, out, (const int64_t *)NULL);
			for (int64_t v : out) {
				if (v < 0)  // an unsigned value past INT64_MAX
					TGLError<ParquetReader>(FORMAT_ERROR, "File %s, row group %d, column %s: value %llu is too large",
					                        m_filename.c_str(), row_group, c.name.c_str(), (unsigned long long)v);
			}
			return;
		case Type::FLOAT:
		case Type::DOUBLE: {
			vector<double> vals;

			c.physical_type == Type::FLOAT ? read_as<float>(row_group, col, vals, (const double *)NULL) : read_as<double>(row_group, col, vals, (const double *)NULL);
			out.resize(vals.size());
			for (size_t i = 0; i < vals.size(); ++i) {
				// [-2^63, 2^63) is where the cast is defined; NaN fails the comparisons
				if (!(vals[i] >= -9223372036854775808.0 && vals[i] < 9223372036854775808.0) || vals[i] != floor(vals[i]))
					TGLError<ParquetReader>(FORMAT_ERROR, "File %s, row group %d, column %s: value %.17g is not an integer",
					                        m_filename.c_str(), row_group, c.name.c_str(), vals[i]);
				out[i] = (int64_t)vals[i];
			}
			return;
		}
		}
	}
	wrong_type(m_filename, row_group, c, "an integer (INT32 or INT64), or FLOAT or DOUBLE holding whole numbers");
}

void ParquetReader::read_strings(int row_group, int col, vector<string> &out)
{
	const ParquetColumn &c = m_columns.at(col);

	if (c.physical_type != Type::BYTE_ARRAY || !c.other_type.empty())
		wrong_type(m_filename, row_group, c, "a string (BYTE_ARRAY)");
	read_as<string>(row_group, col, out, (const string *)NULL);
}
