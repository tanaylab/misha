/*
 * Parquet.cpp
 *
 * Thrift compact protocol: https://github.com/apache/thrift/blob/master/doc/specs/thrift-compact-protocol.md
 * Encodings (PLAIN, RLE / bit-packed hybrid, RLE_DICTIONARY): apache/parquet-format, Encodings.md
 */

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <unordered_map>

#include "Parquet.h"
#include "TGLException.h"
#include "rdbutils.h"

using namespace parquet;

// PLAIN is little-endian; values are copied to the file as they are in memory
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "Parquet output assumes a little-endian host");

//------------------------------ Thrift compact protocol ------------------------------

namespace {

enum { CT_I32 = 5, CT_I64 = 6, CT_BINARY = 8, CT_LIST = 9, CT_STRUCT = 12 };

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
	if (logical_type_string) {
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

//------------------------------ ParquetWriter ------------------------------

// Column chunks with more distinct values (by bit pattern) are written PLAIN: arrow's 1 MB
// dictionary page limit, in floats
static const uint32_t MAX_DICTIONARY_SIZE = 262144;

ParquetWriter::~ParquetWriter()
{
	if (m_fd >= 0)
		::close(m_fd);
}

void ParquetWriter::open(const string &filename, const vector<string> &value_names, bool int64_coords)
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

	SEXP version = rdb::run_in_R("as.character(utils::packageVersion(\"misha\"))", R_BaseEnv);
	if (!Rf_isString(version) || Rf_length(version) != 1)
		TGLError<ParquetWriter>(FORMAT_ERROR, "Failed to read the version of misha");

	m_meta = FileMetaData();
	m_meta.created_by = string("misha ") + CHAR(STRING_ELT(version, 0));
	m_meta.schema.resize(4 + value_names.size());
	m_meta.schema[0].name = "schema";
	m_meta.schema[0].num_children = 3 + value_names.size();
	for (size_t i = 1; i < m_meta.schema.size(); ++i) {
		SchemaElement &col = m_meta.schema[i];
		col.repetition_type = FieldRepetitionType::REQUIRED;
		col.type = i == 1 ? Type::BYTE_ARRAY : i <= 3 ? (int64_coords ? Type::INT64 : Type::INT32) : Type::FLOAT;
		col.name = i <= 3 ? coords[i - 1] : value_names[i - 4];
	}
	m_meta.schema[1].converted_type = ConvertedType::UTF8;
	m_meta.schema[1].logical_type_string = true;

	m_int64_coords = int64_coords;
	m_chrom.clear();
	m_starts.clear();
	m_ends.clear();
	m_values.assign(value_names.size(), vector<uint32_t>());

	m_filename = filename;
	m_fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
	if (m_fd < 0)
		TGLError<ParquetWriter>(FILE_ERROR, "Failed to open file %s for writing: %s", filename.c_str(), strerror(errno));
	m_offset = 0;
	write("PAR1");
}

void ParquetWriter::add_row(const string &chrom, int64_t start, int64_t end, const double *values)
{
	if (!m_int64_coords && (start < INT32_MIN || end > INT32_MAX))
		TGLError<ParquetWriter>(FORMAT_ERROR, "Interval (%s, %lld, %lld) does not fit 32-bit Parquet coordinates", chrom.c_str(), (long long)start, (long long)end);

	if (chrom != m_chrom || m_starts.size() >= MAX_ROW_GROUP_ROWS) {
		flush_row_group();
		m_chrom = chrom;
	}

	m_starts.push_back(start);
	m_ends.push_back(end);
	for (size_t i = 0; i < m_values.size(); ++i) {
		float v = (float)values[i];
		uint32_t bits;
		memcpy(&bits, &v, sizeof(bits));
		m_values[i].push_back(bits);
	}
}

void ParquetWriter::close()
{
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
	if (m_starts.empty())
		return;

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
	chrom.meta_data.has_statistics = true;
	chrom.meta_data.statistics.min_value = chrom.meta_data.statistics.max_value = m_chrom;

	write_coord_column(rg.columns[1], "start", m_starts);
	write_coord_column(rg.columns[2], "end", m_ends);
	for (size_t i = 0; i < m_values.size(); ++i)
		write_fixed_column(rg.columns[i + 3], Type::FLOAT, m_meta.schema[i + 4].name, m_values[i]);

	for (const ColumnChunk &cc : rg.columns) {
		rg.total_byte_size += cc.meta_data.total_uncompressed_size;
		rg.total_compressed_size += cc.meta_data.total_compressed_size;
	}
	m_meta.num_rows += rg.num_rows;
	m_meta.row_groups.push_back(rg);

	m_starts.clear();
	m_ends.clear();
	for (vector<uint32_t> &v : m_values)
		v.clear();
}

void ParquetWriter::write_coord_column(ColumnChunk &cc, const string &name, const vector<int64_t> &coords)
{
	auto minmax = minmax_element(coords.begin(), coords.end());
	Statistics &stats = cc.meta_data.statistics;

	if (m_int64_coords) {
		write_fixed_column(cc, Type::INT64, name, coords);
		stats.min_value = plain(*minmax.first);
		stats.max_value = plain(*minmax.second);
	} else {
		write_fixed_column(cc, Type::INT32, name, vector<int32_t>(coords.begin(), coords.end()));
		stats.min_value = plain((int32_t)*minmax.first);
		stats.max_value = plain((int32_t)*minmax.second);
	}
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
	write_page(header, string((const char *)vals.data(), vals.size() * sizeof(T)));
	md.encodings = { Encoding::PLAIN, Encoding::RLE };
	md.total_uncompressed_size = md.total_compressed_size = m_offset - cc.file_offset;
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
	write_page(header, dict);

	int bit_width = 1;
	while (((uint64_t)1 << bit_width) < dict_size)
		++bit_width;
	string data(1, (char)bit_width);
	put_rle_hybrid(data, indices, bit_width);

	header = PageHeader();
	header.data_page_header.num_values = indices.size();
	header.data_page_header.encoding = Encoding::RLE_DICTIONARY;
	md.data_page_offset = m_offset;
	write_page(header, data);

	md.encodings = { Encoding::PLAIN, Encoding::RLE, Encoding::RLE_DICTIONARY };
	md.total_uncompressed_size = md.total_compressed_size = m_offset - start;
}

void ParquetWriter::write_page(PageHeader &header, const string &payload)
{
	header.uncompressed_page_size = header.compressed_page_size = payload.size();
	string bytes;
	header.encode(bytes);
	write(bytes);
	write(payload);
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
