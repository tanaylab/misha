/*
 * Parquet.h
 *
 * Parquet files without a Parquet library: the file metadata structs, ParquetWriter and
 * ParquetReader.
 */

#ifndef PARQUET_H_
#define PARQUET_H_

#include <cstdint>
#include <string>
#include <vector>

using namespace std;

// The file metadata of parquet.thrift (apache/parquet-format, src/main/thrift/parquet.thrift).
// Enum values and the field ids in Parquet.cpp are copied from that file. Each struct carries
// the fields misha writes and the ones a reader needs; encode() appends the struct in the Thrift
// compact protocol, decode() reads it back (resetting absent fields to their defaults and
// skipping fields it does not know). decode() followed by encode() reproduces what encode()
// wrote.
namespace parquet {

struct ThriftIn;  // a Thrift compact protocol reader, in Parquet.cpp

struct Type { enum type { BOOLEAN = 0, INT32 = 1, INT64 = 2, INT96 = 3, FLOAT = 4, DOUBLE = 5, BYTE_ARRAY = 6, FIXED_LEN_BYTE_ARRAY = 7 }; };
struct ConvertedType { enum type { UTF8 = 0, ENUM = 4, UINT_8 = 11, UINT_64 = 14, INT_64 = 18, JSON = 19 }; };
// The members of the LogicalType union, by field id
struct LogicalType { enum type { STRING = 1, ENUM = 4, INTEGER = 10, JSON = 12 }; };
struct FieldRepetitionType { enum type { REQUIRED = 0, OPTIONAL = 1, REPEATED = 2 }; };
struct Encoding { enum type { PLAIN = 0, PLAIN_DICTIONARY = 2, RLE = 3, BIT_PACKED = 4, DELTA_BINARY_PACKED = 5, DELTA_LENGTH_BYTE_ARRAY = 6, DELTA_BYTE_ARRAY = 7, RLE_DICTIONARY = 8, BYTE_STREAM_SPLIT = 9 }; };
struct CompressionCodec { enum type { UNCOMPRESSED = 0, SNAPPY = 1, GZIP = 2, LZO = 3, BROTLI = 4, LZ4 = 5, ZSTD = 6, LZ4_RAW = 7 }; };
struct PageType { enum type { DATA_PAGE = 0, INDEX_PAGE = 1, DICTIONARY_PAGE = 2, DATA_PAGE_V2 = 3 }; };

// min_value / max_value hold the PLAIN encoding of the value (no length prefix for BYTE_ARRAY)
struct Statistics {
	int64_t null_count = 0;
	string  max_value;
	string  min_value;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct SchemaElement {
	int     type = -1;                  // Type; -1: absent (the root)
	int     repetition_type = -1;       // FieldRepetitionType; -1: absent
	string  name;
	int     num_children = -1;          // -1: absent (a leaf)
	int     converted_type = -1;        // ConvertedType; -1: absent
	int     logical_type = -1;          // LogicalType (the union member's field id); -1: absent. encode() writes only STRING
	bool    integer_signed = true;      // LogicalType INTEGER: isSigned

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct DataPageHeader {
	int32_t num_values = 0;
	int     encoding = Encoding::PLAIN;
	int     definition_level_encoding = Encoding::RLE;
	int     repetition_level_encoding = Encoding::RLE;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct DictionaryPageHeader {
	int32_t num_values = 0;
	int     encoding = Encoding::PLAIN;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

// Only read: misha writes data page v1
struct DataPageHeaderV2 {
	int32_t num_values = 0;
	int32_t num_nulls = 0;
	int32_t num_rows = 0;
	int     encoding = Encoding::PLAIN;
	int32_t definition_levels_byte_length = 0;
	int32_t repetition_levels_byte_length = 0;
	bool    is_compressed = true;

	void decode(ThriftIn &in);
};

struct PageHeader {
	int                  type = PageType::DATA_PAGE;
	int32_t              uncompressed_page_size = 0;
	int32_t              compressed_page_size = 0;
	DataPageHeader       data_page_header;        // written for a DATA_PAGE
	DictionaryPageHeader dictionary_page_header;  // written for a DICTIONARY_PAGE
	DataPageHeaderV2     data_page_header_v2;     // read only

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct ColumnMetaData {
	int            type = Type::FLOAT;
	vector<int>    encodings;
	vector<string> path_in_schema;
	int            codec = CompressionCodec::UNCOMPRESSED;
	int64_t        num_values = 0;
	int64_t        total_uncompressed_size = 0;
	int64_t        total_compressed_size = 0;
	int64_t        data_page_offset = 0;
	int64_t        dictionary_page_offset = -1;  // -1: absent
	bool           has_statistics = false;
	Statistics     statistics;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct ColumnChunk {
	string         file_path;  // read only: the chunk is in another file
	int64_t        file_offset = 0;
	ColumnMetaData meta_data;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

struct RowGroup {
	vector<ColumnChunk> columns;
	int64_t             total_byte_size = 0;
	int64_t             num_rows = 0;
	int64_t             file_offset = 0;
	int64_t             total_compressed_size = 0;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

// encode() also writes column_orders: TYPE_ORDER for every leaf column, which is what tells
// readers that min_value / max_value follow the type's own order.
struct FileMetaData {
	int32_t               version = 2;
	vector<SchemaElement> schema;  // the root, then one leaf per column
	int64_t               num_rows = 0;
	vector<RowGroup>      row_groups;
	string                created_by;

	void encode(string &out) const;
	void decode(ThriftIn &in);
};

}

// Writes a Parquet file of intervals and values: chrom (string), start, end (INT32, or INT64
// with int64_coords) and one FLOAT column per value name, all OPTIONAL (nullable, though no
// value is null) and GZIP-compressed, as arrow writes them. Rows must arrive grouped by
// chromosome; a row group holds one chromosome and at most MAX_ROW_GROUP_ROWS rows.
// The file is complete only after close(); the destructor closes the file without finishing it.
//
// merge() concatenates files written by ParquetWriter with the same schema into one: their row
// groups in the order of parts, under one footer. A part without rows is left out unless no
// part has rows.
class ParquetWriter {
public:
	enum Errors { FILE_ERROR, FORMAT_ERROR };

	static const uint64_t MAX_ROW_GROUP_ROWS = 1048576;  // defined in Parquet.cpp too, for ODR uses

	ParquetWriter() = default;
	ParquetWriter(const ParquetWriter &) = delete;  // owns the file descriptor
	ParquetWriter &operator=(const ParquetWriter &) = delete;
	~ParquetWriter();

	// created_by: the file's created_by, "misha <version>"; name: what messages call the file
	// (the target that filename, a temporary file, stands for)
	void open(const string &filename, const vector<string> &value_names, bool int64_coords, const string &created_by, const string &name);
	void add_row(const string &chrom, int64_t start, int64_t end, const double *values);
	void close();

	// name: what messages call target
	static void merge(const string &target, const vector<string> &parts, const string &name);

private:
	int                      m_fd = -1;
	string                   m_filename;  // what messages call the file
	bool                     m_int64_coords = false;
	uint64_t                 m_offset = 0;
	parquet::FileMetaData    m_meta;

	// the row group being buffered
	string                   m_chrom;
	vector<int64_t>          m_starts;
	vector<int64_t>          m_ends;
	vector<vector<uint32_t>> m_values;  // float bit patterns, one vector per value column

	void create(const string &filename, const string &name);
	void flush_row_group();
	void write_coord_column(parquet::ColumnChunk &cc, const string &name, const vector<int64_t> &coords);
	template <typename T> void write_fixed_column(parquet::ColumnChunk &cc, int type, const string &name, const vector<T> &vals);
	void write_dictionary_chunk(parquet::ColumnMetaData &md, const string &dict, uint32_t dict_size, const vector<uint32_t> &indices);
	uint64_t write_page(parquet::PageHeader &header, const string &payload);
	void write(const string &bytes);
	void append(int fd, const string &filename, uint64_t offset, uint64_t size);
};

struct ParquetColumn {
	string name;
	int    physical_type;  // parquet::Type; -1 if the schema has none
	bool   optional;
	bool   unsigned_int;   // UINT_* or INTEGER(isSigned = false): INT32 and INT64 values are unsigned
	string other_type;     // a converted or logical type that makes the values neither numbers nor
	                       // strings to misha (DECIMAL, DATE, TIMESTAMP, ...); "" if none
};

// Reads the columns of a Parquet file with a flat schema (a root group of leaf columns) with
// pread() only, so that forked processes can read through the file descriptor they inherit
// (there is no file offset to share). Columns are decoded one column chunk (row group x column)
// at a time, and only when asked for: an unsupported type, annotation (ParquetColumn::other_type),
// encoding or codec is reported when that column is read, in an error that names the file, the
// row group and the column.
// Supported: INT32, INT64, FLOAT, DOUBLE, BYTE_ARRAY; PLAIN, PLAIN_DICTIONARY, RLE_DICTIONARY,
// BYTE_STREAM_SPLIT (values), RLE (levels); data pages v1 and v2; UNCOMPRESSED, SNAPPY, GZIP.
class ParquetReader {
public:
	enum Errors { FILE_ERROR, FORMAT_ERROR, UNSUPPORTED };

	ParquetReader() = default;
	ParquetReader(const ParquetReader &) = delete;  // owns the file descriptor
	ParquetReader &operator=(const ParquetReader &) = delete;
	~ParquetReader();

	// name: what messages call the file, if not filename
	void open(const string &filename, const string &name = string());

	const vector<ParquetColumn> &columns() const { return m_columns; }
	int column_index(const string &name) const;  // -1 if absent
	int num_row_groups() const { return (int)m_meta.row_groups.size(); }
	int64_t num_rows(int row_group) const { return m_meta.row_groups[row_group].num_rows; }
	const parquet::FileMetaData &metadata() const { return m_meta; }

	// One value per row of the row group
	void read_doubles(int row_group, int col, vector<double> &out);   // INT32/INT64/FLOAT/DOUBLE, null -> NaN
	void read_int64s(int row_group, int col, vector<int64_t> &out);   // INT32/INT64, or FLOAT/DOUBLE of whole numbers; null is an error
	void read_strings(int row_group, int col, vector<string> &out);   // BYTE_ARRAY of at most 1024 bytes, null is an error

private:
	friend class ParquetWriter;  // merge() copies the row groups through m_fd

	int                   m_fd = -1;
	string                m_filename;  // what messages call the file
	uint64_t              m_size = 0;
	uint64_t              m_footer_offset = 0;  // where the footer starts: the row groups end there
	parquet::FileMetaData m_meta;
	vector<ParquetColumn> m_columns;

	template <typename T, typename O> void read_as(int row_group, int col, vector<O> &out, const O *null_value);
	template <typename T> void read_chunk(int row_group, int col, vector<T> &vals, vector<uint8_t> &defs);
	void pread_exact(void *buf, uint64_t size, uint64_t offset);
};

#endif /* PARQUET_H_ */
