/*
 * Parquet.h
 *
 * Parquet files without a Parquet library: the file metadata structs and ParquetWriter.
 */

#ifndef PARQUET_H_
#define PARQUET_H_

#include <cstdint>
#include <string>
#include <vector>

using namespace std;

// The file metadata of parquet.thrift (apache/parquet-format, src/main/thrift/parquet.thrift).
// Enum values and the field ids in Parquet.cpp are copied from that file. Each struct carries
// the fields misha writes; encode() appends the struct in the Thrift compact protocol.
namespace parquet {

struct Type { enum type { BOOLEAN = 0, INT32 = 1, INT64 = 2, INT96 = 3, FLOAT = 4, DOUBLE = 5, BYTE_ARRAY = 6, FIXED_LEN_BYTE_ARRAY = 7 }; };
struct ConvertedType { enum type { UTF8 = 0 }; };
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
};

struct SchemaElement {
	int     type = -1;                  // Type; -1: absent (the root)
	int     repetition_type = -1;       // FieldRepetitionType; -1: absent
	string  name;
	int     num_children = -1;          // -1: absent (a leaf)
	int     converted_type = -1;        // ConvertedType; -1: absent
	bool    logical_type_string = false; // LogicalType STRING

	void encode(string &out) const;
};

struct DataPageHeader {
	int32_t num_values = 0;
	int     encoding = Encoding::PLAIN;
	int     definition_level_encoding = Encoding::RLE;
	int     repetition_level_encoding = Encoding::RLE;

	void encode(string &out) const;
};

struct DictionaryPageHeader {
	int32_t num_values = 0;
	int     encoding = Encoding::PLAIN;

	void encode(string &out) const;
};

struct PageHeader {
	int                  type = PageType::DATA_PAGE;
	int32_t              uncompressed_page_size = 0;
	int32_t              compressed_page_size = 0;
	DataPageHeader       data_page_header;        // written for a DATA_PAGE
	DictionaryPageHeader dictionary_page_header;  // written for a DICTIONARY_PAGE

	void encode(string &out) const;
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
};

struct ColumnChunk {
	int64_t        file_offset = 0;
	ColumnMetaData meta_data;

	void encode(string &out) const;
};

struct RowGroup {
	vector<ColumnChunk> columns;
	int64_t             total_byte_size = 0;
	int64_t             num_rows = 0;
	int64_t             file_offset = 0;
	int64_t             total_compressed_size = 0;

	void encode(string &out) const;
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
};

}

// Writes a Parquet file of intervals and values: chrom (string), start, end (INT32, or INT64
// with int64_coords) and one FLOAT column per value name, all REQUIRED. Rows must arrive
// grouped by chromosome; a row group holds one chromosome and at most MAX_ROW_GROUP_ROWS rows.
// The file is complete only after close(); the destructor closes the file without finishing it.
class ParquetWriter {
public:
	enum Errors { FILE_ERROR, FORMAT_ERROR };

	static const uint64_t MAX_ROW_GROUP_ROWS = 1048576;

	ParquetWriter() = default;
	ParquetWriter(const ParquetWriter &) = delete;  // owns the file descriptor
	ParquetWriter &operator=(const ParquetWriter &) = delete;
	~ParquetWriter();

	void open(const string &filename, const vector<string> &value_names, bool int64_coords);
	void add_row(const string &chrom, int64_t start, int64_t end, const double *values);
	void close();

private:
	int                      m_fd = -1;
	string                   m_filename;
	bool                     m_int64_coords = false;
	uint64_t                 m_offset = 0;
	parquet::FileMetaData    m_meta;

	// the row group being buffered
	string                   m_chrom;
	vector<int64_t>          m_starts;
	vector<int64_t>          m_ends;
	vector<vector<uint32_t>> m_values;  // float bit patterns, one vector per value column

	void flush_row_group();
	void write_coord_column(parquet::ColumnChunk &cc, const string &name, const vector<int64_t> &coords);
	template <typename T> void write_fixed_column(parquet::ColumnChunk &cc, int type, const string &name, const vector<T> &vals);
	void write_dictionary_chunk(parquet::ColumnMetaData &md, const string &dict, uint32_t dict_size, const vector<uint32_t> &indices);
	void write_page(parquet::PageHeader &header, const string &payload);
	void write(const string &bytes);
};

#endif /* PARQUET_H_ */
