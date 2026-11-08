#ifndef WIG_H_
#define WIG_H_

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GIntervals.h"
#include "UnknownChroms.h"

using namespace std;

//----------------------------------------------------------------------------
// Wig class parses file in WIG or BedGraph formats and returns its data
//
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class Wig {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, BAD_CHROM, BAD_START, BAD_STEP, BAD_SPAN, BAD_VAL, BAD_COORD };

	Wig() {}
	Wig(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms) { init(chromkey, filename, ignore_unknown_chroms); }

	void init(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms);

	// returns true if data for chromosome exists; (float)udata of each interval contains the value
	bool get_data(int chromid, GIntervals &intervals);

	// Chromosome names that appear in the file but do not exist in the genome database.
	const UnknownChroms &get_unknown_chroms() const { return m_unknown_chroms; }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

	// Records read by get_data() that ran past the end of their chromosome: clipped to it, or
	// dropped when they start at or past it; and the first of them, "line N: (chrom, start, end),
	// chromosome size N" ("" if none)
	uint64_t get_num_clipped() const { return m_num_clipped; }
	uint64_t get_num_dropped() const { return m_num_dropped; }
	const string &get_first_clipped() const { return m_first_clipped; }

private:
	enum { CHROM_FIELD, START_FIELD, STEP_FIELD, SPAN_FIELD, NUM_FIELDS };
	enum RecordType { FIXED_STEP_REC, VAR_STEP_REC, VAL_REC, COORD_VAL_REC, BEDGRAPH_REC };

	struct FixedStepRec {
		int     chromid;
		int64_t start;
		int64_t step;
		int64_t span;
	};

	struct VarStepRec {
		int     chromid;
		int64_t span;
	};

	struct CoordValRec {
		int64_t coord;
		float   val;
	};

	struct BedGraph {
		int     chromid;
		int64_t start;
		int64_t end;
		float   val;
	};

	struct Rec {
		RecordType type;
		union {
			FixedStepRec fixed_step;
			VarStepRec   var_step;
			CoordValRec  coord_val;
			float        val;
			BedGraph     bedgraph;
		};
	};

	// A field of the line being parsed: NUL-terminated in the line buffer. len is its full
	// length, which a NUL byte inside the field makes longer than strlen(str).
	struct Field {
		char  *str;
		size_t len;

		bool operator==(const char *s) const { return len == strlen(s) && !memcmp(str, s, len); }
		bool starts_with(const string &prefix) const { return len >= prefix.size() && !memcmp(str, prefix.data(), prefix.size()); }
	};

	static string FIELDS_STRS[NUM_FIELDS];

	GenomeChromKey *m_chromkey;
	vector<float>   m_chrom_data;
	vector<long>    m_chrom_fpos;
	vector<int64_t> m_chrom_lineno;
	BufferedFile    m_bfile;
	bool            m_ignore_unknown_chroms;
	UnknownChroms   m_unknown_chroms;

	// The file is read a block at a time into m_buf and split into fields in place.
	// [m_buf_pos, m_buf_end) is the part not parsed yet; m_buf[0] is at offset m_buf_fpos.
	vector<char>    m_buf;
	size_t          m_buf_pos{0};
	size_t          m_buf_end{0};
	int64_t         m_buf_fpos{0};
	bool            m_file_end{false};
	vector<Field>   m_fields;

	uint64_t        m_num_clipped{0};
	uint64_t        m_num_dropped{0};
	string          m_first_clipped;

	// Notes a record that runs past the end of the chromosome; returns true if it is to be dropped
	bool clip(int chromid, int64_t start, int64_t end, int64_t lineno);

	// The last chromosome name looked up and its id (-1: not in the database)
	string          m_last_chrom;
	int             m_last_chromid{-1};
	bool            m_has_last_chrom{false};

	int64_t tell() const { return m_buf_fpos + m_buf_pos; }
	void    seek(int64_t pos);
	bool    next_line(char *&line, size_t &len, bool &terminated);
	int     split_fields();

	bool read_record(Rec &rec, int64_t &lineno);

	int     str2chromid(const char *str, int64_t lineno);
	int64_t str2span(const char *str, int64_t lineno);
	int64_t str2step(const char *str, int64_t lineno);
	int64_t str2start(const char *str, int64_t lineno);
	int64_t str2coord(const char *str, int64_t base, int64_t lineno);
	float   str2val(const char *str, int64_t lineno);
};


//-------------------------------------IMPLEMENTATION -----------------------------------------

// Consecutive records nearly always name the same chromosome: the last name is looked up once
inline int Wig::str2chromid(const char *str, int64_t lineno)
{
	size_t len = strlen(str);

	if (m_has_last_chrom && len == m_last_chrom.size() && !memcmp(str, m_last_chrom.data(), len))
		return m_last_chromid;

	int chromid;
	try {
		chromid = m_chromkey->chrom2id(str);
	} catch (TGLException &e) {
		if (!m_ignore_unknown_chroms)
			TGLError<Wig>(BAD_CHROM, "WIG file %s, line %ld: %s\n", m_bfile.file_name().c_str(), lineno, e.msg());
		m_unknown_chroms.record(str);
		chromid = -1;
	}
	m_last_chrom.assign(str, len);
	m_last_chromid = chromid;
	m_has_last_chrom = true;
	return chromid;
}

inline int64_t Wig::str2span(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t span = strtoll(str, &endptr, 10);
	if (*endptr || span < 1)
		TGLError<Wig>(BAD_SPAN, "WIG file %s, line %ld: invalid value of span", m_bfile.file_name().c_str(), lineno);
	return span;
}

inline int64_t Wig::str2step(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t step = strtoll(str, &endptr, 10);
	if (*endptr || step < 1)
		TGLError<Wig>(BAD_STEP, "WIG file %s, line %ld: invalid value of step", m_bfile.file_name().c_str(), lineno);
	return step;
}

inline int64_t Wig::str2start(const char *str, int64_t lineno)
{
	char *endptr;
	int64_t start = strtoll(str, &endptr, 10);
	if (*endptr || start < 1)
		TGLError<Wig>(BAD_START, "WIG file %s, line %ld: invalid value of start", m_bfile.file_name().c_str(), lineno);
	return start - 1;
}

inline int64_t Wig::str2coord(const char *str, int64_t base, int64_t lineno)
{
	char *endptr;
	int64_t coord = strtoll(str, &endptr, 10);
	if (*endptr || coord < base)
		TGLError<Wig>(BAD_COORD, "WIG file %s, line %ld: invalid coordinate", m_bfile.file_name().c_str(), lineno);
	return coord - base;
}

inline float Wig::str2val(const char *str, int64_t lineno)
{
	char *endptr;
	float val = strtod(str, &endptr);
	if (*endptr)
		TGLError<Wig>(BAD_VAL, "WIG file %s, line %ld: invalid value", m_bfile.file_name().c_str(), lineno);
	return val;
}

#endif /* WIG_H_ */
