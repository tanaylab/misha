#ifndef GENOMEARRAYSCSV_H_
#define GENOMEARRAYSCSV_H_

#include <cstdint>
#include <string>
#include <vector>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GIntervals.h"
#include "TGLException.h"
#include "UnknownChroms.h"

using namespace std;

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeArraysCsv {
public:
	enum Errors { FILE_ERROR, FORMAT_ERROR };

	// With clip_to_chroms, get_intervals() clips a record that runs past the end of its chromosome
	// to it, or drops it if it starts at or past the end; otherwise such a record is an error
	void init(const char *filename, const GenomeChromKey &chromkey, bool clip_to_chroms = false);

	const vector<string> &get_colnames() const { return m_colnames; }
	const GIntervals &get_intervals(int chromid);
	void  get_sliced_vals(GIntervals::const_iterator iinterval, vector<float> &vals);

	// Chromosome names that appear in the file but do not exist in the genome database.
	const UnknownChroms &get_unknown_chroms() const { return m_unknown_chroms; }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

	// The records get_intervals() clipped and dropped, and the first of them, "line N: (chrom,
	// start, end), chromosome size N" ("" if none)
	uint64_t get_num_clipped() const { return m_num_clipped; }
	uint64_t get_num_dropped() const { return m_num_dropped; }
	const string &get_first_clipped() const { return m_first_clipped; }

protected:
	struct Position {
		long bytes;
		long lineno;

		Position() {}
		Position(long _bytes, long _lineno) : bytes(_bytes), lineno(_lineno) {}
	};

	typedef vector<Position> Positions;
	typedef vector<Positions> ChromsPositions;

	BufferedFile          m_bfile;
	const GenomeChromKey *m_chromkey;
	ChromsPositions       m_chroms_positions;
	GIntervals            m_intervals;
	vector<string>        m_colnames;
	vector<string>        m_fields;
	UnknownChroms         m_unknown_chroms;
	bool                  m_clip_to_chroms{false};
	uint64_t              m_num_clipped{0};
	uint64_t              m_num_dropped{0};
	string                m_first_clipped;

	int read_fields(const Position &pos);
};

#endif

