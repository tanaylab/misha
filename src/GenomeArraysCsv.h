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

	void init(const char *filename, const GenomeChromKey &chromkey);

	const vector<string> &get_colnames() const { return m_colnames; }
	const GIntervals &get_intervals(int chromid);
	void  get_sliced_vals(GIntervals::const_iterator iinterval, vector<float> &vals);

	// Chromosome names that appear in the file but do not exist in the genome database.
	const UnknownChroms &get_unknown_chroms() const { return m_unknown_chroms; }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

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

	int read_fields(const Position &pos);
};

#endif

