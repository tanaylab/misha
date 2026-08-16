#ifndef GENOMEARRAYSCSV_H_
#define GENOMEARRAYSCSV_H_

#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GIntervals.h"
#include "TGLException.h"

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
	// Aliases are resolved by GenomeChromKey::chrom2id, so a name reachable through the
	// alias chain never lands here. Only the first MAX_REPORTED_UNKNOWN_CHROMS distinct
	// names are kept; get_num_unknown_chroms() returns the true number of distinct names.
	const vector<string> &get_unknown_chroms() const { return m_unknown_chroms; }
	uint64_t get_num_unknown_chroms() const { return m_unknown_chroms_seen.size(); }

	// Number of database chromosomes for which the file contains data.
	uint64_t get_num_matched_chroms() const;

protected:
	static const size_t MAX_REPORTED_UNKNOWN_CHROMS = 32;

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
	vector<string>        m_unknown_chroms;
	set<string>           m_unknown_chroms_seen;

	int read_fields(const Position &pos);

	void record_unknown_chrom(const string &chrom);
};

#endif

