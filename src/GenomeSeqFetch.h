/*
 * GenomeSeqFetch.h
 *
 *  Created on: Nov 21, 2010
 *      Author: hoichman
 */

#ifndef GENOMESEQFETCH_H_
#define GENOMESEQFETCH_H_

#include <string>
#include <vector>

#include "BufferedFile.h"
#include "GenomeChromKey.h"
#include "GenomeUtils.h"
#include "GInterval.h"

// Forward declaration
class GenomeIndex;

// -------------------- GenomeSeqFetch  -----------------------
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeSeqFetch {
public:
	enum Errors { FILE_READ_FAILED, INVALID_INTERVAL };

	GenomeSeqFetch();
	~GenomeSeqFetch();

	void set_seqdir(const std::string &dir);
	void read_interval(const GInterval &interval, const GenomeChromKey &chromkey, std::vector<char> &result);

private:
	std::string  m_seqdir;
	int          m_cur_chromid;
	BufferedFile m_bfile;

	// Single-interval cache for sequence-based vtracks
	bool               m_cache_valid;
	GInterval          m_cached_interval;
	std::vector<char>  m_cached_seq;

	// Indexed genome support
	bool               m_indexed_mode;
	GenomeIndex*       m_index;
};

#endif /* GENOMESEQFETCH_H_ */
