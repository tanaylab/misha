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
//
// !!!!!!!!! NEVER READ THROUGH ONE INSTANCE FROM SEVERAL FORKED PROCESSES !!!!!!!!!!!!!!!!
// Forked processes inherit the same open file description, i.e. the same file offset, so two
// children seeking and reading concurrently hand each other's bytes back -- silently, with the
// right length and the wrong bases. Multitasking code must construct its GenomeSeqFetch (and
// anything owning one, such as TrackExprScanner) AFTER distribute_task() returns in the child.
// set_seqdir() opens genome.seq immediately on an indexed database, and the per-chromosome mode
// keeps its file open across reads, so both formats are affected.

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
