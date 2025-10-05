#include <errno.h>

#include "port.h"

#include "GenomeSeqFetch.h"

void GenomeSeqFetch::read_interval(const GInterval &interval, const GenomeChromKey &chromkey, vector<char> &result)
{
	// Check cache first
	if (m_cache_valid &&
	    m_cached_interval.chromid == interval.chromid &&
	    m_cached_interval.start == interval.start &&
	    m_cached_interval.end == interval.end &&
	    m_cached_interval.strand == interval.strand) {
		// Cache hit! Return cached sequence
		result = m_cached_seq;
		return;
	}

	if (m_cur_chromid != interval.chromid) {
		char filename[PATH_MAX];

		m_cur_chromid = interval.chromid;
		m_cache_valid = false;  // Invalidate cache on chromosome change
		snprintf(filename, sizeof(filename), "%s/%s.seq", m_seqdir.c_str(), chromkey.id2chrom(interval.chromid).c_str());
		m_bfile.close();
		m_bfile.open(filename, "rb");

		if (m_bfile.error())
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", filename, strerror(errno));
	}

	interval.verify(chromkey, false);
	result.clear();

	int64_t size = min(interval.end, m_bfile.file_size()) - interval.start;

	if (size < 0)
		return;

	if (!size)
		size++;

	result.resize(size);
	m_bfile.seek(interval.start, SEEK_SET);
	if (m_bfile.read(&*result.begin(), result.size()) != result.size()) {
		if (m_bfile.error())
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", m_bfile.file_name().c_str(), strerror(errno));

		TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed", m_bfile.file_name().c_str());
	}

	// if strand == -1 make the sequence reverse-complementary
	if (interval.strand == -1) {
		for (vector<char>::iterator i = result.begin(); i != result.end(); ++i)
			*i = basepair2complementary(*i);
		reverse(result.begin(), result.end());
	}

	// Update cache with the newly read sequence
	m_cached_interval = interval;
	m_cached_seq = result;
	m_cache_valid = true;
}
