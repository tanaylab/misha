#include "IntervalValidator.h"
#include "GenomeChromKey.h"
#include "rdbutils.h" // For verror
#include <cmath> // For ceil
#include <algorithm> // For max

using namespace std;
using namespace rdb;

IntervalValidator::IntervalValidator(const GenomeChromKey &chromkey) :
	m_chromkey(chromkey)
{
}

void IntervalValidator::restrict_bins(int64_t maxbins, GIntervals &intervals, unsigned binsize) const
{
	for (GIntervals::const_iterator iinterval = intervals.begin(); iinterval != intervals.end(); ++iinterval) {
		int64_t bins = max((int64_t)0, (int64_t)ceil(iinterval->end / binsize) - (int64_t)(iinterval->start / binsize));

		if (bins > maxbins)
			verror("The interval %s [%ld, %ld) covers too wide range of samples that might cause memory allocation failure.\n"
					"(bins covered: %ld, bins limit: %ld)\n", m_chromkey.id2chrom(iinterval->chromid).c_str(), iinterval->start, iinterval->end, bins, maxbins);
	}
}

