#ifndef INTERVALVALIDATOR_H_
#define INTERVALVALIDATOR_H_

#include <cstdint>
#include "GIntervals.h"

class GenomeChromKey;

// Utility class for validating intervals
class IntervalValidator {
public:
	// Constructor takes a reference to chromosome key for chromosome name lookups
	explicit IntervalValidator(const GenomeChromKey &chromkey);

	// Verifies that the number of bins in each interval does not exceed the limit
	// Throws an error if any interval exceeds the bin limit
	void restrict_bins(int64_t maxbins, GIntervals &intervals, unsigned binsize) const;

private:
	const GenomeChromKey &m_chromkey;
};

#endif /* INTERVALVALIDATOR_H_ */

