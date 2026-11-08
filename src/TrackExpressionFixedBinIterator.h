#ifndef TRACKEXPRESSIONFIXEDBINITERATOR_H_
#define TRACKEXPRESSIONFIXEDBINITERATOR_H_

#include "TrackExpressionIterator.h"

//------------------------------ TrackExpressionFixedBinIterator ---------------------------------

class TrackExpressionFixedBinIterator : public TrackExpression1DIterator {
public:
	TrackExpressionFixedBinIterator() : TrackExpression1DIterator(FIXED_BIN), m_binsize(0) {}

	bool begin(int64_t binsize, GIntervalsFetcher1D &scope);
	virtual bool next();

	int64_t get_bin_size() const { return m_binsize; }

	// The bin of the current interval and the end bin of the current scope interval: next() gives
	// the bins up to get_end_bin() - 1 and then goes on to the next scope interval.
	int64_t get_cur_bin() const { return m_cur_bin; }
	int64_t get_end_bin() const { return m_end_bin; }

	// Moves k intervals forward, as k calls of next() would; 1 <= k <= get_end_bin() - get_cur_bin().
	bool advance(int64_t k);

private:
	int64_t     m_binsize;
	int64_t     m_cur_bin;
	int64_t     m_end_bin;
};

#endif /* TRACKEXPRESSIONFIXEDBINITERATOR_H_ */
