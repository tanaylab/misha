#ifndef TRACKEXPRESSIONROWS1DITERATOR_H_
#define TRACKEXPRESSIONROWS1DITERATOR_H_

#include "TrackExpressionIterator.h"

// Walks a given list of 1D intervals in order, repeats included, with no
// sorting, unification or scope intersection (gextract.variants evaluates
// the same bin once per substitution). The row index is reported as the
// scope index, which is how the scanner finds each row's substitution.
class TrackExpressionRows1DIterator : public TrackExpression1DIterator {
public:
	TrackExpressionRows1DIterator() : TrackExpression1DIterator(INTERVALS1D) {}

	bool begin(const GIntervals &rows, GIntervalsFetcher1D &scope) {
		TrackExpression1DIterator::begin(scope);
		m_rows = &rows;
		m_idx = 0;
		return set_cur();
	}

	virtual bool next() { ++m_idx; return set_cur(); }

	virtual int64_t get_cur_scope_idx() const { return (int64_t)m_idx; }
	virtual int64_t get_cur_scope_chrom_idx() const { return -1; }

private:
	const GIntervals *m_rows = nullptr;
	size_t            m_idx = 0;

	bool set_cur() {
		if (m_idx >= m_rows->size()) {
			end();
			return false;
		}
		m_last_interval = m_last_scope_interval = (*m_rows)[m_idx];
		return true;
	}
};

#endif
