/*
 * BinFinder.cpp
 *
 *  Created on: Dec 29, 2010
 *      Author: hoichman
 */

#include "BinFinder.h"

void BinFinder::init(const double *breaks, unsigned num_breaks, bool include_lowest, bool right)
{
	if (num_breaks < 2)
		TGLError<BinFinder>(BAD_NUM_BREAKS, "Invalid number of breaks %d", num_breaks);

	m_binsize = breaks[1] - breaks[0];
	m_include_lowest = include_lowest;
    m_right = right;

	m_breaks.clear();
	m_breaks.reserve(num_breaks);
	m_breaks.push_back(breaks[0]);

	for (unsigned i = 1; i < num_breaks; ++i) {
		if (breaks[i] == breaks[i - 1])
			TGLError<BinFinder>(NOT_UNIQUE_BREAKS, "Breaks are not unique (break[%d]=break[%d]=%g)", i - 1, i, breaks[i]);

		if (breaks[i] < breaks[i - 1])
			TGLError<BinFinder>(UNSORTED_BREAKS, "Breaks are not sorted (break[%d]=%g, break[%d]=%g)", i - 1, breaks[i - 1], i, breaks[i]);

		// we cast double to float to compensate for possible loss of precision
		if ((float)(breaks[i] - breaks[i - 1]) != (float)m_binsize)
			m_binsize = 0;

		m_breaks.push_back(breaks[i]);
	}

	// When breaks contain +/-Inf the per-step diff is Inf and the loop above
	// keeps m_binsize == Inf. The uniform-binsize fast path in val2bin() then
	// computes Inf/Inf = NaN, whose cast to int is undefined behaviour and
	// silently misroutes every value to a single bin. Force binary search.
	if (!std::isfinite(m_binsize))
		m_binsize = 0;
}
