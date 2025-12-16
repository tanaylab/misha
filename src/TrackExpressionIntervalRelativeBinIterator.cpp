#include <math.h>
#include <algorithm>

#include "rdbutils.h"
#include "TrackExpressionIntervalRelativeBinIterator.h"

using namespace std;
using namespace rdb;

bool TrackExpressionIntervalRelativeBinIterator::begin(
	const GIntervals &src_intervals,
	int64_t binsize,
	PartialBinsMode partial_bins_mode,
	GIntervalsFetcher1D &scope)
{
	TrackExpression1DIterator::begin(scope);

	if (binsize <= 0)
		verror("Bin size must be positive (got %ld)", binsize);

	m_src_intervals = (GIntervals *)&src_intervals;
	m_binsize = binsize;
	m_partial_bins_mode = partial_bins_mode;
	m_icur_src_interval = m_src_intervals->begin() - 1;
	m_cur_bin = 0;
	m_num_bins = 0;
	m_cur_src_interval_idx = -1;

	if (m_src_intervals->empty()) {
		end();
		return false;
	}

	return next();
}

void TrackExpressionIntervalRelativeBinIterator::compute_bins_for_current_interval()
{
	int64_t interval_size = m_icur_src_interval->end - m_icur_src_interval->start;

	if (m_partial_bins_mode == EXACT) {
		// Only full bins - truncate division
		m_num_bins = interval_size / m_binsize;
	} else {
		// Include partial last bin (will be clipped) - round up
		m_num_bins = (int64_t)ceil((double)interval_size / m_binsize);
	}

	m_cur_bin = 0;
}

bool TrackExpressionIntervalRelativeBinIterator::advance_to_next_src_interval()
{
	++m_icur_src_interval;
	++m_cur_src_interval_idx;

	while (m_icur_src_interval != m_src_intervals->end()) {
		compute_bins_for_current_interval();

		if (m_num_bins > 0) {
			m_last_interval.chromid = m_icur_src_interval->chromid;
			m_last_scope_interval = *m_icur_src_interval;
			return true;
		}

		// Skip intervals too small for any bins (only happens in EXACT mode)
		++m_icur_src_interval;
		++m_cur_src_interval_idx;
	}

	return false;
}

bool TrackExpressionIntervalRelativeBinIterator::next()
{
	if (isend())
		return false;

	// Try to advance within current source interval
	if (m_cur_src_interval_idx >= 0 && m_cur_bin < m_num_bins - 1) {
		m_cur_bin++;
	} else {
		// Need to move to next source interval
		if (!advance_to_next_src_interval()) {
			end();
			return false;
		}
	}

	// Compute bin coordinates relative to source interval start
	int64_t bin_start = m_icur_src_interval->start + m_cur_bin * m_binsize;
	int64_t bin_end = bin_start + m_binsize;

	// Clip to source interval boundaries (for CLIP mode, last bin may be partial)
	m_last_interval.start = bin_start;
	m_last_interval.end = min(bin_end, m_icur_src_interval->end);
	m_last_interval.chromid = m_icur_src_interval->chromid;

	return true;
}
