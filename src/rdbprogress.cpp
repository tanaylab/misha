/*
 * rdbprogress.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <sys/time.h>

#include "rdbutils.h"

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

#include "rdbprogress.h"

using namespace rdb;

void Progress_reporter::init(uint64_t maxsteps, uint64_t init_report_step, uint64_t report_interval, uint64_t min_report_interval)
{
	m_maxsteps = maxsteps;
	m_report_step = init_report_step;
	m_report_interval = report_interval;
	m_min_report_interval = min_report_interval;

	m_numsteps = 0;
	m_numsteps_from_last_report = 0;
	m_last_progress_reported = -1;
	m_last_report_clock = get_cur_clock();
	m_elapsed_clock = 0;
}

uint64_t Progress_reporter::get_cur_clock()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

void Progress_reporter::report(uint64_t delta_steps_done)
{
	monitor_memusage();

	m_numsteps_from_last_report += delta_steps_done;
	m_numsteps += delta_steps_done;
	if (m_numsteps_from_last_report > m_report_step) {
		uint64_t curclock = get_cur_clock();
		double delta = curclock - m_last_report_clock;

		if (delta)
			m_report_step = (int)(m_report_step * (m_report_interval / delta) + .5);
		else
			m_report_step *= 10;

		if (delta > m_min_report_interval) {
			int progress = m_maxsteps ? (int)(100. * m_numsteps / m_maxsteps) : 0;

			if (m_last_progress_reported < 0 && !m_report_prefix.empty())
				REprintf("%s", m_report_prefix.c_str());

			if (progress != m_last_progress_reported) {
				REprintf("%d%%...", progress);
				update_progress(progress);
			} else
				REprintf(".");
			m_last_progress_reported = progress;
			m_numsteps_from_last_report = 0;
			m_last_report_clock = curclock;
			m_elapsed_clock = (uint64_t)delta;
		}
	}
}

void Progress_reporter::report_last()
{
	if (m_last_progress_reported >= 0) {
		if (m_last_progress_reported != 100)
			REprintf("100%%\n");
		else
			REprintf("\n");
	}
	update_progress(100);
}
