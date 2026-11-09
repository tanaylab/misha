#include <cstdint>
#include <cmath>
#include <limits>
#include <cstdio>

#include "ConfigurationManager.h"
#include "ConfigurationDefaults.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

bool ConfigurationManager::get_multitasking() const
{
	// A call nested in another misha call (made by R code that the outer call runs, such as
	// the FUN of gintervals.mapply) does not multitask: the shared memory, semaphores and kid
	// slots belong to the outermost call, and a nested prepare4multitasking() would reuse them.
	if (RdbInitializer::get_ref_count() > 1)
		return false;

	if (m_multitasking < 0) {
		SEXP r_multitasking = Rf_GetOption1(Rf_install("gmultitasking"));

		if (Rf_isLogical(r_multitasking) && Rf_length(r_multitasking) >= 1) {
			if (LOGICAL(r_multitasking)[0] == NA_LOGICAL)
				verror("options(gmultitasking) must be TRUE or FALSE, not NA");
			m_multitasking = (int)LOGICAL(r_multitasking)[0];
		}
		else
			m_multitasking = false;
	}
	return (bool)m_multitasking;
}

uint64_t ConfigurationManager::get_max_processes() const
{
	if (!m_max_processes) {
		SEXP r_max_processes = Rf_GetOption1(Rf_install("gmax.processes"));

		if (Rf_isReal(r_max_processes))
			m_max_processes = (uint64_t)REAL(r_max_processes)[0];
		else if (Rf_isInteger(r_max_processes))
			m_max_processes = INTEGER(r_max_processes)[0];
		else
			m_max_processes = misha::config::DEFAULT_MAX_PROCESSES;
		if (m_max_processes < 1)
			m_max_processes = misha::config::DEFAULT_MAX_PROCESSES;
	}
	return m_max_processes;
}

uint64_t ConfigurationManager::get_max_processes2core() const
{
	if (!m_max_processes2core) {
		SEXP r_max_processes2core = Rf_GetOption1(Rf_install("gmax.processes2core"));

		if (Rf_isReal(r_max_processes2core))
			m_max_processes2core = (uint64_t)REAL(r_max_processes2core)[0];
		else if (Rf_isInteger(r_max_processes2core))
			m_max_processes2core = INTEGER(r_max_processes2core)[0];
		else
			m_max_processes2core = misha::config::DEFAULT_MAX_PROCESSES2CORE;
		if (m_max_processes2core < 1)
			m_max_processes2core = misha::config::DEFAULT_MAX_PROCESSES2CORE;
	}
	return m_max_processes2core;
}

uint64_t ConfigurationManager::get_min_scope4process() const
{
	if (!m_min_scope4process) {
		SEXP r_min_scope4process = Rf_GetOption1(Rf_install("gmin.scope4process"));

		if (Rf_isReal(r_min_scope4process))
			m_min_scope4process = (uint64_t)REAL(r_min_scope4process)[0];
		else if (Rf_isInteger(r_min_scope4process))
			m_min_scope4process = INTEGER(r_min_scope4process)[0];
		else
			m_min_scope4process = misha::config::DEFAULT_MIN_SCOPE4PROCESS;
	}
	return m_min_scope4process;
}

uint64_t ConfigurationManager::get_min_seqs_work4process() const
{
	if (!m_min_seqs_work4process) {
		SEXP r_min_seqs_work4process = Rf_GetOption1(Rf_install("gmin.seqs.work4process"));

		if (Rf_isReal(r_min_seqs_work4process))
			m_min_seqs_work4process = (uint64_t)REAL(r_min_seqs_work4process)[0];
		else if (Rf_isInteger(r_min_seqs_work4process))
			m_min_seqs_work4process = INTEGER(r_min_seqs_work4process)[0];
		else
			m_min_seqs_work4process = misha::config::DEFAULT_MIN_SEQS_WORK4PROCESS;
	}
	return m_min_seqs_work4process;
}

uint64_t ConfigurationManager::get_max_data_size() const
{
	if (!m_max_data_size) {
		SEXP r_max_data_size = Rf_GetOption1(Rf_install("gmax.data.size"));

		// options(gmax.data.size = NA) gives a logical NA, which is as invalid as NA_real_
		if (Rf_isLogical(r_max_data_size) && Rf_length(r_max_data_size) >= 1 && LOGICAL(r_max_data_size)[0] == NA_LOGICAL)
			verror("options(gmax.data.size) must be a number, not NA");

		if (is_num_vector(r_max_data_size) && Rf_length(r_max_data_size) >= 1) {
			double size = num_elt(r_max_data_size, 0);

			// Casting NA, Inf, a negative value or one of 2^64 or more to uint64_t is undefined.
			// Inf and any value of 2^53 or more mean no limit. A value in [0, 1) counts as 1,
			// as in gtrack.liftover: 0 could not be told apart from "not read yet" above, and a
			// quantile sampling buffer of size 0 crashed R.
			if (std::isnan(size))
				verror("options(gmax.data.size) must be a number, not NA");
			if (size < 0)
				verror(std::isinf(size) ? "options(gmax.data.size) must not be negative (got -Inf)" : "options(gmax.data.size) must not be negative (got %g)", size);
			if (size >= (double)misha::config::MAX_DATA_SIZE_NO_LIMIT)
				m_max_data_size = misha::config::MAX_DATA_SIZE_NO_LIMIT;
			else
				m_max_data_size = size < 1 ? 1 : (uint64_t)size;
		} else
			m_max_data_size = misha::config::MAX_DATA_SIZE_NO_LIMIT;
	}
	return m_max_data_size;
}

uint64_t ConfigurationManager::get_max_mem_usage() const
{
	if (!m_max_mem_usage) {
		SEXP r_max_mem_usage = Rf_GetOption1(Rf_install("gmax.mem.usage"));

		// In KB. As for gmax.data.size: NA (logical or double) and negative values are errors
		// (casting them to uint64_t is undefined), and Inf or any value too large to hold in
		// bytes is no limit. 0 stays allowed: it throttles the kids to one at a time.
		if (Rf_isLogical(r_max_mem_usage) && Rf_length(r_max_mem_usage) >= 1 && LOGICAL(r_max_mem_usage)[0] == NA_LOGICAL)
			verror("options(gmax.mem.usage) must be a number, not NA");

		if (is_num_vector(r_max_mem_usage) && Rf_length(r_max_mem_usage) >= 1) {
			double kb = num_elt(r_max_mem_usage, 0);

			if (std::isnan(kb))
				verror("options(gmax.mem.usage) must be a number, not NA");
			if (kb < 0)
				verror(std::isinf(kb) ? "options(gmax.mem.usage) must not be negative (got -Inf)" : "options(gmax.mem.usage) must not be negative (got %g)", kb);
			m_max_mem_usage = kb * 1000 >= (double)misha::config::UNLIMITED ? misha::config::UNLIMITED : (uint64_t)kb * 1000;
		} else
			m_max_mem_usage = misha::config::UNLIMITED;
	}
	return m_max_mem_usage;
}

uint64_t ConfigurationManager::get_big_intervals_size() const
{
	if (!m_big_intervals_size) {
		SEXP r_big_intervals_size = Rf_GetOption1(Rf_install("gbig.intervals.size"));

		if (Rf_isReal(r_big_intervals_size))
			m_big_intervals_size = (uint64_t)REAL(r_big_intervals_size)[0];
		else if (Rf_isInteger(r_big_intervals_size))
			m_big_intervals_size = INTEGER(r_big_intervals_size)[0];
		else
			m_big_intervals_size = misha::config::UNLIMITED;
		m_big_intervals_size = min(m_big_intervals_size, get_max_data_size());
	}
	return m_big_intervals_size;
}

uint64_t ConfigurationManager::get_quantile_edge_data_size() const
{
	if (!m_quantile_edge_data_size) {
		SEXP r_quantile_edge_data_size = Rf_GetOption1(Rf_install("gquantile.edge.data.size"));

		if (Rf_isReal(r_quantile_edge_data_size))
			m_quantile_edge_data_size = (uint64_t)REAL(r_quantile_edge_data_size)[0];
		else if (Rf_isInteger(r_quantile_edge_data_size))
			m_quantile_edge_data_size = INTEGER(r_quantile_edge_data_size)[0];
		else
			m_quantile_edge_data_size = 0;
	}
	return m_quantile_edge_data_size;
}

double ConfigurationManager::get_multitask_max_records_factor() const
{
	if (m_multitask_max_records_factor <= 0.0) {
		SEXP r_factor = Rf_GetOption1(Rf_install("gmultitask.max.records.factor"));
		double factor = 0.0;

		if (Rf_isReal(r_factor))
			factor = REAL(r_factor)[0];
		else if (Rf_isInteger(r_factor))
			factor = (double)INTEGER(r_factor)[0];
		else
			factor = misha::config::DEFAULT_MULTITASK_MAX_RECORDS_FACTOR;

		if (!std::isfinite(factor) || factor < 1.0)
			factor = 1.0;

		m_multitask_max_records_factor = factor;
	}
	return m_multitask_max_records_factor;
}

uint64_t ConfigurationManager::get_track_chunk_size() const
{
	if (!m_track_chunk_size) {
		SEXP r_track_chunk_size = Rf_GetOption1(Rf_install("gtrack.chunk.size"));

		if (Rf_isReal(r_track_chunk_size))
			m_track_chunk_size = (uint64_t)REAL(r_track_chunk_size)[0];
		else if (Rf_isInteger(r_track_chunk_size))
			m_track_chunk_size = INTEGER(r_track_chunk_size)[0];
		else
			m_track_chunk_size = misha::config::DEFAULT_TRACK_CHUNK_SIZE;
	}
	return m_track_chunk_size;
}

uint64_t ConfigurationManager::get_track_num_chunks() const
{
	if (!m_track_num_chunks) {
		SEXP r_track_num_chunks = Rf_GetOption1(Rf_install("gtrack.num.chunks"));

		if (Rf_isReal(r_track_num_chunks))
			m_track_num_chunks = (uint64_t)REAL(r_track_num_chunks)[0];
		else if (Rf_isInteger(r_track_num_chunks))
			m_track_num_chunks = INTEGER(r_track_num_chunks)[0];
		else
			m_track_num_chunks = misha::config::DEFAULT_TRACK_NUM_CHUNKS;
	}
	return m_track_num_chunks;
}

rdb::MultitaskingMode ConfigurationManager::select_multitasking_mode(bool is_deterministic, uint64_t estimated_size) const
{
	// If multitasking is disabled, use single-threaded mode
	if (!get_multitasking())
		return rdb::MT_MODE_SINGLE;

	// Auto-disable multitasking for small datasets to avoid fork overhead
	// See ConfigurationDefaults.h for rationale on MIN_RECORDS_PER_PROCESS threshold
	uint64_t min_size_for_multitasking = get_max_processes() * misha::config::MIN_RECORDS_PER_PROCESS;
	if (estimated_size < min_size_for_multitasking)
		return rdb::MT_MODE_SINGLE;

	// For deterministic operations, use MMAP mode if size fits in buffer
	if (is_deterministic && estimated_size <= get_max_data_size())
		return rdb::MT_MODE_MMAP;

	// If size is too large, fall back to single-threaded mode
	// (In the future, we might add dynamic buffer allocation or other strategies)
	return rdb::MT_MODE_SINGLE;
}

void ConfigurationManager::verify_max_data_size(uint64_t data_size, const char *data_name, bool check_all_kids,
											  const char *suggestion, const char *unit) const
{
	if (data_size > get_max_data_size()) {
		// gmax.data.size is a number of records (rows of a result, values, bases), not of bytes.
		// data_size is often the number reached so far, not the final one: no value to suggest.
		verror("%s size (at least %llu %s) exceeded the maximum allowed (%llu %s, set by gmax.data.size).\n\n"
			   "Suggestions:\n"
			   "%s"
			   "  - Raise the limit: options(gmax.data.size = <a larger number of %s>), or Inf for no limit\n"
			   "  - Reduce scope (use smaller genomic regions or filters)\n"
			   "  - Process data in smaller chunks\n\n"
			   "Note: gmax.data.size is auto-configured based on your system RAM.\n"
			   "      You can check the current setting with getOption('gmax.data.size')",
			   data_name && *data_name ? data_name : "Result", (unsigned long long)data_size, unit,
			   (unsigned long long)get_max_data_size(), unit, suggestion ? suggestion : "", unit);
	}

	if (check_all_kids)
		update_res_data_size(data_size);
}

void ConfigurationManager::update_res_data_size(uint64_t data_size) const
{
	// This function is defined in rdbutils.cpp
	// It updates cumulative data size tracking for child processes
	rdb::update_res_data_size(data_size);
}
