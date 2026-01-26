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
	if (m_multitasking < 0) {
		SEXP r_multitasking = Rf_GetOption1(Rf_install("gmultitasking"));

		if (Rf_isLogical(r_multitasking))
			m_multitasking = (int)LOGICAL(r_multitasking)[0];
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

		if (Rf_isReal(r_max_data_size))
			m_max_data_size = (uint64_t)REAL(r_max_data_size)[0];
		else if (Rf_isInteger(r_max_data_size))
			m_max_data_size = INTEGER(r_max_data_size)[0];
		else
			m_max_data_size = misha::config::UNLIMITED;
	}
	return m_max_data_size;
}

uint64_t ConfigurationManager::get_max_mem_usage() const
{
	if (!m_max_mem_usage) {
		SEXP r_max_mem_usage = Rf_GetOption1(Rf_install("gmax.mem.usage"));

		if (Rf_isReal(r_max_mem_usage))
			m_max_mem_usage = (uint64_t)REAL(r_max_mem_usage)[0] * 1000;
		else if (Rf_isInteger(r_max_mem_usage))
			m_max_mem_usage = INTEGER(r_max_mem_usage)[0] * 1000;
		else
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

void ConfigurationManager::verify_max_data_size(uint64_t data_size, const char *data_name, bool check_all_kids) const
{
	if (data_size > get_max_data_size()) {
		// Calculate a suggested size with headroom (see ConfigurationDefaults.h)
		uint64_t suggested_size = (uint64_t)(data_size * misha::config::SUGGESTED_SIZE_MULTIPLIER);

		// Format sizes as human-readable strings
		char current_size_str[100], needed_size_str[100], suggested_size_str[100];
		if (get_max_data_size() >= 1e9)
			snprintf(current_size_str, sizeof(current_size_str), "%.1f GB", get_max_data_size() / 1e9);
		else if (get_max_data_size() >= 1e6)
			snprintf(current_size_str, sizeof(current_size_str), "%.1f MB", get_max_data_size() / 1e6);
		else
			snprintf(current_size_str, sizeof(current_size_str), "%lu bytes", (unsigned long)get_max_data_size());

		if (data_size >= 1e9)
			snprintf(needed_size_str, sizeof(needed_size_str), "%.1f GB", data_size / 1e9);
		else if (data_size >= 1e6)
			snprintf(needed_size_str, sizeof(needed_size_str), "%.1f MB", data_size / 1e6);
		else
			snprintf(needed_size_str, sizeof(needed_size_str), "%lu bytes", (unsigned long)data_size);

		if (suggested_size >= 1e9)
			snprintf(suggested_size_str, sizeof(suggested_size_str), "%.0f", (double)suggested_size);
		else if (suggested_size >= 1e6)
			snprintf(suggested_size_str, sizeof(suggested_size_str), "%.0f", (double)suggested_size);
		else
			snprintf(suggested_size_str, sizeof(suggested_size_str), "%lu", (unsigned long)suggested_size);

		verror("%s size (%s) exceeded the maximum allowed (%s).\n\n"
			   "Suggestions:\n"
			   "  1. Increase buffer size:\n"
			   "     options(gmax.data.size = %s)\n"
			   "  2. Reduce scope (use smaller genomic regions or filters)\n"
			   "  3. Process data in smaller chunks\n\n"
			   "Note: gmax.data.size is auto-configured based on your system RAM.\n"
			   "      You can check the current setting with getOption('gmax.data.size')",
			   data_name, needed_size_str, current_size_str, suggested_size_str);
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
