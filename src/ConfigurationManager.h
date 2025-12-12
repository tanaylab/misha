#ifndef CONFIGURATIONMANAGER_H_
#define CONFIGURATIONMANAGER_H_

#include <cstdint>
#include <R.h>
#include <Rinternals.h>

// Undefine R macros that conflict with C++ standard library
#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif
#ifdef warning
#undef warning
#endif

namespace rdb {
	enum MultitaskingMode : int;
}

class ConfigurationManager {
public:
	ConfigurationManager(SEXP envir) : m_envir(envir) {}

	// Returns true if multitasking is switched on
	bool get_multitasking() const;

	// Returns absolute maximal number of concurrently opened processes for parallel computation
	uint64_t get_max_processes() const;

	// Returns the maximal number of concurrently opened processes per core for parallel computation
	uint64_t get_max_processes2core() const;

	// Returns minimal scope range per process for parallel computation
	uint64_t get_min_scope4process() const;

	// Returns minimal sequence workload per process for parallel computation (for gseq.pwm, gseq.kmer)
	uint64_t get_min_seqs_work4process() const;

	// Returns the upper limit for data size
	uint64_t get_max_data_size() const;

	// Returns the upper limit for memory usage
	uint64_t get_max_mem_usage() const;

	// Returns the threshold for creating a big intervals set
	uint64_t get_big_intervals_size() const;

	// Returns the size of the buffer used to store highest/lowest values for high-precision computation of quantiles
	uint64_t get_quantile_edge_data_size() const;

	// Returns the chunk size of 2D track
	uint64_t get_track_chunk_size() const;

	// Returns the number of chunks of 2D track
	uint64_t get_track_num_chunks() const;

	// Selects the appropriate multitasking mode based on estimated result size
	// is_deterministic: true if result size can be known precisely before running
	// estimated_size: estimated number of result records (intervals, values, etc.)
	rdb::MultitaskingMode select_multitasking_mode(bool is_deterministic, uint64_t estimated_size) const;

	// Verifies that the data size does not exceed the maximum allowed.
	// If check_all_kids == true then the limit is checked cumulatively for all child processes.
	void verify_max_data_size(uint64_t data_size, const char *data_name = "Result", bool check_all_kids = true) const;

private:
	[[maybe_unused]] SEXP m_envir;

	// Mutable cache for lazy-loaded configuration values
	mutable int                   m_multitasking{-1};
	mutable uint64_t              m_max_data_size{0};
	mutable uint64_t              m_max_mem_usage{0};
	mutable uint64_t              m_big_intervals_size{0};
	mutable uint64_t              m_max_processes{0};
	mutable uint64_t              m_max_processes2core{0};
	mutable uint64_t              m_min_scope4process{0};
	mutable uint64_t              m_min_seqs_work4process{0};
	mutable uint64_t              m_quantile_edge_data_size{0};
	mutable uint64_t              m_track_chunk_size{0};
	mutable uint64_t              m_track_num_chunks{0};

	// Helper function to update cumulative data size for child processes
	void update_res_data_size(uint64_t data_size) const;
};

#endif /* CONFIGURATIONMANAGER_H_ */

