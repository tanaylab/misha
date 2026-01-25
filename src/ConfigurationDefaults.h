/*
 * ConfigurationDefaults.h
 *
 * Centralized configuration defaults for the misha package.
 * This file defines the single source of truth for all default configuration values
 * used by the C++ layer.
 *
 * NOTE: If you modify a default here, ensure the corresponding value in
 * R/config-defaults.R is also updated to maintain consistency between R and C++.
 */

#ifndef CONFIGURATIONDEFAULTS_H_
#define CONFIGURATIONDEFAULTS_H_

#include <cstdint>
#include <limits>

namespace misha {
namespace config {

// ============================================================================
// Default values (must match R/config-defaults.R)
// ============================================================================

// Threshold for converting interval sets to disk-based "big" format (1 million)
constexpr uint64_t DEFAULT_BIG_INTERVALS_SIZE = 1000000;

// Max memory for child processes (in KB = 10 GB)
constexpr uint64_t DEFAULT_MAX_MEM_USAGE_KB = 10000000;

// Evaluation buffer size for vectorized R expression evaluation
constexpr uint64_t DEFAULT_BUF_SIZE = 1000;

// Min genomic scope per parallel process (base pairs)
constexpr uint64_t DEFAULT_MIN_SCOPE4PROCESS = 10000;

// Min sequence workload per parallel process (for gseq.pwm, gseq.kmer)
constexpr uint64_t DEFAULT_MIN_SEQS_WORK4PROCESS = 100000;

// Buffer size for high-precision quantile edge values
constexpr uint64_t DEFAULT_QUANTILE_EDGE_DATA_SIZE = 100000;

// Chunk size for 2D tracks (base pairs)
constexpr uint64_t DEFAULT_TRACK_CHUNK_SIZE = 100000;

// Number of chunks for 2D tracks (0 = auto)
constexpr uint64_t DEFAULT_TRACK_NUM_CHUNKS = 0;

// Max concurrent processes (absolute limit)
constexpr uint64_t DEFAULT_MAX_PROCESSES = 64;

// Max processes per CPU core
constexpr uint64_t DEFAULT_MAX_PROCESSES2CORE = 4;

// ============================================================================
// Thresholds and multipliers (with documented rationale)
// ============================================================================

/*
 * Minimum records per process before multitasking overhead exceeds benefit.
 *
 * Rationale:
 * - Fork overhead is approximately ~1ms per child process
 * - Processing 1000 records takes approximately ~10ms
 * - This gives ~10% overhead, which is an acceptable tradeoff for parallelism
 *
 * Empirical validation: Testing with 89 processes was slower at 20K total records,
 * confirming that the dynamic threshold (max_processes * MIN_RECORDS_PER_PROCESS)
 * effectively prevents unnecessary parallelization overhead.
 */
constexpr uint64_t MIN_RECORDS_PER_PROCESS = 1000;

/*
 * Multiplier for suggested size in error messages (1.5 = 50% headroom).
 *
 * When a user exceeds gmax.data.size, we suggest a new value that is 50% larger
 * than the actual needed size. This provides headroom for:
 * - Slight variations in data size between runs
 * - Additional overhead from metadata and temporary structures
 */
constexpr double SUGGESTED_SIZE_MULTIPLIER = 1.5;

/*
 * Sentinel value for "unlimited" or "unknown" size.
 * Used when:
 * - An R option is not set (defaults to unlimited)
 * - Size estimation cannot determine a valid estimate
 */
constexpr uint64_t UNLIMITED = std::numeric_limits<uint64_t>::max();

}  // namespace config
}  // namespace misha

#endif /* CONFIGURATIONDEFAULTS_H_ */
