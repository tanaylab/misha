#ifndef CHAININTERVALCONVERTER_H_
#define CHAININTERVALCONVERTER_H_

#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <string>

// Forward declarations
namespace rdb {
	class IntervUtils;
	class ChainIntervals;
}

// Converter class for chain interval conversions between R and C++
class ChainIntervalConverter {
public:
	// Constructor takes a reference to IntervUtils for accessing conversion methods
	explicit ChainIntervalConverter(rdb::IntervUtils &iu);

	// Converts R chain intervals (data frame) to C++ ChainIntervals
	// rchain: R data frame with chain interval columns
	// chain_intervs: output vector of chain intervals
	// src_id2chrom: output vector mapping source chromosome IDs to names
	void convert_rchain_intervs(SEXP rchain, rdb::ChainIntervals &chain_intervs, std::vector<std::string> &src_id2chrom);

	// Converts C++ ChainIntervals to R data frame
	// chain_intervs: input vector of chain intervals
	// src_id2chrom: input vector mapping source chromosome IDs to names
	// Returns: R data frame with chain interval columns
	SEXP convert_chain_intervs(const rdb::ChainIntervals &chain_intervs, std::vector<std::string> &src_id2chrom);

private:
	rdb::IntervUtils &m_iu;
};

#endif /* CHAININTERVALCONVERTER_H_ */

