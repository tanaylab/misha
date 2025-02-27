#ifndef GENOME_SEQ_SCORER_H_
#define GENOME_SEQ_SCORER_H_

#include <string>
#include <memory>
#include <limits>

#include "rdbutils.h"

#include "GenomeSeqFetch.h"
#include "GInterval.h"
#include "GenomeChromKey.h"

// Base class for genomic sequence scoring tools
class GenomeSeqScorer {
public:
    GenomeSeqScorer(const std::string& genome_root, bool extend = true, char strand = 0);
    virtual ~GenomeSeqScorer() = default;

    // Pure virtual function to be implemented by derived classes
    virtual float score_interval(const GInterval& interval, const GenomeChromKey& chromkey) = 0;

protected:
    // Helper method to calculate appropriate interval based on extension setting
    GInterval calculate_expanded_interval(const GInterval& interval, const GenomeChromKey& chromkey, int64_t pattern_length);
    
    bool m_extend;
    char m_strand; // 0 = both strands, 1 = forward only, -1 = reverse only
    GenomeSeqFetch m_seqfetch;
};

#endif // GENOME_SEQ_SCORER_H_