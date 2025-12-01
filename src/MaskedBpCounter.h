#ifndef MASKED_BP_COUNTER_H_
#define MASKED_BP_COUNTER_H_

#include <string>
#include "GenomeSeqScorer.h"

class MaskedBpCounter : public GenomeSeqScorer
{
public:
    enum CountMode
    {
        COUNT,    // Count total masked base pairs
        FRACTION  // Calculate fraction of masked base pairs
    };

    MaskedBpCounter(const std::string &genome_root, CountMode mode = COUNT);

    // Constructor with shared GenomeSeqFetch for caching
    MaskedBpCounter(GenomeSeqFetch* shared_seqfetch, CountMode mode = COUNT);

    // Implement the virtual function from the base class
    float score_interval(const GInterval &interval, const GenomeChromKey &chromkey) override;

    CountMode get_mode() const { return m_mode; }

private:
    CountMode m_mode;
};

#endif // MASKED_BP_COUNTER_H_
