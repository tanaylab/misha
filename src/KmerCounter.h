#ifndef KMER_COUNTER_H_
#define KMER_COUNTER_H_

#include <string>
#include <memory>

#include "GenomeSeqScorer.h"

class KmerCounter : public GenomeSeqScorer
{
public:
    enum CountMode
    {
        SUM,     // Count total occurrences of kmer
        FRACTION // Calculate fraction of kmer in sequence
    };

    KmerCounter(const std::string &kmer, const std::string &genome_root,
                CountMode mode = SUM, bool extend = true, char strand = 0);

    // Constructor with shared GenomeSeqFetch for caching
    KmerCounter(const std::string &kmer, GenomeSeqFetch* shared_seqfetch,
                CountMode mode = SUM, bool extend = true, char strand = 0);

    // Implement the virtual function from the base class
    float score_interval(const GInterval &interval, const GenomeChromKey &chromkey) override;

private:
    std::string m_kmer;
    CountMode m_mode;

    // Helper method to count kmers in an interval with a specific strand
    size_t count_in_interval(const GInterval &fetch_interval, const GenomeChromKey &chromkey,
                            const GInterval &original_interval, size_t &count);
};

#endif // KMER_COUNTER_H_