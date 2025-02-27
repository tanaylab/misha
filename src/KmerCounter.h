#ifndef KMER_COUNTER_H_
#define KMER_COUNTER_H_

#include <string>
#include <memory>
#include <algorithm>
#include <cctype>
#include <limits>
#include <R.h>
#include <Rinternals.h>
#include "rdbutils.h"
#include "GenomeSeqFetch.h"
#include "GInterval.h"
#include "GenomeChromKey.h"

class KmerCounter
{
public:
    enum CountMode
    {
        SUM,     // Count total occurrences of kmer
        FRACTION // Calculate fraction of kmer in sequence
    };

    KmerCounter(const std::string &kmer, const std::string &genome_root, CountMode mode = SUM)
        : m_kmer(kmer), m_mode(mode)
    {
        // Validate kmer
        if (m_kmer.empty())
        {
            rdb::verror("Kmer string cannot be empty");
        }

        // Convert kmer to uppercase
        std::transform(m_kmer.begin(), m_kmer.end(), m_kmer.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });

        m_seqfetch.set_seqdir(genome_root + "/seq");
    }

    // Count kmers in a given interval
    float count_interval(const GInterval &interval, const GenomeChromKey &chromkey)
    {
        // Safety check - shouldn't happen due to constructor validation
        if (m_kmer.empty())
        {
            return std::numeric_limits<float>::quiet_NaN();
        }

        // Fetch the sequence for the interval
        std::vector<char> seq;
        try
        {
            m_seqfetch.read_interval(interval, chromkey, seq);

            // If sequence is too short to contain even one kmer
            if (seq.size() < m_kmer.length())
            {
                return (m_mode == FRACTION) ? 0.0f : 0.0f;
            }

            std::string target(seq.begin(), seq.end());

            // Convert target to uppercase for case-insensitive matching
            std::transform(target.begin(), target.end(), target.begin(),
                           [](unsigned char c)
                           { return std::toupper(c); });

            // Count occurrences of kmer in the sequence
            size_t count = 0;
            size_t pos = 0;

            while ((pos = target.find(m_kmer, pos)) != std::string::npos)
            {
                count++;
                pos += 1; // Move to the next possible match position
            }

            if (m_mode == FRACTION)
            {
                // Number of possible positions where kmers can start
                size_t possible_positions = target.length() - m_kmer.length() + 1;
                return static_cast<float>(count) / possible_positions;
            }
            else
            {
                // Return the raw count
                return static_cast<float>(count);
            }
        }
        catch (TGLException &e)
        {
            return std::numeric_limits<float>::quiet_NaN();
        }
        catch (std::exception &e)
        {
            // Catch any other standard exceptions
            return std::numeric_limits<float>::quiet_NaN();
        }
    }

private:
    std::string m_kmer;
    CountMode m_mode;
    GenomeSeqFetch m_seqfetch;
};

#endif // KMER_COUNTER_H_