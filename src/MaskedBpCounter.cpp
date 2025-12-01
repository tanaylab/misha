#include "MaskedBpCounter.h"
#include <algorithm>
#include <cctype>
#include <limits>

MaskedBpCounter::MaskedBpCounter(const std::string &genome_root, CountMode mode)
    : GenomeSeqScorer(genome_root, false, 0), m_mode(mode)
{
    // Note: extend=false, strand=0 passed to base class
    // Masking is position-specific, no extension needed
}

MaskedBpCounter::MaskedBpCounter(GenomeSeqFetch* shared_seqfetch, CountMode mode)
    : GenomeSeqScorer(shared_seqfetch, false, 0), m_mode(mode)
{
    // Note: extend=false, strand=0 passed to base class
}

float MaskedBpCounter::score_interval(const GInterval &interval, const GenomeChromKey &chromkey)
{
    // Read sequence for the exact interval (no extension)
    std::vector<char> seq;

    try
    {
        m_seqfetch_ptr->read_interval(interval, chromkey, seq);

        if (seq.empty())
        {
            return 0.0f;
        }

        // Count lowercase (masked) base pairs
        size_t masked_count = 0;
        for (char c : seq)
        {
            if (std::islower(static_cast<unsigned char>(c)))
            {
                masked_count++;
            }
        }

        // Return based on mode
        if (m_mode == FRACTION)
        {
            // Return fraction of masked bases
            return seq.size() > 0 ? static_cast<float>(masked_count) / seq.size() : 0.0f;
        }
        else
        {
            // Return raw count
            return static_cast<float>(masked_count);
        }
    }
    catch (TGLException &e)
    {
        // Return NaN on error
        return std::numeric_limits<float>::quiet_NaN();
    }
    catch (std::exception &e)
    {
        return std::numeric_limits<float>::quiet_NaN();
    }
}
