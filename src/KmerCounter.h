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

    KmerCounter(const std::string &kmer, const std::string &genome_root, CountMode mode = SUM, bool extend = true)
        : m_kmer(kmer), m_mode(mode), m_extend(extend)
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

        // Calculate appropriate interval based on extension setting
        GInterval fetch_interval = interval;

        // If extend is true, extend the interval to allow for kmers that span the boundary
        if (m_extend && m_kmer.length() > 1)
        {
            // Explicitly cast both arguments to the same type (int64_t) to avoid template deduction issues
            // Using static_cast to ensure proper conversion
            int64_t kmer_length_minus_one = static_cast<int64_t>(m_kmer.length() - 1);
            int64_t start_pos = static_cast<int64_t>(fetch_interval.start);
            int64_t zero = 0;

            fetch_interval.start = static_cast<decltype(fetch_interval.start)>(
                std::max(start_pos - kmer_length_minus_one, zero));

            // For end coordinate, also ensure both arguments are same type
            int64_t end_pos = static_cast<int64_t>(fetch_interval.end);
            int64_t chrom_size = static_cast<int64_t>(chromkey.get_chrom_size(interval.chromid));

            fetch_interval.end = static_cast<decltype(fetch_interval.end)>(
                std::min(end_pos + kmer_length_minus_one, chrom_size));
        }

        // Fetch the sequence for the interval
        std::vector<char> seq;
        try
        {
            m_seqfetch.read_interval(fetch_interval, chromkey, seq);

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

            // Calculate the relative position in the fetched sequence that corresponds
            // to the start of the original interval
            size_t original_start_pos = 0;
            if (fetch_interval.start < interval.start)
            {
                original_start_pos = interval.start - fetch_interval.start;
            }

            // Calculate the end position in the fetched sequence that corresponds
            // to the end of the original interval
            size_t original_end_pos = target.length();
            if (fetch_interval.end > interval.end)
            {
                original_end_pos = original_end_pos - (fetch_interval.end - interval.end);
            }

            // Ensure positions are within valid range
            original_end_pos = std::min(original_end_pos, target.length());

            // Calculate the number of possible positions where a kmer can start
            // within the original interval (not the extended one)
            size_t possible_positions = original_end_pos > original_start_pos ? original_end_pos - original_start_pos : 0;

            // Only count kmers whose start position falls within the original interval
            for (size_t pos = original_start_pos;
                 pos < original_end_pos && pos <= target.length() - m_kmer.length();
                 pos++)
            {
                if (target.compare(pos, m_kmer.length(), m_kmer) == 0)
                {
                    count++;
                }
            }

            if (m_mode == FRACTION)
            {
                // Return the fraction of positions that contain the kmer
                // Denominator is adjusted to count only positions where a kmer could fully fit
                size_t valid_positions = possible_positions > m_kmer.length() - 1 ? possible_positions - (m_kmer.length() - 1) : 0;
                return valid_positions > 0 ? static_cast<float>(count) / valid_positions : 0.0f;
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
    bool m_extend;
    GenomeSeqFetch m_seqfetch;
};

#endif // KMER_COUNTER_H_