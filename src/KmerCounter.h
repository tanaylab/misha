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

    KmerCounter(const std::string &kmer, const std::string &genome_root, CountMode mode = SUM, bool extend = true, char strand = 0)
        : m_kmer(kmer), m_mode(mode), m_extend(extend), m_strand(strand)
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

        // Calculate total counts from both strands if necessary
        size_t total_count = 0;
        size_t possible_positions = 0;

        // Count the forward strand if requested
        if (m_strand == 0 || m_strand == 1)
        {
            // Forward strand processing (strand = 1)
            GInterval forward_interval = fetch_interval;
            forward_interval.strand = 1;

            size_t forward_count = 0;
            possible_positions = count_in_interval(forward_interval, chromkey, interval, forward_count);
            total_count += forward_count;
        }

        // Count the reverse strand if requested
        if (m_strand == 0 || m_strand == -1)
        {
            // Reverse strand processing (strand = -1)
            GInterval reverse_interval = fetch_interval;
            reverse_interval.strand = -1;

            size_t reverse_count = 0;
            possible_positions = count_in_interval(reverse_interval, chromkey, interval, reverse_count);
            total_count += reverse_count;
        }

        // Return appropriate value based on mode
        if (m_mode == FRACTION)
        {
            // Return the fraction of positions that contain the kmer
            // Denominator is adjusted to count only positions where a kmer could fully fit
            size_t valid_positions = possible_positions > m_kmer.length() - 1 ? possible_positions - (m_kmer.length() - 1) : 0;
            return valid_positions > 0 ? static_cast<float>(total_count) / valid_positions : 0.0f;
        }
        else
        {
            // Return the raw count
            return static_cast<float>(total_count);
        }
    }

private:
    std::string m_kmer;
    CountMode m_mode;
    bool m_extend;
    char m_strand; // 0 = both strands, 1 = forward only, -1 = reverse only
    GenomeSeqFetch m_seqfetch;

    // Helper method to count kmers in an interval with a specific strand
    size_t count_in_interval(const GInterval &fetch_interval, const GenomeChromKey &chromkey,
                             const GInterval &original_interval, size_t &count)
    {
        count = 0;
        std::vector<char> seq;
        try
        {
            m_seqfetch.read_interval(fetch_interval, chromkey, seq);

            // If sequence is too short to contain even one kmer
            if (seq.size() < m_kmer.length())
            {
                return 0;
            }

            std::string target(seq.begin(), seq.end());

            // Convert target to uppercase for case-insensitive matching
            std::transform(target.begin(), target.end(), target.begin(),
                           [](unsigned char c)
                           { return std::toupper(c); });

            // Calculate the relative position in the fetched sequence that corresponds
            // to the start of the original interval
            size_t original_start_pos = 0;
            if (fetch_interval.start < original_interval.start)
            {
                original_start_pos = original_interval.start - fetch_interval.start;
            }

            // Calculate the end position in the fetched sequence that corresponds
            // to the end of the original interval
            size_t original_end_pos = target.length();
            if (fetch_interval.end > original_interval.end)
            {
                original_end_pos = original_end_pos - (fetch_interval.end - original_interval.end);
            }

            // Ensure positions are within valid range
            original_end_pos = std::min(original_end_pos, target.length());

            // Calculate the number of possible positions where a kmer can start
            // within the original interval (not the extended one)
            size_t possible_positions = original_end_pos > original_start_pos ? original_end_pos - original_start_pos : 0;

            // Count occurrences of kmer in the sequence
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

            return possible_positions;
        }
        catch (TGLException &e)
        {
            return 0;
        }
        catch (std::exception &e)
        {
            // Catch any other standard exceptions
            return 0;
        }
    }
};

#endif // KMER_COUNTER_H_