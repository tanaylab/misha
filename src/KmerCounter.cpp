#include "KmerCounter.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>

KmerCounter::KmerCounter(const std::string &kmer, const std::string &genome_root,
                        CountMode mode, bool extend, char strand)
    : GenomeSeqScorer(genome_root, extend, strand), m_kmer(kmer), m_mode(mode)
{
    // Validate kmer
    if (m_kmer.empty())
    {
        rdb::verror("Kmer string cannot be empty");
    }

    // Convert kmer to uppercase
    std::transform(m_kmer.begin(), m_kmer.end(), m_kmer.begin(),
                   [](unsigned char c) { return std::toupper(c); });
}

KmerCounter::KmerCounter(const std::string &kmer, GenomeSeqFetch* shared_seqfetch,
                        CountMode mode, bool extend, char strand)
    : GenomeSeqScorer(shared_seqfetch, extend, strand), m_kmer(kmer), m_mode(mode)
{
    // Validate kmer
    if (m_kmer.empty())
    {
        rdb::verror("Kmer string cannot be empty");
    }

    // Convert kmer to uppercase
    std::transform(m_kmer.begin(), m_kmer.end(), m_kmer.begin(),
                   [](unsigned char c) { return std::toupper(c); });
}

float KmerCounter::score_interval(const GInterval &interval, const GenomeChromKey &chromkey)
{
    // Safety check - shouldn't happen due to constructor validation
    if (m_kmer.empty())
    {
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Calculate appropriate interval based on extension setting
    GInterval fetch_interval = calculate_expanded_interval(interval, chromkey, m_kmer.length());

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
        if (m_strand == 0)
        {
            // For both strands, we need to multiply by 2
            valid_positions *= 2;
        }
        return valid_positions > 0 ? static_cast<float>(total_count) / valid_positions : 0.0f;
    }
    else
    {
        // Return the raw count
        return static_cast<float>(total_count);
    }
}

size_t KmerCounter::count_in_interval(const GInterval &fetch_interval, const GenomeChromKey &chromkey,
                                     const GInterval &original_interval, size_t &count)
{
    count = 0;
    std::vector<char> seq;
    try
    {
        m_seqfetch_ptr->read_interval(fetch_interval, chromkey, seq);

        // If sequence is too short to contain even one kmer
        if (seq.size() < m_kmer.length())
        {
            return 0;
        }

        std::string target(seq.begin(), seq.end());

        // Convert target to uppercase for case-insensitive matching
        std::transform(target.begin(), target.end(), target.begin(),
                       [](unsigned char c) { return std::toupper(c); });

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
        const char* target_data = target.data();
        const char* kmer_data = m_kmer.data();
        const size_t kmer_len = m_kmer.length();

        for (size_t pos = original_start_pos;
             pos < original_end_pos && pos <= target.length() - kmer_len;
             pos++)
        {
            // Use memcmp for faster fixed-length comparison
            if (std::memcmp(target_data + pos, kmer_data, kmer_len) == 0)
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