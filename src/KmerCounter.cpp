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

double KmerCounter::score_interval(const GInterval &interval, const GenomeChromKey &chromkey)
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

    // Count the forward strand if requested
    if (m_strand == 0 || m_strand == 1)
    {
        // Forward strand processing (strand = 1)
        GInterval forward_interval = fetch_interval;
        forward_interval.strand = 1;

        size_t forward_count = 0;
        count_in_interval(forward_interval, chromkey, interval, forward_count);
        total_count += forward_count;
    }

    // Count the reverse strand if requested
    if (m_strand == 0 || m_strand == -1)
    {
        // Reverse strand processing (strand = -1)
        GInterval reverse_interval = fetch_interval;
        reverse_interval.strand = -1;

        size_t reverse_count = 0;
        count_in_interval(reverse_interval, chromkey, interval, reverse_count);
        total_count += reverse_count;
    }

    // Return appropriate value based on mode
    if (m_mode == FRACTION)
    {
        size_t valid_positions = possible_positions(interval, chromkey);
        return valid_positions > 0 ? static_cast<float>(total_count) / valid_positions : 0.0f;
    }
    else
    {
        // Return the raw count
        return static_cast<float>(total_count);
    }
}

size_t KmerCounter::possible_positions(const GInterval &interval, const GenomeChromKey &chromkey) const
{
    // With extend the read runs k-1 bases past the end (clipped at the
    // chromosome end), so every anchor in the interval can fit; without it the
    // last k-1 anchors cannot.
    const int64_t k = static_cast<int64_t>(m_kmer.length());
    const int64_t read_end = std::min(static_cast<int64_t>(interval.end) + read_extension(),
                                      static_cast<int64_t>(chromkey.get_chrom_size(interval.chromid)));
    const int64_t fit_end = std::min(static_cast<int64_t>(interval.end), read_end - k + 1);
    const size_t per_strand = fit_end > interval.start ? static_cast<size_t>(fit_end - interval.start) : 0;
    return m_strand == 0 ? 2 * per_strand : per_strand;
}

void KmerCounter::count_in_interval(const GInterval &fetch_interval, const GenomeChromKey &chromkey,
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
            return;
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
    }
    catch (TGLException &e)
    {
        count = 0;
    }
    catch (std::exception &e)
    {
        // Catch any other standard exceptions
        count = 0;
    }
}

float KmerCounter::score_from_positions(const std::vector<size_t> &fwd_positions,
                                        const std::vector<size_t> &rev_positions,
                                        const GInterval &original_interval,
                                        const GInterval &fetch_interval,
                                        const GenomeChromKey &chromkey) const
{
    // Calculate the bounds within the fetched sequence that correspond to the original interval
    size_t original_start_pos = 0;
    if (fetch_interval.start < original_interval.start) {
        original_start_pos = original_interval.start - fetch_interval.start;
    }

    size_t original_end_pos = fetch_interval.end - fetch_interval.start;
    if (fetch_interval.end > original_interval.end) {
        original_end_pos = original_end_pos - (fetch_interval.end - original_interval.end);
    }

    size_t kmer_len = m_kmer.length();

    // Count positions based on strand parameter
    size_t total_count = 0;

    // Count forward strand if needed
    if (m_strand == 0 || m_strand == 1) {
        size_t fwd_count = 0;
        for (size_t pos : fwd_positions) {
            if (pos >= original_start_pos && pos < original_end_pos) {
                fwd_count++;
            }
        }
        total_count += fwd_count;
    }

    // Count reverse strand if needed
    if (m_strand == 0 || m_strand == -1) {
        // rev_positions index the reverse complement of fetch_interval. Anchor a
        // match at its genomic start, forward offset fetch_len - pos - kmer_len,
        // as the forward strand and count_in_interval() do.
        const size_t fetch_len = fetch_interval.end - fetch_interval.start;
        size_t rev_count = 0;
        for (size_t pos : rev_positions) {
            const size_t fwd_pos = fetch_len - pos - kmer_len;
            if (fwd_pos >= original_start_pos && fwd_pos < original_end_pos) {
                rev_count++;
            }
        }
        total_count += rev_count;
    }

    // Return based on mode
    if (m_mode == FRACTION) {
        const size_t total_possible_positions = possible_positions(original_interval, chromkey);
        return total_possible_positions > 0 ? static_cast<float>(total_count) / total_possible_positions : 0.0f;
    } else {
        return static_cast<float>(total_count);
    }
}