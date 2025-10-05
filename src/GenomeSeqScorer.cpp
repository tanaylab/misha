#include "GenomeSeqScorer.h"
#include <algorithm>

GenomeSeqScorer::GenomeSeqScorer(const std::string &genome_root, bool extend, char strand)
    : m_extend(extend), m_strand(strand), m_seqfetch_ptr(&m_seqfetch_owned)
{
    m_seqfetch_owned.set_seqdir(genome_root + "/seq");
}

GenomeSeqScorer::GenomeSeqScorer(GenomeSeqFetch* shared_seqfetch, bool extend, char strand)
    : m_extend(extend), m_strand(strand), m_seqfetch_ptr(shared_seqfetch)
{
    // Using shared seqfetch, no need to initialize m_seqfetch_owned
}

GInterval GenomeSeqScorer::calculate_expanded_interval(const GInterval &interval, const GenomeChromKey &chromkey, int64_t pattern_length)
{
    GInterval expanded_interval = interval;

    // If extend is true, extend the interval to allow for patterns that span the boundary
    if (m_extend && pattern_length > 1)
    {
        int64_t pattern_length_minus_one = static_cast<int64_t>(pattern_length - 1);

        if (m_strand == -1)
        {
            // For reverse strand, we need to expand the start
            int64_t start_pos = static_cast<int64_t>(expanded_interval.start);
            int64_t zero = 0;

            expanded_interval.start = static_cast<decltype(expanded_interval.start)>(
                std::max(start_pos - pattern_length_minus_one, zero));
        }
        else if (m_strand == 1)
        {
            // For forward strand, we need to expand the end
            int64_t end_pos = static_cast<int64_t>(expanded_interval.end);
            int64_t chrom_size = static_cast<int64_t>(chromkey.get_chrom_size(interval.chromid));

            expanded_interval.end = static_cast<decltype(expanded_interval.end)>(
                std::min(end_pos + pattern_length_minus_one, chrom_size));
        }
        else
        {
            // For both strands, expand both directions
            int64_t start_pos = static_cast<int64_t>(expanded_interval.start);
            int64_t end_pos = static_cast<int64_t>(expanded_interval.end);
            int64_t zero = 0;
            int64_t chrom_size = static_cast<int64_t>(chromkey.get_chrom_size(interval.chromid));

            expanded_interval.start = static_cast<decltype(expanded_interval.start)>(
                std::max(start_pos - pattern_length_minus_one, zero));

            expanded_interval.end = static_cast<decltype(expanded_interval.end)>(
                std::min(end_pos + pattern_length_minus_one, chrom_size));
        }
    }

    return expanded_interval;
}