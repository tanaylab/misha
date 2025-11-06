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
    // For all strand modes, we extend END only because:
    // - At genomic position i, we always need sequence [i, i+motif_len)
    // - For reverse strand, we compute RC of [i, i+motif_len) at the same genomic positions
    // - Therefore, a motif at position (end-1) needs sequence up to (end-1)+motif_len
    if (m_extend && pattern_length > 1)
    {
        int64_t pattern_length_minus_one = static_cast<int64_t>(pattern_length - 1);
        int64_t end_pos = static_cast<int64_t>(expanded_interval.end);
        int64_t chrom_size = static_cast<int64_t>(chromkey.get_chrom_size(interval.chromid));

        expanded_interval.end = static_cast<decltype(expanded_interval.end)>(
            std::min(end_pos + pattern_length_minus_one, chrom_size));
    }

    return expanded_interval;
}