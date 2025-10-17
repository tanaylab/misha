#include <errno.h>
#include <sys/stat.h>

#include "port.h"

#include "GenomeSeqFetch.h"
#include "GenomeIndex.h"

GenomeSeqFetch::GenomeSeqFetch()
    : m_cur_chromid(-1), m_cache_valid(false), m_indexed_mode(false), m_index(nullptr) {}

GenomeSeqFetch::~GenomeSeqFetch() {
    if (m_index) {
        delete m_index;
    }
}

void GenomeSeqFetch::set_seqdir(const std::string &dir) {
    m_seqdir = dir;

    // Check for indexed format (genome.idx exists)
    std::string idx_path = m_seqdir + "/genome.idx";
    struct stat st;

    if (stat(idx_path.c_str(), &st) == 0) {
        // Indexed mode: load index and open main genome.seq file
        m_indexed_mode = true;
        m_index = new GenomeIndex();
        m_index->load(idx_path);

        // Open the consolidated genome.seq file
        std::string seq_path = m_seqdir + "/genome.seq";
        m_bfile.open(seq_path.c_str(), "rb");
        if (m_bfile.error()) {
            TGLError<GenomeSeqFetch>(FILE_READ_FAILED,
                "Failed to open genome.seq: %s", strerror(errno));
        }
    } else {
        // Legacy per-chromosome mode
        m_indexed_mode = false;
        m_index = nullptr;
    }
}

void GenomeSeqFetch::read_interval(const GInterval &interval, const GenomeChromKey &chromkey, vector<char> &result)
{
	// Check cache first
	if (m_cache_valid &&
	    m_cached_interval.chromid == interval.chromid &&
	    m_cached_interval.start == interval.start &&
	    m_cached_interval.end == interval.end &&
	    m_cached_interval.strand == interval.strand) {
		// Cache hit! Return cached sequence
		result = m_cached_seq;
		return;
	}

	result.clear();

	if (m_indexed_mode) {
		// INDEXED MODE: Read from consolidated genome.seq using index
		const ContigIndexEntry *entry = m_index->get_entry(interval.chromid);
		if (!entry) {
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED,
				"Contig ID %d not found in genome index", interval.chromid);
		}

		// Validate interval bounds
		if (interval.start < 0 || interval.end < 0) {
			TGLError<GenomeSeqFetch>(INVALID_INTERVAL,
				"Invalid interval coordinates: start=%ld end=%ld",
				interval.start, interval.end);
		}

		if (interval.end > (int64_t)entry->length) {
			TGLError<GenomeSeqFetch>(INVALID_INTERVAL,
				"Interval end %ld exceeds contig length %lu for contig %s",
				interval.end, entry->length, entry->name.c_str());
		}

		// Calculate absolute position in genome.seq
		uint64_t abs_offset = entry->offset + interval.start;
		int64_t size = interval.end - interval.start;

		if (size <= 0) {
			if (size == 0)
				size = 1;
			else
				return;
		}

		// Seek and read from consolidated file
		result.resize(size);
		m_bfile.seek(abs_offset, SEEK_SET);
		if (m_bfile.read(&*result.begin(), result.size()) != result.size()) {
			if (m_bfile.error())
				TGLError<GenomeSeqFetch>(FILE_READ_FAILED,
					"Failed to read sequence data from genome.seq: %s", strerror(errno));
			TGLError<GenomeSeqFetch>(FILE_READ_FAILED,
				"Failed to read sequence data from genome.seq");
		}

	} else {
		// LEGACY MODE: Per-chromosome files
		if (m_cur_chromid != interval.chromid) {
			char filename[PATH_MAX];
			const char *chrom_name = chromkey.id2chrom(interval.chromid).c_str();

			m_cur_chromid = interval.chromid;
			m_cache_valid = false;  // Invalidate cache on chromosome change

			// Try chromosome name as-is first
			snprintf(filename, sizeof(filename), "%s/%s.seq", m_seqdir.c_str(), chrom_name);
			m_bfile.close();
			m_bfile.open(filename, "rb");

			// If file doesn't exist and name doesn't start with "chr", try with "chr" prefix (backward compatibility)
			if (m_bfile.error() && strncmp(chrom_name, "chr", 3) != 0) {
				snprintf(filename, sizeof(filename), "%s/chr%s.seq", m_seqdir.c_str(), chrom_name);
				m_bfile.open(filename, "rb");
			}

			if (m_bfile.error())
				TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", filename, strerror(errno));
		}

		interval.verify(chromkey, false);

		int64_t size = min(interval.end, m_bfile.file_size()) - interval.start;

		if (size < 0)
			return;

		if (!size)
			size++;

		// Resize output buffer to required size; reserve a bit more to avoid reallocations
		result.resize(size);
		m_bfile.seek(interval.start, SEEK_SET);
		if (m_bfile.read(&*result.begin(), result.size()) != result.size()) {
			if (m_bfile.error())
				TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed: %s", m_bfile.file_name().c_str(), strerror(errno));

			TGLError<GenomeSeqFetch>(FILE_READ_FAILED, "Reading sequence file %s failed", m_bfile.file_name().c_str());
		}
	}

	// Apply reverse-complement if needed (same for both modes)
	if (interval.strand == -1) {
		size_t len = result.size();
		for (size_t i = 0; i < len / 2; ++i) {
			char tmp = basepair2complementary(result[len - 1 - i]);
			result[len - 1 - i] = basepair2complementary(result[i]);
			result[i] = tmp;
		}
		// Handle middle element for odd-length sequences
		if (len % 2 == 1) {
			result[len / 2] = basepair2complementary(result[len / 2]);
		}
	}

	// Update cache only for reasonably small sequences to avoid expensive copies
	if (result.size() <= 131072) { // <=128KB
		m_cached_interval = interval;
		m_cached_seq = result;
		m_cache_valid = true;
	} else {
		m_cache_valid = false;
	}
}
