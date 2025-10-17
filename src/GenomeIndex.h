/*
 * GenomeIndex.h
 *
 * Index for single-file genome storage
 * Maps contig ID to offset/length in genome.seq
 */

#ifndef GENOMEINDEX_H_
#define GENOMEINDEX_H_

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Index entry for a single contig
struct ContigIndexEntry {
    uint32_t chromid;      // Chromosome ID from GenomeChromKey
    uint64_t offset;       // Byte offset in genome.seq
    uint64_t length;       // Length in bases
    uint64_t reserved;     // Reserved for future use (compression flags, etc.)
    string   name;         // Contig name

    ContigIndexEntry() : chromid(0), offset(0), length(0), reserved(0) {}
    ContigIndexEntry(uint32_t id, uint64_t off, uint64_t len, const string &n)
        : chromid(id), offset(off), length(len), reserved(0), name(n) {}
};

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!
class GenomeIndex {
public:
    enum Errors {
        FILE_READ_FAILED,
        INVALID_FORMAT,
        VERSION_MISMATCH,
        CHECKSUM_FAILED
    };

    GenomeIndex();
    ~GenomeIndex();

    // Load index from file
    void load(const string &index_path);

    // Get index entry for a chromosome
    const ContigIndexEntry* get_entry(uint32_t chromid) const;

    // Check if index is loaded
    bool is_loaded() const { return m_loaded; }

    // Get all entries (for iteration)
    const vector<ContigIndexEntry>& get_all_entries() const { return m_entries; }

    // Get number of contigs
    size_t get_num_contigs() const { return m_entries.size(); }

private:
    bool m_loaded;
    vector<ContigIndexEntry> m_entries;
    unordered_map<uint32_t, size_t> m_chromid_to_index;

    // Magic header and version
    static const char MAGIC_HEADER[8];
    static const uint32_t INDEX_VERSION = 1;

    // Helper: validate magic header
    bool validate_header(FILE *fp);

    // Helper: compute checksum of index entries
    uint64_t compute_checksum(const vector<ContigIndexEntry> &entries);
};

#endif /* GENOMEINDEX_H_ */
