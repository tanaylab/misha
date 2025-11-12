/*
 * GenomeSeqMultiImport.cpp
 *
 * Multi-FASTA import for indexed genome format
 * Creates genome.seq + genome.idx from multi-FASTA file
 */

#include <algorithm>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <string>

#include "rdbutils.h"
#include "BufferedFile.h"
#include "GenomeIndex.h"
#include "CRC64.h"

using namespace std;
using namespace rdb;

namespace {

// Sanitize FASTA header to extract clean contig name
string sanitize_fasta_header(const string &header) {
    // Remove leading '>' if present
    string clean = header;
    if (!clean.empty() && clean[0] == '>') {
        clean = clean.substr(1);
    }

    // Trim leading/trailing whitespace
    size_t start = 0;
    while (start < clean.size() && isspace(clean[start])) {
        start++;
    }
    size_t end = clean.size();
    while (end > start && isspace(clean[end - 1])) {
        end--;
    }
    clean = clean.substr(start, end - start);

    // Extract first token (before whitespace)
    size_t space_pos = clean.find_first_of(" \t");
    if (space_pos != string::npos) {
        clean = clean.substr(0, space_pos);
    }

    // For headers with pipes (e.g., "gi|12345|ref|NC_000001.1|"), extract the meaningful part
    // Strategy: keep the last non-empty segment before a pipe or the whole string
    size_t last_pipe = clean.find_last_of('|');
    if (last_pipe != string::npos && last_pipe + 1 < clean.size()) {
        // There's content after the last pipe - use it
        clean = clean.substr(last_pipe + 1);
    } else if (last_pipe == clean.size() - 1) {
        // Trailing pipe - remove it and try again
        clean = clean.substr(0, last_pipe);
        // Look for a segment with actual content (not just database identifiers)
        // Try to find "ref|something" or similar patterns
        size_t prev_pipe = clean.find_last_of('|');
        if (prev_pipe != string::npos && prev_pipe + 1 < clean.size()) {
            string last_segment = clean.substr(prev_pipe + 1);
            // Check if it looks like a meaningful identifier (not just numbers)
            bool has_alpha = false;
            for (char c : last_segment) {
                if (isalpha(c)) {
                    has_alpha = true;
                    break;
                }
            }
            if (has_alpha || last_segment.find('.') != string::npos) {
                clean = last_segment;
            }
        }
    }

    // Remove common prefixes that might remain
    const char *prefixes[] = {"lcl|", "gi|", "ref|", "gnl|", "gb|", "emb|", "dbj|", "pir|", "prf|", "sp|", "tr|", NULL};
    bool changed = true;
    while (changed) {
        changed = false;
        for (const char **prefix = prefixes; *prefix != NULL; prefix++) {
            size_t prefix_len = strlen(*prefix);
            if (clean.size() > prefix_len &&
                clean.compare(0, prefix_len, *prefix) == 0) {
                clean = clean.substr(prefix_len);
                changed = true;
                break;
            }
        }
    }

    // Replace problematic characters with underscore
    for (size_t i = 0; i < clean.size(); i++) {
        char c = clean[i];
        if (!isalnum(c) && c != '_' && c != '-' && c != '.') {
            clean[i] = '_';
        }
    }

    // Ensure name is not empty
    if (clean.empty()) {
        clean = "contig";
    }

    return clean;
}

// Compute checksum for index entries (same algorithm as GenomeIndex)
uint64_t compute_index_checksum(const vector<ContigIndexEntry> &entries) {
    misha::CRC64 crc64;
    uint64_t checksum = crc64.init_incremental();

    for (const auto &entry : entries) {
        // Hash chromid
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.chromid, sizeof(entry.chromid));

        // Hash name
        if (!entry.name.empty()) {
            checksum = crc64.compute_incremental(checksum,
                (const unsigned char*)entry.name.c_str(), entry.name.size());
        }

        // Hash offset and length
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.offset, sizeof(entry.offset));
        checksum = crc64.compute_incremental(checksum,
            (const unsigned char*)&entry.length, sizeof(entry.length));
    }

    return crc64.finalize_incremental(checksum);
}

// Write genome.idx file
void write_index_file(const string &index_path,
                      const vector<ContigIndexEntry> &entries) {
    FILE *fp = fopen(index_path.c_str(), "wb");
    if (!fp) {
        verror("Failed to create index file %s: %s",
               index_path.c_str(), strerror(errno));
    }

    // Write magic header "MISHAIDX"
    const char magic_header[8] = {'M','I','S','H','A','I','D','X'};
    if (fwrite(magic_header, 1, 8, fp) != 8) {
        fclose(fp);
        verror("Failed to write index header to %s", index_path.c_str());
    }

    // Write version
    uint32_t version = 1;  // INDEX_VERSION
    if (fwrite(&version, sizeof(version), 1, fp) != 1) {
        fclose(fp);
        verror("Failed to write version to %s", index_path.c_str());
    }

    // Write number of contigs
    uint32_t num_contigs = entries.size();
    if (fwrite(&num_contigs, sizeof(num_contigs), 1, fp) != 1) {
        fclose(fp);
        verror("Failed to write contig count to %s", index_path.c_str());
    }

    // Compute and write checksum
    uint64_t checksum = compute_index_checksum(entries);
    if (fwrite(&checksum, sizeof(checksum), 1, fp) != 1) {
        fclose(fp);
        verror("Failed to write checksum to %s", index_path.c_str());
    }

    // Write contig entries
    for (const auto &entry : entries) {
        // Write chromid
        if (fwrite(&entry.chromid, sizeof(entry.chromid), 1, fp) != 1) {
            fclose(fp);
            verror("Failed to write chromid to %s", index_path.c_str());
        }

        // Write name length and name
        uint16_t name_length = entry.name.size();
        if (fwrite(&name_length, sizeof(name_length), 1, fp) != 1) {
            fclose(fp);
            verror("Failed to write name length to %s", index_path.c_str());
        }

        if (name_length > 0) {
            if (fwrite(entry.name.c_str(), 1, name_length, fp) != name_length) {
                fclose(fp);
                verror("Failed to write contig name to %s", index_path.c_str());
            }
        }

        // Write offset, length, reserved
        if (fwrite(&entry.offset, sizeof(entry.offset), 1, fp) != 1 ||
            fwrite(&entry.length, sizeof(entry.length), 1, fp) != 1 ||
            fwrite(&entry.reserved, sizeof(entry.reserved), 1, fp) != 1) {
            fclose(fp);
            verror("Failed to write entry data to %s", index_path.c_str());
        }
    }

    fclose(fp);
}

} // anonymous namespace

extern "C" {

// Import multi-FASTA file to indexed format
// Arguments:
//   _fasta: path to multi-FASTA file
//   _seq: path to output genome.seq file
//   _index: path to output genome.idx file
//   _sort: logical, whether to sort chromosomes alphabetically (default TRUE for backward compatibility)
//   _envir: environment
// Returns:
//   Data frame with columns: name (character), size (numeric)
SEXP gseq_multifasta_import(SEXP _fasta, SEXP _seq, SEXP _index, SEXP _sort, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;

        // Validate arguments
        if (!Rf_isString(_fasta) || Rf_length(_fasta) != 1)
            verror("fasta argument is not a string");

        if (!Rf_isString(_seq) || Rf_length(_seq) != 1)
            verror("seq argument is not a string");

        if (!Rf_isString(_index) || Rf_length(_index) != 1)
            verror("index argument is not a string");

        // Parse sort argument (default to TRUE for backward compatibility)
        bool sort_chromosomes = true;
        if (!Rf_isNull(_sort)) {
            if (!Rf_isLogical(_sort) || Rf_length(_sort) != 1)
                verror("sort argument must be a logical value");
            sort_chromosomes = Rf_asLogical(_sort);
        }

        const char *fasta_fname = CHAR(STRING_ELT(_fasta, 0));
        const char *seq_fname = CHAR(STRING_ELT(_seq, 0));
        const char *index_fname = CHAR(STRING_ELT(_index, 0));

        // Open files
        BufferedFile fasta_file;
        BufferedFile seq_file;

        if (fasta_file.open(fasta_fname, "r"))
            verror("Failed to open FASTA file %s: %s",
                   fasta_fname, strerror(errno));

        if (seq_file.open(seq_fname, "w"))
            verror("Failed to create sequence file %s: %s",
                   seq_fname, strerror(errno));

        // Parse multi-FASTA and build index
        vector<ContigIndexEntry> entries;
        vector<char> seq_buffer;
        string current_header;
        uint64_t current_offset = 0;
        int contig_index = -1;
        bool in_sequence = false;

        // Read line by line
        string line;
        vector<char> line_buf;
        int c;

        while ((c = fasta_file.getc()) != EOF) {
            if (c == '\n') {
                // Process accumulated line
                line.assign(line_buf.begin(), line_buf.end());
                line_buf.clear();

                if (line.empty()) {
                    continue;
                }

                if (line[0] == '>') {
                    // New contig header

                    // Finalize previous contig if any
                    if (in_sequence && contig_index >= 0) {
                        // Flush sequence buffer
                        if (seq_buffer.size() > 0) {
                            seq_file.write(&seq_buffer.front(), seq_buffer.size());
                            if (seq_file.error()) {
                                verror("Failed to write sequence data: %s",
                                       strerror(errno));
                            }

                            // Accumulate final chunk (earlier chunks already accumulated during periodic flushes)
                            uint64_t contig_length = seq_buffer.size();
                            entries.back().length += contig_length;  // Use += not = to preserve earlier chunks

                            current_offset += contig_length;
                            seq_buffer.clear();
                        }
                    }

                    // Start new contig
                    contig_index++;
                    string contig_name = sanitize_fasta_header(line);

                    ContigIndexEntry entry;
                    entry.chromid = contig_index;  // Temporary chromid (will be reassigned after sorting)
                    entry.offset = current_offset;
                    entry.length = 0;  // Will be updated when contig is complete
                    entry.reserved = 0;
                    entry.name = contig_name;

                    entries.push_back(entry);
                    in_sequence = true;

                } else if (line[0] == ';') {
                    // Comment line, skip
                    continue;

                } else if (in_sequence) {
                    // Sequence data
                    for (char ch : line) {
                        if (isalpha(ch) || ch == '-') {
                            seq_buffer.push_back(ch);

                            // Flush buffer periodically to avoid excessive memory
                            if (seq_buffer.size() >= 1048576) {  // 1 MB
                                seq_file.write(&seq_buffer.front(), seq_buffer.size());
                                if (seq_file.error()) {
                                    verror("Failed to write sequence data: %s",
                                           strerror(errno));
                                }
                                entries.back().length += seq_buffer.size();
                                current_offset += seq_buffer.size();
                                seq_buffer.clear();
                            }
                        } else if (!isspace(ch)) {
                            verror("Invalid character '%c' in FASTA sequence", ch);
                        }
                    }
                }
            } else {
                line_buf.push_back(c);
            }
        }

        // Process final line if file doesn't end with newline
        if (!line_buf.empty()) {
            line.assign(line_buf.begin(), line_buf.end());
            line_buf.clear();

            if (!line.empty()) {
                if (line[0] == '>') {
                    // Final line is a header (contig with no sequence)
                    // Finalize previous contig if any
                    if (in_sequence && contig_index >= 0 && seq_buffer.size() > 0) {
                        seq_file.write(&seq_buffer.front(), seq_buffer.size());
                        if (seq_file.error()) {
                            verror("Failed to write sequence data: %s", strerror(errno));
                        }
                        entries.back().length += seq_buffer.size();
                        current_offset += seq_buffer.size();
                        seq_buffer.clear();
                    }

                    // Start new contig (will have zero length)
                    contig_index++;
                    string contig_name = sanitize_fasta_header(line);

                    ContigIndexEntry entry;
                    entry.chromid = contig_index;
                    entry.offset = current_offset;
                    entry.length = 0;
                    entry.reserved = 0;
                    entry.name = contig_name;

                    entries.push_back(entry);
                    in_sequence = true;

                } else if (line[0] != ';' && in_sequence) {
                    // Final line is sequence data
                    for (char ch : line) {
                        if (isalpha(ch) || ch == '-') {
                            seq_buffer.push_back(ch);
                        } else if (!isspace(ch)) {
                            verror("Invalid character '%c' in FASTA sequence", ch);
                        }
                    }
                }
            }
        }

        // Check for read errors
        if (fasta_file.error()) {
            verror("Error reading FASTA file %s: %s",
                   fasta_fname, strerror(errno));
        }

        // Finalize last contig
        if (in_sequence && contig_index >= 0 && seq_buffer.size() > 0) {
            seq_file.write(&seq_buffer.front(), seq_buffer.size());
            if (seq_file.error()) {
                verror("Failed to write sequence data: %s", strerror(errno));
            }
            entries.back().length += seq_buffer.size();
            seq_buffer.clear();
        }

        // Close sequence file
        seq_file.close();

        // Check if we got any contigs
        if (entries.empty()) {
            verror("No contigs found in FASTA file %s", fasta_fname);
        }

        // Check for duplicate sanitized names
        std::unordered_map<string, int> name_counts;
        for (const auto &entry : entries) {
            name_counts[entry.name]++;
        }
        for (const auto &pair : name_counts) {
            if (pair.second > 1) {
                verror("Duplicate contig name '%s' after sanitization (%d occurrences). "
                       "Please ensure FASTA headers produce unique names.",
                       pair.first.c_str(), pair.second);
            }
        }

        // Optionally sort entries alphabetically by name
        // This maintains backward compatibility (default: sort=TRUE)
        // For database conversion, we want to preserve FASTA order (sort=FALSE)
        if (sort_chromosomes) {
            // Sort entries by name and reassign chromids to match sorted order
            std::sort(entries.begin(), entries.end(),
                      [](const ContigIndexEntry &a, const ContigIndexEntry &b) {
                          return a.name < b.name;
                      });

            // Reassign chromids to match sorted order (0, 1, 2, ...)
            for (size_t i = 0; i < entries.size(); i++) {
                entries[i].chromid = i;
            }
        }
        // If not sorting, chromids are already assigned in FASTA order (0, 1, 2, ...)

        // Write index file
        write_index_file(index_fname, entries);

        // Return data frame with contig names and sizes
        SEXP name_col = PROTECT(Rf_allocVector(STRSXP, entries.size()));
        SEXP size_col = PROTECT(Rf_allocVector(REALSXP, entries.size()));

        for (size_t i = 0; i < entries.size(); i++) {
            SET_STRING_ELT(name_col, i, Rf_mkChar(entries[i].name.c_str()));
            REAL(size_col)[i] = (double)entries[i].length;
        }

        // Create data frame
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(result, 0, name_col);
        SET_VECTOR_ELT(result, 1, size_col);

        // Set column names
        SEXP col_names = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(col_names, 0, Rf_mkChar("name"));
        SET_STRING_ELT(col_names, 1, Rf_mkChar("size"));
        Rf_setAttrib(result, R_NamesSymbol, col_names);

        // Set class to data.frame
        SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(row_names)[0] = NA_INTEGER;
        INTEGER(row_names)[1] = -(int)entries.size();
        Rf_setAttrib(result, R_RowNamesSymbol, row_names);

        SEXP df_class = PROTECT(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(df_class, 0, Rf_mkChar("data.frame"));
        Rf_setAttrib(result, R_ClassSymbol, df_class);

        UNPROTECT(6);

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
