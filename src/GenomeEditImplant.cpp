/*
 * GenomeEditImplant.cpp
 *
 * C++ fast path for ggenome.implant: read a reference FASTA, apply
 * perturbations (donor sequence replacements), write output FASTA + .fai.
 *
 * Replaces the pure-R streaming path that was bottlenecked by readLines().
 */

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <map>
#include <string>
#include <vector>

#include "rdbutils.h"

using namespace std;
using namespace rdb;

namespace {

struct Perturbation {
    int start;   // 0-based
    int end;     // 0-based, exclusive
    const char *donor; // pointer into R string storage (valid for call duration)
};

struct FaiEntry {
    string name;
    long long length;   // sequence length in bases
    long long offset;   // byte offset of first base in output file
    int linebases;      // bases per line
    int linewidth;      // bytes per line (linebases + 1 for newline)
};

} // anonymous namespace

extern "C" {

SEXP C_ggenome_implant(SEXP _genome_fasta, SEXP _output,
                       SEXP _pert_chrom, SEXP _pert_start, SEXP _pert_end,
                       SEXP _pert_seq, SEXP _line_width)
{
    try {
        // --- unpack R arguments ---
        const char *genome_fasta = CHAR(STRING_ELT(_genome_fasta, 0));
        const char *output_path  = CHAR(STRING_ELT(_output, 0));
        int n_perts = Rf_length(_pert_chrom);
        int line_width = Rf_asInteger(_line_width);

        // --- build perturbation index: chrom -> sorted-desc-by-start ---
        map<string, vector<Perturbation>> pert_map;
        for (int i = 0; i < n_perts; i++) {
            Perturbation p;
            p.start = INTEGER(_pert_start)[i];
            p.end   = INTEGER(_pert_end)[i];
            p.donor = CHAR(STRING_ELT(_pert_seq, i));
            string chrom(CHAR(STRING_ELT(_pert_chrom, i)));
            pert_map[chrom].push_back(p);
        }
        // Sort each chrom's perturbations by start descending
        for (auto &kv : pert_map) {
            sort(kv.second.begin(), kv.second.end(),
                 [](const Perturbation &a, const Perturbation &b) {
                     return a.start > b.start;
                 });
        }

        // --- open files with large buffers ---
        FILE *fin = fopen(genome_fasta, "r");
        if (!fin) {
            Rf_error("Cannot open input FASTA: %s", genome_fasta);
        }

        FILE *fout = fopen(output_path, "w");
        if (!fout) {
            fclose(fin);
            Rf_error("Cannot open output file: %s", output_path);
        }

        // 4MB I/O buffers
        const size_t BUF_SIZE = 4 * 1024 * 1024;
        vector<char> inbuf(BUF_SIZE);
        vector<char> outbuf(BUF_SIZE);
        setvbuf(fin,  inbuf.data(),  _IOFBF, BUF_SIZE);
        setvbuf(fout, outbuf.data(), _IOFBF, BUF_SIZE);

        // --- streaming state ---
        vector<char> seq;        // accumulated sequence for current chrom
        seq.reserve(256 * 1024 * 1024); // 256MB pre-alloc for large chroms
        string current_chrom;
        vector<FaiEntry> fai_entries;
        long long byte_offset = 0;

        // Line buffer for fgets
        const size_t LINE_BUF = 65536;
        vector<char> line(LINE_BUF);

        // --- helper lambda: flush one chromosome ---
        auto flush_chrom = [&]() {
            if (current_chrom.empty()) return;

            long long seq_len = (long long)seq.size();

            // Apply perturbations (reverse-sorted by start, so coordinates stay valid)
            auto it = pert_map.find(current_chrom);
            if (it != pert_map.end()) {
                for (const auto &p : it->second) {
                    if (p.start < 0 || p.end > seq_len) {
                        fclose(fin);
                        fclose(fout);
                        Rf_error("Interval out of bounds: %s:%d-%d (chromosome length: %lld)",
                                 current_chrom.c_str(), p.start, p.end, seq_len);
                    }
                    int donor_len = p.end - p.start;
                    memcpy(&seq[p.start], p.donor, donor_len);
                }
            }

            // Write header
            fprintf(fout, ">%s\n", current_chrom.c_str());
            long long header_bytes = (long long)current_chrom.size() + 2; // '>' + name + '\n'

            // Write sequence lines
            long long pos = 0;
            while (pos < seq_len) {
                long long chunk = min((long long)line_width, seq_len - pos);
                fwrite(&seq[pos], 1, chunk, fout);
                fputc('\n', fout);
                pos += chunk;
            }

            // Collect .fai metadata
            FaiEntry entry;
            entry.name = current_chrom;
            entry.length = seq_len;
            entry.offset = byte_offset + header_bytes;
            entry.linebases = (seq_len > 0) ? min((long long)line_width, seq_len) : 0;
            entry.linewidth = (seq_len > 0) ? entry.linebases + 1 : 0;
            fai_entries.push_back(entry);

            // Update byte offset
            long long n_full_lines = seq_len / line_width;
            long long remainder = seq_len % line_width;
            long long data_bytes = n_full_lines * (line_width + 1);
            if (remainder > 0) data_bytes += remainder + 1;
            byte_offset += header_bytes + data_bytes;

            seq.clear();
        };

        // --- read input FASTA ---
        while (fgets(line.data(), LINE_BUF, fin)) {
            if (line[0] == '>') {
                // Flush previous chromosome
                flush_chrom();

                // Parse new chromosome name: strip '>' and trailing whitespace/newline
                char *name_start = line.data() + 1;
                // Find end: first whitespace or newline
                char *p = name_start;
                while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
                current_chrom.assign(name_start, p - name_start);
            } else {
                // Sequence line: append with inline toupper, strip newline
                for (char *p = line.data(); *p && *p != '\n' && *p != '\r'; p++) {
                    seq.push_back((char)toupper((unsigned char)*p));
                }
            }
        }
        // Flush last chromosome
        flush_chrom();

        fclose(fin);

        // --- write .fai file ---
        string fai_path = string(output_path) + ".fai";
        FILE *ffai = fopen(fai_path.c_str(), "w");
        if (!ffai) {
            fclose(fout);
            Rf_error("Cannot open .fai file for writing: %s", fai_path.c_str());
        }
        for (const auto &e : fai_entries) {
            fprintf(ffai, "%s\t%lld\t%lld\t%d\t%d\n",
                    e.name.c_str(), e.length, e.offset, e.linebases, e.linewidth);
        }
        fclose(ffai);
        fclose(fout);

        // --- return .fai as R data frame ---
        int n = (int)fai_entries.size();

        SEXP name_col   = PROTECT(Rf_allocVector(STRSXP,  n));
        SEXP length_col = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP offset_col = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP lbase_col  = PROTECT(Rf_allocVector(INTSXP,  n));
        SEXP lwidth_col = PROTECT(Rf_allocVector(INTSXP,  n));

        for (int i = 0; i < n; i++) {
            SET_STRING_ELT(name_col, i, Rf_mkChar(fai_entries[i].name.c_str()));
            REAL(length_col)[i] = (double)fai_entries[i].length;
            REAL(offset_col)[i] = (double)fai_entries[i].offset;
            INTEGER(lbase_col)[i]  = fai_entries[i].linebases;
            INTEGER(lwidth_col)[i] = fai_entries[i].linewidth;
        }

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
        SET_VECTOR_ELT(result, 0, name_col);
        SET_VECTOR_ELT(result, 1, length_col);
        SET_VECTOR_ELT(result, 2, offset_col);
        SET_VECTOR_ELT(result, 3, lbase_col);
        SET_VECTOR_ELT(result, 4, lwidth_col);

        SEXP col_names = PROTECT(Rf_allocVector(STRSXP, 5));
        SET_STRING_ELT(col_names, 0, Rf_mkChar("name"));
        SET_STRING_ELT(col_names, 1, Rf_mkChar("length"));
        SET_STRING_ELT(col_names, 2, Rf_mkChar("offset"));
        SET_STRING_ELT(col_names, 3, Rf_mkChar("linebases"));
        SET_STRING_ELT(col_names, 4, Rf_mkChar("linewidth"));
        Rf_setAttrib(result, R_NamesSymbol, col_names);

        SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(row_names)[0] = NA_INTEGER;
        INTEGER(row_names)[1] = -n;
        Rf_setAttrib(result, R_RowNamesSymbol, row_names);

        SEXP df_class = PROTECT(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(df_class, 0, Rf_mkChar("data.frame"));
        Rf_setAttrib(result, R_ClassSymbol, df_class);

        UNPROTECT(9);
        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
    return R_NilValue;
}

} // extern "C"
