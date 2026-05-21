#include <cstdint>
#include <cstring>
#include <stdio.h>

#include <memory>
#include <unordered_map>
#include <string>
#include <vector>

#include <sys/wait.h>
#include <zlib.h>

#include "rdbinterval.h"
#include "rdbutils.h"
#include "rdbprogress.h"

#include "BufferedFile.h"
#include "HashFunc.h"

#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "TrackExpressionScanner.h"

using namespace std;
using namespace rdb;

namespace {

// ----- ByteSource: byte-stream abstraction over the import file ------------
//
// One of three concrete sources is selected by open_source() based on a
// magic-byte sniff. The FSM in gtrackimport_mappedseq reads via src->getc()
// and checks src->error() at EOF; PipeSource additionally exposes the child
// process exit status via finish() after EOF.

class ByteSource {
public:
    virtual ~ByteSource() = default;
    virtual int getc() = 0;
    virtual bool error() const = 0;
};

class PlainSource : public ByteSource {
    BufferedFile buf_;
public:
    explicit PlainSource(const std::string &path) {
        buf_.open(path.c_str(), "r");
        if (buf_.error())
            TGLError("Failed to open %s: %s", path.c_str(), strerror(errno));
    }
    int getc() override { return buf_.getc(); }
    bool error() const override { return buf_.error() != 0; }
};

class GzipSource : public ByteSource {
    gzFile gz_;
public:
    explicit GzipSource(const std::string &path) : gz_(nullptr) {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_)
            TGLError("Failed to open gzipped %s: %s", path.c_str(), strerror(errno));
    }
    ~GzipSource() override { if (gz_) gzclose(gz_); }
    int getc() override { return gzgetc(gz_); }
    bool error() const override {
        if (!gz_) return true;
        int err = 0;
        gzerror(gz_, &err);
        return err != Z_OK && err != Z_STREAM_END;
    }
};

// PipeSource wraps a popen() child (used for BAM auto-detect via
// `samtools view <path>`). finish() must be called after the parser has
// drained the pipe so the caller can inspect the child's exit status; if
// finish() is not called, the destructor pcloses as a fallback (discarding
// the status). 127 from pclose means "command not found" (samtools missing).
class PipeSource : public ByteSource {
    FILE *fp_;
    std::string cmd_;
    bool finished_;
    int finish_status_;
public:
    explicit PipeSource(const std::string &cmd)
        : fp_(nullptr), cmd_(cmd), finished_(false), finish_status_(-1)
    {
        fp_ = popen(cmd.c_str(), "r");
        if (!fp_)
            TGLError("Failed to popen '%s': %s", cmd.c_str(), strerror(errno));
    }
    ~PipeSource() override {
        if (fp_) pclose(fp_);
    }
    int getc() override { return ::getc(fp_); }
    bool error() const override { return fp_ ? ferror(fp_) != 0 : true; }
    int finish() {
        if (finished_) return finish_status_;
        finished_ = true;
        finish_status_ = pclose(fp_);
        fp_ = nullptr;
        return finish_status_;
    }
    const std::string &cmd() const { return cmd_; }
};

// Single-quote escape for paths handed to /bin/sh via popen: every ' becomes
// '\'' and the whole string is wrapped in '...'. Safe against shell
// metacharacters in filesystem paths.
static std::string shellquote_single(const std::string &s) {
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
        if (c == '\'') out.append("'\\''");
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

// Peek up to nbuf magic bytes from path. Returns count read (0 on open
// failure).
static size_t read_magic_bytes(const std::string &path, unsigned char *buf, size_t nbuf) {
    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) return 0;
    size_t n = fread(buf, 1, nbuf, fp);
    fclose(fp);
    return n;
}

// Select the right source for a path:
//   - BAM (bgzip magic 1f 8b 08 04): popen("samtools view <quoted-path>")
//   - plain gzip (magic 1f 8b ...):  GzipSource via zlib
//   - else:                          PlainSource via BufferedFile
static std::unique_ptr<ByteSource> open_source(const std::string &path) {
    unsigned char magic[4] = {0, 0, 0, 0};
    size_t n = read_magic_bytes(path, magic, 4);
    if (n == 4 && magic[0] == 0x1f && magic[1] == 0x8b &&
                  magic[2] == 0x08 && magic[3] == 0x04) {
        // bgzip / BAM
        const std::string cmd = "samtools view " + shellquote_single(path);
        return std::unique_ptr<ByteSource>(new PipeSource(cmd));
    }
    if (n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
        // plain gzip
        return std::unique_ptr<ByteSource>(new GzipSource(path));
    }
    return std::unique_ptr<ByteSource>(new PlainSource(path));
}

}  // namespace

extern "C" {

SEXP gtrackimport_mappedseq(SEXP _track, SEXP _infile, SEXP _pileup, SEXP _binsize, SEXP _cols_order, SEXP _remove_dups, SEXP _envir)
{
	enum { SEQ_COL, CHROM_COL, COORD_COL, STRAND_COL, NUM_COLS };
	const char *COL_NAMES[NUM_COLS] = { "sequence", "chromosome", "coordinate", "strand" };

	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_infile) || Rf_length(_infile) != 1)
			verror("File argument is not a string");

		if (Rf_length(_pileup) != 1 || ((!Rf_isReal(_pileup) || REAL(_pileup)[0] != (int)REAL(_pileup)[0]) && !Rf_isInteger(_pileup)))
			verror("Pileup argument is not an integer");

		if (Rf_length(_binsize) != 1 || ((!Rf_isReal(_binsize) || REAL(_binsize)[0] != (int)REAL(_binsize)[0]) && !Rf_isInteger(_binsize)))
			verror("Binsize argument is not an integer");

		if (!Rf_isNull(_cols_order) && (Rf_length(_cols_order) != NUM_COLS || (!Rf_isReal(_cols_order) && !Rf_isInteger(_cols_order))))
			verror("cols.order argument must be a vector with %d numeric values", NUM_COLS);

		if (!Rf_isNull(_cols_order) && Rf_isReal(_cols_order)) {
			for (int i = 0; i < NUM_COLS; i++) {
				if (REAL(_cols_order)[i] != (int)REAL(_cols_order)[i])
					verror("cols.order is not an integer");
			}
		}

		if (Rf_length(_remove_dups) > 1 || !Rf_isLogical(_remove_dups))
			verror("remove.dups argument must be a logical value");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *infilename = CHAR(STRING_ELT(_infile, 0));
		int pileup = Rf_isReal(_pileup) ? (int)REAL(_pileup)[0] : INTEGER(_pileup)[0];
		double binsize = Rf_isReal(_binsize) ? (int)REAL(_binsize)[0] : INTEGER(_binsize)[0];
		int cols_order[NUM_COLS];
		bool remove_dups = LOGICAL(_remove_dups)[0];
		bool is_sam_format = Rf_isNull(_cols_order);

		if (is_sam_format) { // SAM format
			cols_order[SEQ_COL] = 10;
			cols_order[CHROM_COL] = 3;
			cols_order[COORD_COL] = 4;
			cols_order[STRAND_COL] = 2;
		} else {
			for (int i = 0; i < NUM_COLS; i++)
				cols_order[i] = Rf_isReal(_cols_order) ? (int)REAL(_cols_order)[i] : INTEGER(_cols_order)[i];
		}

		if (pileup < 0)
			verror("Pileup cannot be negative");

		if (pileup == 0 && binsize >= 0)
			verror("Invalid binsize.\nSparse track is created when pileup is zero. Binsize must be set to -1 then.");

		if (pileup > 0 && binsize <= 0)
			verror("Invalid binsize.\nDense track is created when pileup is greater than zero. Binsize must be a positive integer then.");

		for (int i = 0; i < NUM_COLS; i++) {
			if (cols_order[i] <= 0)
				verror("Invalid columns order: %s column's order is %d", COL_NAMES[i]);

			for (int j = i + 1; j < NUM_COLS; j++) {
				if (cols_order[i] == cols_order[j])
					verror("Invalid columns order: %s column has the same order as %s column", COL_NAMES[i], COL_NAMES[j]);
			}
		}

		IntervUtils iu(_envir);
		GIntervals all_genome_intervs;
		unordered_map<string, int> str2chrom;
		int64_t genome_len = 0;

		iu.get_all_genome_intervs(all_genome_intervs);
		for (GIntervals::const_iterator iinterv = all_genome_intervs.begin(); iinterv != all_genome_intervs.end(); ++iinterv) {
			str2chrom[iu.id2chrom(iinterv->chromid)] = iinterv - all_genome_intervs.begin();
			genome_len += iinterv->end;
		}

		unsigned num_chroms = all_genome_intervs.size();
		vector<unsigned> num_mapped(num_chroms, 0);
		vector<unsigned> num_dups(num_chroms, 0);
		vector< vector<int64_t> > coords(2 * num_chroms);

		string dirname = create_track_dir(_envir, track);
		int total_unmapped = 0;
		std::unique_ptr<ByteSource> src = open_source(infilename);

		int col = 1;
		int active_col_idx = -1;
		int c;
		string str[NUM_COLS];
		// int line = 1;
		int pos = 0;

		Progress_reporter progress;
		progress.init(genome_len, 1000000);

		for (int i = 0; i < NUM_COLS; i++) {
			if (cols_order[i] == 1) {
				active_col_idx = i;
				break;
			}
		}

		while (1) {
			c = src->getc();

			// in SAM file skip the lines that start with @
			if (!pos && is_sam_format && c == '@') {
				while (1) {
					c = src->getc();
					if (c == '\n' || c == EOF)
						break;
				}

				if (c == EOF) 
					break;
				// ++line;
				continue;
			}
			++pos;

			if (c == '\n' || c == EOF || c == '\t') {
				if (c == '\n' || c == EOF) {
					int num_nonempty_strs = 0;
					bool mapped = false;

					pos = 0;
					for (int i = 0; i < NUM_COLS; i++) {
						if (!str[i].empty())
							num_nonempty_strs++;
					}

					while (num_nonempty_strs == NUM_COLS) {
						unordered_map<string, int>::iterator istr2chrom;
						int chrom_idx;
						int64_t coord;
						char *endptr;

						if ((istr2chrom = str2chrom.find(str[CHROM_COL])) == str2chrom.end())
							break;
						chrom_idx = istr2chrom->second;

						coord = strtoll(str[COORD_COL].c_str(), &endptr, 10);
						if (*endptr || coord < 0 || coord >= all_genome_intervs[chrom_idx].end)
							break;

						if (is_sam_format) {
							uint64_t num = strtoll(str[STRAND_COL].c_str(), &endptr, 0);
							if (*endptr)
								break;
							str[STRAND_COL] = num & 0x10 ? "-" : "+";
						}

						if (str[STRAND_COL] == "+" || str[STRAND_COL] == "F")
							coords[chrom_idx].push_back(coord);
						else if (str[STRAND_COL] == "-" || str[STRAND_COL] == "R")
							coords[num_chroms + chrom_idx].push_back(coord + str[SEQ_COL].size());
						else
							break;

						mapped = true;
						++num_mapped[chrom_idx];
						break;
					}

					if (!mapped && num_nonempty_strs)
						total_unmapped++;

					if (c == EOF)
						break;

					if (num_nonempty_strs > 0) {
						for (int i = 0; i < NUM_COLS; i++)
							str[i].clear();
					}
					col = 1;
					// line++;
				} else
					col++;

				active_col_idx = -1;
				for (int i = 0; i < NUM_COLS; i++) {
					if (cols_order[i] == col) {
						active_col_idx = i;
						break;
					}
				}
			} else if (active_col_idx >= 0)
				str[active_col_idx].push_back(c);

			check_interrupt();
			progress.report(1);
		}

		if (src->error())
			verror("Error while reading file %s: %s", infilename, strerror(errno));

		if (auto *pipe_src = dynamic_cast<PipeSource *>(src.get())) {
			int raw_status = pipe_src->finish();
			if (raw_status != 0) {
				int code = WIFEXITED(raw_status) ? WEXITSTATUS(raw_status) : -1;
				if (code == 127) {
					verror("BAM input detected at %s but samtools is not on PATH. "
					       "Install samtools (e.g. `apt-get install samtools` or "
					       "`conda install -c bioconda samtools`) or pre-convert: "
					       "`samtools view %s > %s.sam`.",
					       infilename, infilename, infilename);
				} else {
					verror("samtools view %s exited with code %d (raw status %d). "
					       "Run `samtools view %s | head` to see the underlying error.",
					       infilename, code, raw_status, infilename);
				}
			}
		}

		for (unsigned ichrom = 0; ichrom < num_chroms; ichrom++) {
			char filename[FILENAME_MAX];
			snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), iu.id2chrom(all_genome_intervs[ichrom].chromid).c_str());

			// dense track
			if (pileup) {
				GenomeTrackFixedBin gtrack;
				gtrack.init_write(filename, (unsigned)binsize, all_genome_intervs[ichrom].chromid);

				vector<float> trackvals((uint64_t)ceil(all_genome_intervs[ichrom].end / binsize), 0);

				for (int strand = 0; strand < 2; strand++) {
					vector<int64_t> &cur_coords = coords[strand * num_chroms + ichrom];
					sort(cur_coords.begin(), cur_coords.end());

					for (vector<int64_t>::const_iterator icoord = cur_coords.begin(); icoord != cur_coords.end(); ++icoord) {
						if (remove_dups && icoord != cur_coords.begin() && *icoord == *(icoord - 1)) {
							++num_dups[ichrom];
							continue;
						}

						int64_t from_coord = max(strand ? *icoord - pileup : *icoord, (int64_t)0);
						int64_t to_coord = min(strand ? *icoord : *icoord + pileup, all_genome_intervs[ichrom].end);
						int64_t from_bin = (int64_t)(from_coord / binsize);
						int64_t to_bin = (int64_t)ceil(to_coord / binsize) - 1;

						// If from/to bin equals to the last bin in the chromosome, then we should replace use the length of the tail rather than binsize;
						// yet we don't want to introduce another "if" statement + complications. So let the last bin be inaccurate.
						if (from_bin >= to_bin)
							trackvals[from_bin] += (to_coord - from_coord) / binsize;
						else {
							trackvals[from_bin] += from_bin + 1 - from_coord / binsize;
							trackvals[to_bin] += to_coord / binsize - to_bin;
							for (int64_t bin = from_bin + 1; bin < to_bin; ++bin)
								trackvals[bin]++;
						}
					}
				}

				for (vector<float>::const_iterator itrackval = trackvals.begin(); itrackval != trackvals.end(); ++itrackval) {
					gtrack.write_next_bin(*itrackval);
					check_interrupt();
				}
			}
			// sparse tracks
			else {
				GenomeTrackSparse gtrack;
				gtrack.init_write(filename, all_genome_intervs[ichrom].chromid);

				vector<int64_t> *pcoords[2] = { &coords[ichrom], &coords[num_chroms + ichrom] };
				vector<int64_t>::const_iterator icoords[2] = { pcoords[0]->begin(), pcoords[1]->begin() };

				sort(pcoords[0]->begin(), pcoords[0]->end());
				sort(pcoords[1]->begin(), pcoords[1]->end());

				while (icoords[0] != pcoords[0]->end() || icoords[1] != pcoords[1]->end()) {
					float val = 0;
					long coord = -1;

					if (icoords[0] != pcoords[0]->end() && (icoords[1] == pcoords[1]->end() || *icoords[1] >= *icoords[0])) {
						val = max(val + !remove_dups, (float)1);
						coord = *icoords[0];
						for (++icoords[0]; icoords[0] != pcoords[0]->end() && *icoords[0] == coord; ++icoords[0]) {
							++num_dups[ichrom];
							val += !remove_dups;
						}
					}

					if (icoords[1] != pcoords[1]->end() && (coord == -1 || *icoords[1] == coord)) {
						val = max(val + !remove_dups, (float)1);
						coord = *icoords[1];
						for (++icoords[1]; icoords[1] != pcoords[1]->end() && *icoords[1] == coord; ++icoords[1]) {
							++num_dups[ichrom];
							val += !remove_dups;
						}
					}

					check_interrupt();

					gtrack.write_next_interval(GInterval(all_genome_intervs[ichrom].chromid, coord, coord + 1, 0), val);
				}
			}

			progress.report(all_genome_intervs[ichrom].end);
		}

		progress.report_last();

		SEXP answer;
		SEXP chrom_stat, total_stat;
		SEXP chroms, chroms_idx, mapped, dups;
		SEXP col_names;
		SEXP row_names;

		int64_t total_mapped = 0;
		int64_t total_dups = 0;

		rprotect(chrom_stat = RSaneAllocVector(VECSXP, 3));
        rprotect(chroms_idx = RSaneAllocVector(INTSXP, num_chroms));
        rprotect(mapped = RSaneAllocVector(REALSXP, num_chroms));
        rprotect(dups = RSaneAllocVector(REALSXP, num_chroms));
        rprotect(chroms = RSaneAllocVector(STRSXP, num_chroms));
        rprotect(col_names = RSaneAllocVector(STRSXP, 3));
        rprotect(row_names = RSaneAllocVector(INTSXP, num_chroms));

		for (unsigned i = 0; i < num_chroms; i++) {
			INTEGER(chroms_idx)[i] = all_genome_intervs[i].chromid + 1;
			SET_STRING_ELT(chroms, i, Rf_mkChar(iu.id2chrom(i).c_str()));
			REAL(mapped)[i] = num_mapped[i];
			REAL(dups)[i] = num_dups[i];
			INTEGER(row_names)[i] = i + 1;

			total_mapped += num_mapped[i];
			total_dups += num_dups[i];
		}

		SET_STRING_ELT(col_names, 0, Rf_mkChar("chrom"));
		SET_STRING_ELT(col_names, 1, Rf_mkChar("mapped"));
		SET_STRING_ELT(col_names, 2, Rf_mkChar("dups"));

        Rf_setAttrib(chroms_idx, R_LevelsSymbol, chroms);
        Rf_setAttrib(chroms_idx, R_ClassSymbol, Rf_mkString("factor"));

        SET_VECTOR_ELT(chrom_stat, 0, chroms_idx);
        SET_VECTOR_ELT(chrom_stat, 1, mapped);
        SET_VECTOR_ELT(chrom_stat, 2, dups);

        Rf_setAttrib(chrom_stat, R_NamesSymbol, col_names);
        Rf_setAttrib(chrom_stat, R_ClassSymbol, Rf_mkString("data.frame"));
        Rf_setAttrib(chrom_stat, R_RowNamesSymbol, row_names);

		rprotect(total_stat = RSaneAllocVector(REALSXP, 4));
		REAL(total_stat)[0] = total_mapped + total_unmapped + total_dups;
		REAL(total_stat)[1] = total_mapped;
		REAL(total_stat)[2] = total_unmapped;
		REAL(total_stat)[3] = total_dups;
		Rf_setAttrib(total_stat, R_NamesSymbol, RSaneAllocVector(STRSXP, 3));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 0, Rf_mkChar("total"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 1, Rf_mkChar("total.mapped"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 2, Rf_mkChar("total.unmapped"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 3, Rf_mkChar("total.dups"));

		rprotect(answer = RSaneAllocVector(VECSXP, 2));

		SET_VECTOR_ELT(answer, 0, total_stat);
		SET_VECTOR_ELT(answer, 1, chrom_stat);

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
