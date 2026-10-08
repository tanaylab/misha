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

// ByteSource is a byte-stream abstraction over the import file. open_source()
// returns a PlainSource (default), GzipSource (1f 8b magic), or PipeSource
// (bgzip / BAM magic -> samtools view subprocess). The FSM reads via
// src->getc() and checks src->error() at EOF; error() returns false at clean
// EOF for all three sources. After EOF the FSM calls check_finalize() to give
// PipeSource a chance to surface its child's exit status.

class ByteSource {
public:
	virtual ~ByteSource() = default;
	virtual int getc() = 0;
	virtual bool error() const = 0;
	virtual void check_finalize(const char *infilename) { (void)infilename; }
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
// `samtools view <path>`). The caller invokes check_finalize() once parsing
// is done so a non-zero child exit (notably 127 == "command not found"
// for missing samtools) is surfaced as a clear error.
class PipeSource : public ByteSource {
	FILE *fp_;
	bool finished_;
public:
	explicit PipeSource(const std::string &cmd) : fp_(nullptr), finished_(false) {
		fp_ = popen(cmd.c_str(), "r");
		if (!fp_)
			TGLError("Failed to popen '%s': %s", cmd.c_str(), strerror(errno));
	}
	~PipeSource() override {
		// Dtors are implicitly noexcept; any throw during stack unwinding
		// would call std::terminate. The pclose itself can't throw, and
		// REprintf is well-behaved on R 4.x, but be defensive.
		try {
			if (!fp_) return;
			// We get here only on the throw path - the success path goes
			// through check_finalize() which already cleared fp_. Print
			// the child's exit status so a samtools crash mid-stream
			// still leaves a breadcrumb after the FSM's parser error.
			int raw_status = pclose(fp_);
			if (raw_status != 0) {
				int code = WIFEXITED(raw_status) ? WEXITSTATUS(raw_status) : -1;
				REprintf("note: samtools view exited with code %d (raw status %d) "
				         "after the import error above.\n", code, raw_status);
			}
		} catch (...) {
			// swallow - destructors must not propagate exceptions
		}
	}
	// getc_unlocked avoids the per-byte stdio mutex; PipeSource isn't shared.
	int getc() override { return getc_unlocked(fp_); }
	bool error() const override { return fp_ ? ferror(fp_) != 0 : true; }
	void check_finalize(const char *infilename) override {
		if (finished_) return;
		finished_ = true;
		int raw_status = pclose(fp_);
		fp_ = nullptr;
		if (raw_status == 0) return;
		int code = WIFEXITED(raw_status) ? WEXITSTATUS(raw_status) : -1;
		if (code == 127) {
			TGLError("BAM input detected at %s but samtools is not on PATH. "
			         "Install samtools (e.g. `apt-get install samtools` or "
			         "`conda install -c bioconda samtools`) or pre-convert: "
			         "`samtools view %s > %s.sam`.",
			         infilename, infilename, infilename);
		} else {
			TGLError("samtools view %s exited with code %d (raw status %d). "
			         "Run `samtools view %s | head` to see the underlying error.",
			         infilename, code, raw_status, infilename);
		}
	}
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

// Peek up to nbuf bytes from path. Returns count read (0 on open failure).
// We deliberately use raw fopen/fread rather than BufferedFile because the
// latter allocates a 128KB buffer for a 4-byte sniff.
static size_t read_magic_bytes(const std::string &path, unsigned char *buf, size_t nbuf) {
	FILE *fp = fopen(path.c_str(), "rb");
	if (!fp) return 0;
	size_t n = fread(buf, 1, nbuf, fp);
	fclose(fp);
	return n;
}

// True if the decompressed stream starts with the BAM magic "BAM\1".
static bool has_bam_payload(const std::string &path) {
	gzFile gz = gzopen(path.c_str(), "rb");
	if (!gz) return false;
	char magic[4];
	int n = gzread(gz, magic, 4);
	gzclose(gz);
	return n == 4 && !memcmp(magic, "BAM\1", 4);
}

static std::unique_ptr<ByteSource> open_source(const std::string &path) {
	unsigned char magic[4] = {0, 0, 0, 0};
	size_t n = read_magic_bytes(path, magic, 4);
	// bgzip uses gzip magic with FLG.FEXTRA (0x04) and method 0x08; that
	// distinguishes BAM/bgzip from plain gzip (which has FLG 0x00 or 0x08).
	// A bgzipped text file (e.g. a 10x fragments.tsv.gz) is read with zlib.
	if (n == 4 && magic[0] == 0x1f && magic[1] == 0x8b &&
	              magic[2] == 0x08 && magic[3] == 0x04 && has_bam_payload(path)) {
		return std::make_unique<PipeSource>("samtools view " + shellquote_single(path));
	}
	if (n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b)
		return std::make_unique<GzipSource>(path);
	return std::make_unique<PlainSource>(path);
}

// Reference length of a CIGAR (M, D, N, = and X operations); -1 for "*", a malformed CIGAR or an absurd length.
static int64_t cigar_ref_span(const std::string &cigar) {
	const int64_t MAX_LEN = (int64_t)1 << 40; // keeps len and span far from overflow
	int64_t span = 0, len = 0;
	bool has_len = false;
	for (char ch : cigar) {
		if (ch >= '0' && ch <= '9') {
			len = len * 10 + (ch - '0');
			if (len >= MAX_LEN)
				return -1;
			has_len = true;
			continue;
		}
		if (!has_len)
			return -1;
		if (ch == 'M' || ch == 'D' || ch == 'N' || ch == '=' || ch == 'X') {
			span += len;
			if (span >= MAX_LEN)
				return -1;
		}
		else if (ch != 'I' && ch != 'S' && ch != 'H' && ch != 'P')
			return -1;
		len = 0;
		has_len = false;
	}
	return has_len ? -1 : span;
}

}  // namespace

extern "C" {

SEXP gtrackimport_mappedseq(SEXP _track, SEXP _infile, SEXP _pileup, SEXP _binsize, SEXP _cols_order, SEXP _remove_dups,
                            SEXP _paired, SEXP _min_mapq, SEXP _max_fraglen, SEXP _one_based, SEXP _envir)
{
	// The first NUM_USER_COLS columns are the ones cols.order positions; the rest are read only from SAM
	// (MAPQ, CIGAR, PNEXT, TLEN) or from fragment files (END). A column with order 0 is not read.
	enum { SEQ_COL, CHROM_COL, COORD_COL, STRAND_COL, NUM_USER_COLS, MAPQ_COL = NUM_USER_COLS, CIGAR_COL, PNEXT_COL, TLEN_COL, END_COL, NUM_COLS };
	const char *COL_NAMES[NUM_COLS] = { "sequence", "chromosome", "coordinate", "strand", "mapq", "cigar", "pnext", "tlen", "end" };
	// SAM flags: unmapped (0x4), secondary (0x100), QC fail (0x200) and supplementary (0x800) records are never imported
	const uint64_t SAM_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x800;

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

		if (!Rf_isNull(_cols_order) && (Rf_length(_cols_order) != NUM_USER_COLS || (!Rf_isReal(_cols_order) && !Rf_isInteger(_cols_order))))
			verror("cols.order argument must be a vector with %d numeric values", NUM_USER_COLS);

		if (!Rf_isNull(_cols_order) && Rf_isReal(_cols_order)) {
			for (int i = 0; i < NUM_USER_COLS; i++) {
				if (REAL(_cols_order)[i] != (int)REAL(_cols_order)[i])
					verror("cols.order is not an integer");
			}
		}

		if (Rf_length(_remove_dups) > 1 || !Rf_isLogical(_remove_dups))
			verror("remove.dups argument must be a logical value");

		if (Rf_length(_paired) != 1 || !Rf_isLogical(_paired) || LOGICAL(_paired)[0] == NA_LOGICAL)
			verror("paired argument must be TRUE or FALSE");

		if (Rf_length(_min_mapq) != 1 || ((!Rf_isReal(_min_mapq) || REAL(_min_mapq)[0] != (int)REAL(_min_mapq)[0]) && !Rf_isInteger(_min_mapq)))
			verror("min.mapq argument is not an integer");

		if (Rf_length(_max_fraglen) != 1 || ((!Rf_isReal(_max_fraglen) || REAL(_max_fraglen)[0] != (int64_t)REAL(_max_fraglen)[0]) && !Rf_isInteger(_max_fraglen)))
			verror("max.fraglen argument is not an integer");

		if (Rf_length(_one_based) != 1 || !Rf_isLogical(_one_based) || LOGICAL(_one_based)[0] == NA_LOGICAL)
			verror("one.based argument must be TRUE or FALSE");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *infilename = CHAR(STRING_ELT(_infile, 0));
		int pileup = Rf_isReal(_pileup) ? (int)REAL(_pileup)[0] : INTEGER(_pileup)[0];
		double binsize = Rf_isReal(_binsize) ? (int)REAL(_binsize)[0] : INTEGER(_binsize)[0];
		int cols_order[NUM_COLS] = { 0 };
		bool remove_dups = LOGICAL(_remove_dups)[0];
		bool paired = LOGICAL(_paired)[0];
		int min_mapq = Rf_isReal(_min_mapq) ? (int)REAL(_min_mapq)[0] : INTEGER(_min_mapq)[0];
		int64_t max_fraglen = Rf_isReal(_max_fraglen) ? (int64_t)REAL(_max_fraglen)[0] : INTEGER(_max_fraglen)[0];
		bool is_sam_format = Rf_isNull(_cols_order);
		// paired = TRUE on a non-SAM file: tab-delimited fragments, 0-based half-open chrom/start/end in columns 1-3
		// (BED, 10x fragments.tsv.gz). Fragment files are taken as already deduplicated.
		bool is_frag_format = paired && !is_sam_format;
		// SAM POS is 1-based by the spec; a tab-delimited coordinate is 1-based only when one.based is set
		bool one_based = is_sam_format || LOGICAL(_one_based)[0];
		// Tab-delimited files without one.based keep the historical placement: the coordinate is taken as is,
		// and a reverse read is recorded one base past its 5' end (its end-exclusive coordinate).
		bool legacy = !is_sam_format && !is_frag_format && !one_based;
		int64_t rev_end_off = legacy ? 0 : 1; // recorded reverse-read coordinate + rev_end_off = end-exclusive coordinate

		if (is_frag_format && LOGICAL(_one_based)[0])
			verror("one.based is not used for fragment files: their coordinates are 0-based, half-open");

		if (is_sam_format) { // SAM format
			cols_order[SEQ_COL] = 10;
			cols_order[CHROM_COL] = 3;
			cols_order[COORD_COL] = 4;
			cols_order[STRAND_COL] = 2;
			cols_order[MAPQ_COL] = 5;
			cols_order[CIGAR_COL] = 6;
			cols_order[PNEXT_COL] = 8;
			cols_order[TLEN_COL] = 9;
		} else if (is_frag_format) {
			cols_order[CHROM_COL] = 1;
			cols_order[COORD_COL] = 2;
			cols_order[END_COL] = 3;
		} else {
			for (int i = 0; i < NUM_USER_COLS; i++) {
				cols_order[i] = Rf_isReal(_cols_order) ? (int)REAL(_cols_order)[i] : INTEGER(_cols_order)[i];
				if (cols_order[i] <= 0)
					verror("Invalid columns order: %s column's order is %d", COL_NAMES[i], cols_order[i]);
			}
		}

		if (pileup < 0)
			verror("Pileup cannot be negative");

		if (min_mapq < 0)
			verror("min.mapq cannot be negative");

		if (min_mapq > 0 && !is_sam_format)
			verror("min.mapq requires SAM or BAM input (a SAM file needs cols.order = NULL)");

		if (paired) {
			if (pileup)
				verror("pileup is not used with paired = TRUE: each fragment covers its own span. Set pileup to 0.");

			if (binsize <= 0)
				verror("Invalid binsize.\nA dense track is created when paired is TRUE. Binsize must be a positive integer then.");

			if (max_fraglen <= 0)
				verror("max.fraglen must be positive");
		} else {
			if (pileup == 0 && binsize >= 0)
				verror("Invalid binsize.\nSparse track is created when pileup is zero. Binsize must be set to -1 then.");

			if (pileup > 0 && binsize <= 0)
				verror("Invalid binsize.\nDense track is created when pileup is greater than zero. Binsize must be a positive integer then.");
		}

		int num_used_cols = 0;
		for (int i = 0; i < NUM_COLS; i++) {
			if (!cols_order[i])
				continue;
			++num_used_cols;

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
		vector< vector< pair<int64_t, int64_t> > > frags(paired ? num_chroms : 0);
		// lines starting with this char are headers / comments
		int comment_char = is_sam_format ? '@' : is_frag_format ? '#' : 0;

		string dirname = create_track_dir(_envir, track);
		int64_t total_unmapped = 0;
		int64_t total_filtered = 0;
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

			// CRLF line endings: drop the CR
			if (c == '\r')
				continue;

			// skip SAM headers (@) and fragment file comments (#)
			if (!pos && comment_char && c == comment_char) {
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
					bool filtered = false; // a valid record left out by a filter (flags, MAPQ, fragment length)
					bool second_mate = false; // a paired SAM import counts each pair once, by its first mate

					pos = 0;
					for (int i = 0; i < NUM_COLS; i++) {
						if (!str[i].empty())
							num_nonempty_strs++;
					}

					while (num_nonempty_strs == num_used_cols) {
						unordered_map<string, int>::iterator istr2chrom;
						int chrom_idx;
						int64_t coord;
						char *endptr;
						uint64_t flag = 0;
						int64_t span = -1; // reference span of a SAM read, from its CIGAR

						if (is_sam_format) {
							flag = strtoull(str[STRAND_COL].c_str(), &endptr, 0);
							if (*endptr)
								break;
							if (paired && (flag & 0x1) && (flag & 0x80)) {
								second_mate = true;
								break;
							}
						}

						if ((istr2chrom = str2chrom.find(str[CHROM_COL])) == str2chrom.end())
							break;
						chrom_idx = istr2chrom->second;
						int64_t chrom_end = all_genome_intervs[chrom_idx].end;

						coord = strtoll(str[COORD_COL].c_str(), &endptr, 10);
						if (*endptr || coord < (one_based ? 1 : 0)) // SAM POS 0 means no position
							break;
						if (one_based)
							--coord; // 0-based leftmost base
						if (coord >= chrom_end)
							break;

						if (is_frag_format) {
							int64_t end = strtoll(str[END_COL].c_str(), &endptr, 10);
							if (*endptr || end <= coord)
								break;
							if (end - coord > max_fraglen) {
								filtered = true;
								break;
							}
							frags[chrom_idx].emplace_back(coord, min(end, chrom_end));
							mapped = true;
							++num_mapped[chrom_idx];
							break;
						}

						if (is_sam_format) {
							if (flag & 0x4)
								break;

							if (flag & SAM_SKIP_FLAGS) {
								filtered = true;
								break;
							}

							if (min_mapq) {
								long mapq = strtol(str[MAPQ_COL].c_str(), &endptr, 10);
								if (*endptr)
									break;
								if (mapq < min_mapq) {
									filtered = true;
									break;
								}
							}

							// htslib (and so a BAM import) takes a mapped record without a CIGAR ('*') as unmapped;
							// a malformed CIGAR, or one with no aligned reference base, is unusable too
							if ((span = cigar_ref_span(str[CIGAR_COL])) <= 0)
								break;

							if (paired) {
								// One fragment per proper pair, taken from the first mate (as MACS3 BAMPE does):
								// [min(POS, PNEXT), + |TLEN|) in 1-based POS terms. Unpaired reads are filtered.
								if ((flag & 0x3) != 0x3 || (flag & 0x8)) {
									filtered = true;
									break;
								}
								int64_t pnext = strtoll(str[PNEXT_COL].c_str(), &endptr, 10);
								if (*endptr)
									break;
								int64_t tlen = strtoll(str[TLEN_COL].c_str(), &endptr, 10);
								if (*endptr)
									break;
								tlen = llabs(tlen);
								int64_t start = min(coord, pnext - 1);
								if (!tlen || tlen > max_fraglen || start < 0) {
									filtered = true;
									break;
								}
								frags[chrom_idx].emplace_back(start, min(start + tlen, chrom_end));
								mapped = true;
								++num_mapped[chrom_idx];
								break;
							}

							str[STRAND_COL] = flag & 0x10 ? "-" : "+";
						}

						if (str[STRAND_COL] == "+" || str[STRAND_COL] == "F")
							coords[chrom_idx].push_back(coord);
						else if (str[STRAND_COL] == "-" || str[STRAND_COL] == "R") {
							// the 5' end of a reverse read is its rightmost aligned base (one past it in legacy mode)
							if (span <= 0) // a tab-delimited file: no CIGAR
								span = str[SEQ_COL].size();
							// the recorded point must lie on the chromosome; a legacy dense track clips the read instead
							int64_t room = chrom_end - coord;
							if (legacy ? !pileup && span >= room : span > room)
								break;
							coords[num_chroms + chrom_idx].push_back(coord + span - rev_end_off);
						} else
							break;

						mapped = true;
						++num_mapped[chrom_idx];
						break;
					}

					if (filtered)
						total_filtered++;
					else if (!mapped && !second_mate && num_nonempty_strs)
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
			verror("Error while reading file %s", infilename);

		// For PipeSource this drains the child and surfaces non-zero exit
		// (including 127 = "samtools not on PATH"); no-op for other sources.
		src->check_finalize(infilename);

		for (unsigned ichrom = 0; ichrom < num_chroms; ichrom++) {
			char filename[FILENAME_MAX];
			snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), iu.id2chrom(all_genome_intervs[ichrom].chromid).c_str());

			// dense track
			if (pileup || paired) {
				GenomeTrackFixedBin gtrack;
				gtrack.init_write(filename, (unsigned)binsize, all_genome_intervs[ichrom].chromid);

				vector<float> trackvals((uint64_t)ceil(all_genome_intervs[ichrom].end / binsize), 0);

				// adds the fraction of each bin covered by [from_coord, to_coord)
				auto add_coverage = [&](int64_t from_coord, int64_t to_coord) {
					// a reverse read longer than pileup can end past the chromosome: nothing left to add
					if (from_coord >= to_coord)
						return;
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
				};

				if (paired) {
					vector< pair<int64_t, int64_t> > &cur_frags = frags[ichrom];
					sort(cur_frags.begin(), cur_frags.end());

					for (auto ifrag = cur_frags.begin(); ifrag != cur_frags.end(); ++ifrag) {
						if (remove_dups && !is_frag_format && ifrag != cur_frags.begin() && *ifrag == *(ifrag - 1)) {
							++num_dups[ichrom];
							continue;
						}
						add_coverage(ifrag->first, ifrag->second);
					}
					vector< pair<int64_t, int64_t> >().swap(cur_frags);
				}

				for (int strand = 0; strand < 2; strand++) {
					vector<int64_t> &cur_coords = coords[strand * num_chroms + ichrom];
					sort(cur_coords.begin(), cur_coords.end());

					for (vector<int64_t>::const_iterator icoord = cur_coords.begin(); icoord != cur_coords.end(); ++icoord) {
						if (remove_dups && icoord != cur_coords.begin() && *icoord == *(icoord - 1)) {
							++num_dups[ichrom];
							continue;
						}

						int64_t end_coord = strand ? *icoord + rev_end_off : *icoord + pileup;
						add_coverage(max(end_coord - pileup, (int64_t)0), min(end_coord, all_genome_intervs[ichrom].end));
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

		rprotect(total_stat = RSaneAllocVector(REALSXP, 5));
		REAL(total_stat)[0] = total_mapped + total_unmapped + total_filtered; // total_mapped includes the duplicates
		REAL(total_stat)[1] = total_mapped;
		REAL(total_stat)[2] = total_unmapped;
		REAL(total_stat)[3] = total_dups;
		REAL(total_stat)[4] = total_filtered;
		Rf_setAttrib(total_stat, R_NamesSymbol, RSaneAllocVector(STRSXP, 5));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 0, Rf_mkChar("total"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 1, Rf_mkChar("total.mapped"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 2, Rf_mkChar("total.unmapped"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 3, Rf_mkChar("total.dups"));
		SET_STRING_ELT(Rf_getAttrib(total_stat, R_NamesSymbol), 4, Rf_mkChar("total.filtered"));

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
