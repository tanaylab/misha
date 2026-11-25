/*
 * GenomeTrackImportFigwig.cpp
 *
 * gtrack.import_figwig: the tracks of a figwig file as misha tracks, one per figwig track. A
 * track is read as maximal runs of bases with the same value. A dense track's bin is the mean of
 * its bases with a value, as gtrackimportwig computes it from the records of a bigWig; a sparse
 * track holds one interval per run of a value other than 0 and NaN. PyMisha's
 * gtrack_import_figwig makes the same runs and the same bins.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <set>

#include "FigWig.h"
#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "UnknownChroms.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// GenomeTrackImportWig.cpp
void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);
void report_skipped_chroms(const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);

namespace {

// Tracks imported at a time: each holds an open track file with a 1 MiB buffer
const size_t MAX_BATCH_TRACKS = 64;

// A value as the float a track file holds; beyond float's range is an infinity, not undefined
inline float to_float(double v)
{
	if (v > numeric_limits<float>::max())
		return numeric_limits<float>::infinity();
	if (v < -numeric_limits<float>::max())
		return -numeric_limits<float>::infinity();
	return (float)v;
}

// One chromosome of one track. Takes the values of bases [0, valid) in order, as runs or base by
// base, and joins equal neighbours into maximal runs: dense tracks compare the values as read
// (double), sparse tracks as the floats they are stored as.
class TrackSink {
public:
	TrackSink(const char *filename, int chromid, uint64_t chromsize, unsigned binsize, double defval) :
		m_chromid(chromid), m_chromsize(chromsize), m_binsize(binsize), m_defval(defval)
	{
		if (binsize) {
			m_dense.reset(new GenomeTrackFixedBin);
			m_dense->init_write(filename, binsize, chromid);
			m_bins.reserve(10000);
		} else {
			m_sparse.reset(new GenomeTrackSparse);
			m_sparse->init_write(filename, chromid);
		}
	}

	void add_run(uint64_t start, uint64_t end, double v)
	{
		if (start == end)
			return;
		if (m_open && (m_binsize ? (v == m_cur) : (to_float(v) == to_float(m_cur)))) {
			m_end = end;
			return;
		}
		flush();
		m_open = true;
		m_start = start;
		m_end = end;
		m_cur = v;
	}

	void add_values(uint64_t start, const double *vals, uint64_t n)
	{
		uint64_t i = 0;
		while (i < n) {
			uint64_t j = i + 1;
			if (m_binsize) {
				while (j < n && vals[j] == vals[i])
					++j;
			} else {
				float f = to_float(vals[i]);
				while (j < n && to_float(vals[j]) == f)
					++j;
			}
			add_run(start + i, start + j, vals[i]);
			i = j;
		}
	}

	void finish()
	{
		flush();
		if (m_dense) {
			while (m_bin_start < m_chromsize)
				next_bin();
			write_bins();
			m_dense->flush_writes();
		}
	}

private:
	int                                 m_chromid;
	uint64_t                            m_chromsize;
	unsigned                            m_binsize;
	double                              m_defval;
	unique_ptr<GenomeTrackFixedBin>     m_dense;
	unique_ptr<GenomeTrackSparse>       m_sparse;

	// the run being built
	bool     m_open{false};
	uint64_t m_start{0};
	uint64_t m_end{0};
	double   m_cur{numeric_limits<double>::quiet_NaN()};

	// dense: the bin being summed
	uint64_t      m_bin_start{0};
	double        m_sum{0};
	uint64_t      m_covered{0};
	vector<float> m_bins;

	void flush()
	{
		if (!m_open)
			return;
		m_open = false;
		if (m_dense)
			bin_run(m_start, m_end, m_cur);
		else {
			float f = to_float(m_cur);
			if (f != 0 && !std::isnan(f))
				m_sparse->write_next_interval(GInterval(m_chromid, m_start, m_end, 0), f);
		}
	}

	// gtrackimportwig's DenseBinner, in double: a bin is the sum of value x overlap over the runs
	// with a value, the rest of the bin counting as defval unless it is NaN, over the bases counted
	void bin_run(uint64_t start, uint64_t end, double v)
	{
		while (m_bin_start + m_binsize <= start)
			next_bin();
		while (1) {
			uint64_t bin_end = min(m_bin_start + m_binsize, m_chromsize);
			uint64_t overlap_len = min(end, bin_end) - max(start, m_bin_start);

			if (!std::isnan(v)) {
				m_sum += v * overlap_len;
				m_covered += overlap_len;
			}
			if (end <= bin_end)
				break;
			next_bin();
		}
	}

	void next_bin()
	{
		uint64_t bin_len = min(m_bin_start + m_binsize, m_chromsize) - m_bin_start;

		if (m_covered < bin_len && !std::isnan(m_defval)) {
			m_sum += m_defval * (bin_len - m_covered);
			m_covered = bin_len;
		}
		m_bins.push_back(m_covered ? to_float(m_sum / m_covered) : numeric_limits<float>::quiet_NaN());
		m_sum = 0;
		m_covered = 0;
		m_bin_start += m_binsize;
		if (m_bins.size() >= 10000)
			write_bins();
	}

	void write_bins()
	{
		if (!m_bins.empty())
			m_dense->write_next_bins(m_bins.data(), m_bins.size());
		m_bins.clear();
	}
};

SEXP mk_utf8(const string &s)
{
	return RSaneMkCharLenCE(s.data(), s.size(), CE_UTF8);
}

} // namespace

extern "C" {

// What a figwig file holds: its tracks (name, dtype, group and attributes: a named character
// vector per track, NA for a null value), chromosomes, genome, order and tile
SEXP figwig_info(SEXP _file, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_file) || Rf_length(_file) != 1)
			verror("File argument is not a string");

		FigWig fw;
		fw.open(CHAR(STRING_ELT(_file, 0)));

		const vector<FigWig::Track> &tracks = fw.tracks();
		const vector<FigWig::Chrom> &chroms = fw.chroms();
		const char *fields[] = { "name", "dtype", "group", "attrs", "chroms", "sizes", "genome", "order", "tile" };
		const int num_fields = sizeof(fields) / sizeof(fields[0]);

		SEXP answer = rprotect_ptr(RSaneAllocVector(VECSXP, num_fields));
		SEXP names = rprotect_ptr(RSaneAllocVector(STRSXP, num_fields));
		SEXP name = rprotect_ptr(RSaneAllocVector(STRSXP, tracks.size()));
		SEXP dtype = rprotect_ptr(RSaneAllocVector(STRSXP, tracks.size()));
		SEXP group = rprotect_ptr(RSaneAllocVector(INTSXP, tracks.size()));
		SEXP attrs = rprotect_ptr(RSaneAllocVector(VECSXP, tracks.size()));

		vector<string> left_out;  // attributes that hold a NUL byte, which an R string cannot

		for (size_t i = 0; i < tracks.size(); ++i) {
			const FigWig::Track &t = tracks[i];
			SET_STRING_ELT(name, i, mk_utf8(t.name));
			SET_STRING_ELT(dtype, i, mk_utf8(t.dtype));
			INTEGER(group)[i] = t.group;

			vector<const FigWig::Attr *> kept;
			for (const FigWig::Attr &attr : t.attrs) {
				if (attr.key.find('\0') != string::npos || attr.text.find('\0') != string::npos)
					left_out.push_back("Attribute \"" + string(attr.key.c_str()) + "\" of figwig track \"" + t.name + "\"");
				else
					kept.push_back(&attr);
			}

			SEXP vals = rprotect_ptr(RSaneAllocVector(STRSXP, kept.size()));
			SEXP keys = rprotect_ptr(RSaneAllocVector(STRSXP, kept.size()));
			for (size_t k = 0; k < kept.size(); ++k) {
				SET_STRING_ELT(keys, k, mk_utf8(kept[k]->key));
				SET_STRING_ELT(vals, k, kept[k]->null ? NA_STRING : mk_utf8(kept[k]->text));
			}
			Rf_setAttrib(vals, R_NamesSymbol, keys);
			SET_VECTOR_ELT(attrs, i, vals);
			runprotect(2);  // vals and keys, now held by attrs
		}

		SEXP chrom_names = rprotect_ptr(RSaneAllocVector(STRSXP, chroms.size()));
		SEXP sizes = rprotect_ptr(RSaneAllocVector(REALSXP, chroms.size()));
		for (size_t c = 0; c < chroms.size(); ++c) {
			SET_STRING_ELT(chrom_names, c, mk_utf8(chroms[c].name));
			REAL(sizes)[c] = chroms[c].size;
		}

		SEXP genome = rprotect_ptr(RSaneAllocVector(STRSXP, 1));
		SET_STRING_ELT(genome, 0, fw.has_genome() ? mk_utf8(fw.genome()) : NA_STRING);
		SEXP order = rprotect_ptr(RSaneAllocVector(STRSXP, 1));
		SET_STRING_ELT(order, 0, RSaneMkChar(fw.track_order() ? "tracks" : "positions"));
		SEXP tile = rprotect_ptr(RSaneAllocVector(INTSXP, 1));
		INTEGER(tile)[0] = fw.tile();

		SEXP items[] = { name, dtype, group, attrs, chrom_names, sizes, genome, order, tile };
		for (int i = 0; i < num_fields; ++i) {
			SET_VECTOR_ELT(answer, i, items[i]);
			SET_STRING_ELT(names, i, RSaneMkChar(fields[i]));
		}
		Rf_setAttrib(answer, R_NamesSymbol, names);
		for (const string &what : left_out) {
			string msg = what + " of " + CHAR(STRING_ELT(_file, 0)) + " was left out: it holds a NUL byte, which an attribute cannot";
			add_pending_diagnostic(_envir, "warning", msg.c_str());
		}
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

// Writes the per-chromosome files of figwig tracks _tracks (0-based indices) of _file, for every
// chromosome of the database, into the existing empty directories _dirs. _srcname is what
// messages call the file.
SEXP gtrackimportfigwig(SEXP _dirs, SEXP _file, SEXP _srcname, SEXP _tracks, SEXP _binsize, SEXP _defval, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_dirs) || !Rf_isInteger(_tracks) || Rf_length(_dirs) != Rf_length(_tracks) || !Rf_length(_tracks))
			verror("dirs and tracks must be vectors of the same, positive, length");
		if (!Rf_isString(_file) || Rf_length(_file) != 1 || !Rf_isString(_srcname) || Rf_length(_srcname) != 1)
			verror("File argument is not a string");
		if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
			verror("Binsize argument is not a number");
		if ((!Rf_isReal(_defval) && !Rf_isInteger(_defval)) || Rf_length(_defval) != 1)
			verror("Defval argument is not a number");

		double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
		double defval = Rf_isReal(_defval) ? REAL(_defval)[0] : INTEGER(_defval)[0];
		if (!(dbinsize >= 0) || dbinsize >= 4294967296.0 || dbinsize != floor(dbinsize))
			verror("Invalid value of binsize argument: %g", dbinsize);
		unsigned binsize = (unsigned)dbinsize;
		const char *srcname = CHAR(STRING_ELT(_srcname, 0));

		FigWig fw;
		fw.open(CHAR(STRING_ELT(_file, 0)));

		vector<string> dirs;
		vector<int> tracks;
		for (int i = 0; i < Rf_length(_tracks); ++i) {
			int t = INTEGER(_tracks)[i];
			if (t < 0 || t >= (int)fw.tracks().size())
				verror("Track index %d is out of range", t);
			tracks.push_back(t);
			dirs.push_back(CHAR(STRING_ELT(_dirs, i)));
		}

		IntervUtils iu(_envir);
		const GenomeChromKey &chromkey = iu.get_chromkey();
		const vector<FigWig::Chrom> &fw_chroms = fw.chroms();

		// the file's chromosome of each database chromosome (-1: none)
		vector<int> source(chromkey.get_num_chroms(), -1);
		UnknownChroms unknown;
		uint64_t num_matched = 0;
		uint64_t num_mismatched = 0;
		string first_mismatch;

		for (size_t c = 0; c < fw_chroms.size(); ++c) {
			int chromid = chromkey.find_chromid(fw_chroms[c].name);
			if (chromid < 0) {
				unknown.record(fw_chroms[c].name);
				continue;
			}
			if (source[chromid] >= 0)
				verror("Chromosomes %s and %s of %s both name chromosome %s of the database", fw_chroms[source[chromid]].name.c_str(),
				       fw_chroms[c].name.c_str(), srcname, chromkey.id2chrom(chromid).c_str());
			source[chromid] = c;
			++num_matched;
			if (fw_chroms[c].size != chromkey.get_chrom_size(chromid) && !num_mismatched++) {
				char buf[1000];
				snprintf(buf, sizeof(buf), "%s, %llu bp in the file, %llu bp in the database", chromkey.id2chrom(chromid).c_str(),
				         (unsigned long long)fw_chroms[c].size, (unsigned long long)chromkey.get_chrom_size(chromid));
				first_mismatch = buf;
			}
		}
		verify_chroms_matched(chromkey, srcname, unknown, num_matched, _envir);

		uint64_t tile = fw.tile();
		uint64_t num_batches = (tracks.size() + MAX_BATCH_TRACKS - 1) / MAX_BATCH_TRACKS;
		Progress_reporter progress;
		progress.init(num_batches * chromkey.get_num_chroms(), 1);

		vector<double> buf;
		vector<uint32_t> positions;
		vector<double> values;
		char filename[FILENAME_MAX];

		for (size_t b0 = 0; b0 < tracks.size(); b0 += MAX_BATCH_TRACKS) {
			size_t b1 = min(tracks.size(), b0 + MAX_BATCH_TRACKS);

			// the tracks of the batch by group, for position order
			map<int, vector<size_t>> by_group;
			for (size_t i = b0; i < b1; ++i)
				by_group[fw.tracks()[tracks[i]].group].push_back(i);

			for (int chromid = 0; chromid < (int)chromkey.get_num_chroms(); ++chromid) {
				uint64_t chromsize = chromkey.get_chrom_size(chromid);
				int c = source[chromid];
				uint64_t valid = c < 0 ? 0 : min(chromsize, fw_chroms[c].size);
				vector<unique_ptr<TrackSink>> sinks(b1 - b0);

				for (size_t i = b0; i < b1; ++i) {
					snprintf(filename, sizeof(filename), "%s/%s", dirs[i].c_str(), GenomeTrack::get_1d_filename(chromkey, chromid).c_str());
					sinks[i - b0].reset(new TrackSink(filename, chromid, chromsize, binsize, defval));
				}

				if (valid && !fw.track_order()) {
					for (const auto &grp : by_group) {
						vector<int> slots;
						for (size_t i : grp.second)
							slots.push_back(fw.tracks()[tracks[i]].slot);
						buf.resize(slots.size() * tile);
						for (uint64_t t = 0; t * tile < valid; ++t) {
							uint64_t start = t * tile;
							uint64_t n = min(tile, valid - start);
							if (fw.read_tile(grp.first, fw_chroms[c].first_tile + t, slots, buf.data())) {
								for (size_t k = 0; k < slots.size(); ++k)
									sinks[grp.second[k] - b0]->add_values(start, buf.data() + k * tile, n);
							} else {
								for (size_t k = 0; k < slots.size(); ++k)
									sinks[grp.second[k] - b0]->add_run(start, start + n, 0.);
							}
							if (!(t % 64))
								check_interrupt();
						}
					}
				} else if (valid) {
					for (size_t i = b0; i < b1; ++i) {
						TrackSink &sink = *sinks[i - b0];
						uint64_t pos = 0;
						uint64_t num_blocks = fw.num_blocks(tracks[i], c);
						for (uint64_t iblock = 0; iblock < num_blocks && pos < valid; ++iblock) {
							fw.read_block(tracks[i], c, iblock, positions, values);
							for (size_t k = 0; k < positions.size() && positions[k] < valid; ++k) {
								sink.add_run(pos, positions[k], 0.);
								sink.add_run(positions[k], positions[k] + 1, values[k]);
								pos = positions[k] + 1;
							}
							check_interrupt();
						}
						sink.add_run(pos, valid, 0.);
					}
				}

				for (auto &sink : sinks)
					sink->finish();
				progress.report(1);
				check_interrupt();
			}
		}
		progress.report_last();

		// Only now, when nothing can fail, is it safe to leave pending messages for .gcall()
		if (num_mismatched) {
			char msg[2000];
			snprintf(msg, sizeof(msg),
			         "%llu chromosome(s) of %s differ in length from the database (the first: %s): bases past the database's end "
			         "were dropped, and bases past the file's end have no value",
			         (unsigned long long)num_mismatched, srcname, first_mismatch.c_str());
			add_pending_diagnostic(_envir, "warning", msg);
		}
		report_skipped_chroms(srcname, unknown, num_matched, _envir);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const exception &e) {
		rerror("%s", e.what());
	}
	return R_NilValue;
}

}
