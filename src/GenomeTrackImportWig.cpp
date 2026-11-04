#include <cstdint>
#include <cmath>
#include <cstring>
#include <set>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "BigWig.h"
#include "GenomeArraysCsv.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackSparse.h"
#include "Wig.h"

using namespace std;
using namespace rdb;

// Formats names as "a, b, c (and N more)". "N+" when the count itself was capped.
static string format_chrom_names(const vector<string> &names, uint64_t total, bool truncated)
{
	string res;

	for (size_t i = 0; i < names.size(); ++i) {
		if (i)
			res += ", ";
		res += names[i];
	}

	if (total > names.size()) {
		char buf[100];
		snprintf(buf, sizeof(buf), "%s(and %llu%s more)", names.empty() ? "" : " ",
		         (unsigned long long)(total - names.size()), truncated ? "+" : "");
		res += buf;
	}
	return res;
}

static string format_db_chroms(const GenomeChromKey &chromkey)
{
	vector<string> names;
	size_t shown = min((size_t)chromkey.get_num_chroms(), (size_t)UnknownChroms::MAX_REPORTED);

	for (size_t i = 0; i < shown; ++i)
		names.push_back(chromkey.id2chrom((int)i));
	return format_chrom_names(names, chromkey.get_num_chroms(), false);
}

static string format_alias_hint(SEXP envir)
{
	string hint = "Chromosome aliases are resolved automatically (chr1 <-> 1, M <-> MT, plus ";

	hint += get_groot(envir);
	hint += "/chrom_aliases.tsv if present). If the file's names should map to database chromosomes, "
	        "add them to that file (tab separated: canonical<TAB>alias<TAB>source) and reload the database "
	        "with gdb.init(rescan = TRUE), or rename the chromosomes in the file.";
	return hint;
}

// A file that named chromosomes and matched none of them contributes nothing to the track:
// whatever the track would end up holding (all NaN, or defval everywhere) comes from nowhere
// in the file. Refuse before writing anything.
//
// A file that matched some of them is a different matter, and is left to
// report_skipped_chroms below. The distinction the importer can actually make here is
// between names in the file that resolved and names that did not - not between database
// chromosomes the file covered and ones it did not. A chr1-only bedGraph dropped nothing,
// so it is silent; it is unmatched names that prove data in the file never reached the track.
void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname,
                                  const UnknownChroms &unknown, uint64_t num_matched, SEXP envir)
{
	if (unknown.empty() || num_matched)
		return;

	verror("None of the chromosome names in %s exist in the genome database: no data from the file would reach the track.\n"
	       "File chromosomes: %s.\n"
	       "Database chromosomes: %s.\n"
	       "%s",
	       fname, format_chrom_names(unknown.names(), unknown.num(), unknown.truncated()).c_str(),
	       format_db_chroms(chromkey).c_str(), format_alias_hint(envir).c_str());
}

// Reports the names dropped by an import that did produce data. Which channel it goes to
// depends on what was dropped, because the two cases are not equally alarming:
//
//   - a scaffold / patch / unplaced contig (chrUn_*, *_random, GL000220.1): a database built
//     from the main chromosomes deliberately does not have it, no alias can fix it, and this
//     is what every whole-genome bigWig looks like. A message: visible and in order, and it
//     cannot be aggregated away into "there were 50 or more warnings".
//   - a primary chromosome (chr7, X): the database was supposed to have that one, so the
//     naming is probably wrong. A warning, with the alias hint. Mitochondria are the first
//     kind, not this one - see UnknownChroms.
void report_skipped_chroms(const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir)
{
	if (unknown.empty())
		return;

	char msg[10000];
	const char *severity;

	if (unknown.primary_names().empty()) {
		snprintf(msg, sizeof(msg),
		         "%llu%s chromosome name(s) in %s do not exist in the genome database and were skipped: %s. "
		         "Data for the remaining %llu chromosome(s) was imported.",
		         (unsigned long long)unknown.num(), unknown.truncated() ? "+" : "", fname,
		         format_chrom_names(unknown.names(), unknown.num(), unknown.truncated()).c_str(),
		         (unsigned long long)num_matched);
		severity = "message";
	} else {
		snprintf(msg, sizeof(msg),
		         "%llu%s chromosome name(s) in %s do not exist in the genome database and were skipped, among them "
		         "primary chromosome(s): %s. Data for the remaining %llu chromosome(s) was imported. %s",
		         (unsigned long long)unknown.num(), unknown.truncated() ? "+" : "", fname,
		         format_chrom_names(unknown.primary_names(), unknown.num_primary(), false).c_str(),
		         (unsigned long long)num_matched, format_alias_hint(envir).c_str());
		severity = "warning";
	}

	// Queued for .gcall() to raise once the call has returned (see add_pending_diagnostic
	// for why it cannot be raised here), and only once the import can no longer fail.
	add_pending_diagnostic(envir, severity, msg);
}

// Turns records into the bins of a dense track as they arrive. Records must be sorted and must
// not overlap. A bin is the mean over its bases: a base under a record takes the record's value,
// any other base (or one under a NaN record) takes defval, and a bin left with no value is NaN.
class DenseBinner {
public:
	DenseBinner(GenomeTrackFixedBin &track, uint64_t binsize, uint64_t chromsize, double defval, Progress_reporter &progress) :
		m_track(track), m_binsize(binsize), m_chromsize(chromsize), m_defval(defval), m_progress(progress) { m_bins.reserve(10000); }

	void add(uint64_t start, uint64_t end, float v)
	{
		// The parsers guarantee this; a record that broke it would loop past the chromosome end.
		if (start < m_bin_start || start > end || end > m_chromsize)
			verror("Internal error: import record %llu-%llu is out of order or outside the chromosome",
			       (unsigned long long)start, (unsigned long long)end);

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

	void finish()
	{
		while (m_bin_start < m_chromsize)
			next_bin();
		write();
	}

private:
	GenomeTrackFixedBin &m_track;
	uint64_t             m_binsize;
	uint64_t             m_chromsize;
	double               m_defval;
	Progress_reporter   &m_progress;
	uint64_t             m_bin_start{0};
	double               m_sum{0};
	uint64_t             m_covered{0};
	vector<float>        m_bins;

	void next_bin()
	{
		uint64_t bin_len = min(m_bin_start + m_binsize, m_chromsize) - m_bin_start;

		if (m_covered < bin_len && !std::isnan(m_defval)) {
			m_sum += m_defval * (bin_len - m_covered);
			m_covered = bin_len;
		}

		m_bins.push_back(m_covered ? (float)(m_sum / m_covered) : numeric_limits<float>::quiet_NaN());
		m_sum = 0;
		m_covered = 0;
		m_bin_start += m_binsize;

		if (m_bins.size() >= 10000) {
			write();
			m_progress.report(0);
			check_interrupt();
		}
	}

	void write()
	{
		if (!m_bins.empty())
			m_track.write_next_bins(m_bins.data(), m_bins.size());
		m_bins.clear();
	}
};

extern "C" {

// The chromosome check of the importers above, for gtrack.import of a BED file: its records are passed to
// gtrack_create_dense / gtrack_create_sparse, which accept only names the database has. Returns for each of the
// (distinct) names in _chroms whether the database has it, aliases resolved. With _report FALSE it fails as
// verify_chroms_matched does when none of the names exists. With _report TRUE, called once the track is created,
// it queues the report of the names that do not exist (report_skipped_chroms).
SEXP gimport_check_chroms(SEXP _chroms, SEXP _srcname, SEXP _report, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_chroms))
			verror("Chromosomes argument is not a vector of strings");

		if (!Rf_isString(_srcname) || Rf_length(_srcname) != 1)
			verror("Source name argument is not a string");

		if (!Rf_isLogical(_report) || Rf_length(_report) != 1)
			verror("Report argument is not a logical value");

		IntervUtils iu(_envir);
		const GenomeChromKey &chromkey = iu.get_chromkey();
		const char *srcname = CHAR(STRING_ELT(_srcname, 0));
		UnknownChroms unknown;
		set<int> matched;
		SEXP answer = rprotect_ptr(RSaneAllocVector(LGLSXP, Rf_length(_chroms)));

		for (int i = 0; i < Rf_length(_chroms); ++i) {
			const char *chrom = CHAR(STRING_ELT(_chroms, i));
			int chromid = chromkey.find_chromid(chrom);

			LOGICAL(answer)[i] = chromid >= 0;
			if (chromid >= 0)
				matched.insert(chromid);
			else
				unknown.record(chrom);
		}

		if (LOGICAL(_report)[0])
			report_skipped_chroms(srcname, unknown, matched.size(), _envir);
		else
			verify_chroms_matched(chromkey, srcname, unknown, matched.size(), _envir);
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}

	return R_NilValue;
}

SEXP gtrackimportwig(SEXP _track, SEXP _wig, SEXP _srcname, SEXP _binsize, SEXP _defvalue, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_wig) || Rf_length(_wig) != 1)
			verror("Wig argument is not a string");

		if (!Rf_isString(_srcname) || Rf_length(_srcname) != 1)
			verror("Source name argument is not a string");

		if ((!Rf_isReal(_binsize) && !Rf_isInteger(_binsize)) || Rf_length(_binsize) != 1)
			verror("Binsize argument is not a number");

		if ((!Rf_isReal(_defvalue) && !Rf_isInteger(_defvalue)) || Rf_length(_defvalue) != 1)
			verror("Defvalue argument is not a number");

		const char *track = CHAR(STRING_ELT(_track, 0));
		const char *fname = CHAR(STRING_ELT(_wig, 0));
		// What the user asked to import. gtrack.import unzips a .gz / .zip into a
		// temporary file that is deleted before the call returns, so messages that name
		// fname would point at a path that no longer exists (and, in gtrack.import_set,
		// disagree with the file name reported as failed).
		const char *srcname = CHAR(STRING_ELT(_srcname, 0));
		double dbinsize = Rf_isReal(_binsize) ? REAL(_binsize)[0] : INTEGER(_binsize)[0];
		unsigned binsize = (unsigned)dbinsize;
		double defvalue = Rf_isReal(_defvalue) ? REAL(_defvalue)[0] : INTEGER(_defvalue)[0];

		if (dbinsize < 0 || binsize != dbinsize)
			verror("Invalid value of binsize argument: %g\n", dbinsize);

		string dirname = create_track_dir(_envir, track);
		IntervUtils iu(_envir);
		Wig wig;
		BigWig bigwig;
		GenomeArraysCsv csv;
		GIntervals data;
		vector<float> vals;
		char filename[FILENAME_MAX];
		bool is_bigwig = BigWig::is_bigwig(fname);
		bool is_csv = false;
		string wig_err_msg;

		if (is_bigwig)
			bigwig.init(iu.get_chromkey(), fname);
		else {
			try {
				wig.init(iu.get_chromkey(), fname, true);
			} catch (TGLException &e) {
				if (e.type() != typeid(Wig) || e.code() == Wig::FILE_ERROR)
					throw e;
				wig_err_msg = e.msg();
				is_csv = true;
			}
		}

		if (is_csv) {
			try {
				csv.init(fname, iu.get_chromkey());
				if (csv.get_colnames().size() != 1)
					verror("More than one value column appears in file %s", fname);
			} catch (TGLException &e) {
				if (e.type() != typeid(GenomeArraysCsv) || e.code() == GenomeArraysCsv::FILE_ERROR)
					throw e;
				verror("Unrecognized format of file %s.\nWIG parser error: %s\nCSV parser error: %s",
				       fname, wig_err_msg.c_str(), e.msg());
			}

		}

		const UnknownChroms &unknown_chroms = is_bigwig ? bigwig.get_unknown_chroms() :
			is_csv ? csv.get_unknown_chroms() : wig.get_unknown_chroms();
		uint64_t num_matched_chroms = is_bigwig ? bigwig.get_num_matched_chroms() :
			is_csv ? csv.get_num_matched_chroms() : wig.get_num_matched_chroms();

		verify_chroms_matched(iu.get_chromkey(), srcname, unknown_chroms, num_matched_chroms, _envir);

		GIntervals all_genome_intervs;
		iu.get_all_genome_intervs(all_genome_intervs);

		Progress_reporter progress;
		progress.init(iu.get_chromkey().get_num_chroms(), 1);

		for (int chromid = 0; chromid < (int)iu.get_chromkey().get_num_chroms(); ++chromid) {
			uint64_t chromsize = iu.get_chromkey().get_chrom_size(chromid);
			GenomeTrackFixedBin dense;
			GenomeTrackSparse sparse;
			DenseBinner binner(dense, binsize, chromsize, defvalue, progress);
			uint64_t num_recs = 0;

			check_interrupt();

			snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), chromid).c_str());
			if (binsize)
				dense.init_write(filename, binsize, chromid);
			else
				sparse.init_write(filename, chromid);

			auto add = [&](const GIntervals &intervals) {
				for (GIntervals::const_iterator iinterval = intervals.begin(); iinterval != intervals.end(); ++iinterval) {
					float v;

					if (is_csv) {
						csv.get_sliced_vals(iinterval, vals);
						v = vals.front();
					} else
						memcpy(&v, &iinterval->udata, sizeof(float));

					if (binsize)
						binner.add(iinterval->start, iinterval->end, v);
					else {
						sparse.write_next_interval(*iinterval, v);
						if (++num_recs % 10000 == 0) {
							progress.report(0);
							check_interrupt();
						}
					}
				}
			};

			if (is_csv)
				add(csv.get_intervals(chromid));
			else if (is_bigwig) {
				// one data block at a time: a whole chromosome of records can take gigabytes
				for (uint64_t iblock = 0; bigwig.get_data(chromid, iblock, data); ++iblock)
					add(data);
			} else {
				wig.get_data(chromid, data);
				add(data);
			}

			if (binsize)
				binner.finish();
			progress.report(1);
		}
		progress.report_last();

		// Only now, when nothing downstream can throw, is it safe to leave a pending
		// message behind for .gcall() to raise.
		report_skipped_chroms(srcname, unknown_chroms, num_matched_chroms, _envir);
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    } catch (const std::exception &e) {
        rerror("%s", e.what());
    }

	return R_NilValue;
}

}
