#include <cstdint>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <set>
#include <string>
#include <vector>

#include "BinFinder.h"
#include "GenomeTrackFixedBin.h"
#include "TrackIndex.h"

#include "rdbinterval.h"
#include "rdbutils.h"
#include "TrackExpressionFixedBinIterator.h"
#include "TrackExpressionScanner.h"

using namespace std;
using namespace rdb;

namespace {

// gtrack.modify is the only writer in the package that used to edit a track in
// place ("rb+" straight into the live file). An interrupt, an expression error
// or a full disk left the track durably half-old/half-new under its real name,
// structurally valid, with nothing to mark it as damaged.
//
// It is now staged like every other writer: the data files the modification
// touches are copied into a hidden staging directory, the new values are
// written there, and the result is committed by renaming the staged files back
// over the originals. Until the commit the live track is untouched, so any
// failure (including SIGKILL) leaves the previous track intact; the R wrapper
// removes the staging directory.
//
// Whole-track atomicity is only as good as the number of files: on an indexed
// database there is exactly one data file (track.dat), so the commit is a
// single rename() and is atomic. On a per-chromosome track the commit is one
// rename per touched chromosome, performed back to back with no interruptible
// work in between; each file is still individually consistent.
class Stager {
public:
	Stager(const string &track_dir, const string &stage_dir) :
		m_track_dir(track_dir), m_stage_dir(stage_dir), m_indexed(GenomeTrack::get_track_index(track_dir) != nullptr) {}

	bool enabled() const { return !m_stage_dir.empty(); }

	// Returns the path that should be opened for update for this chromosome.
	// Copies whatever backing file the chromosome lives in on first use.
	string stage(const string &chrom_filename)
	{
		if (m_indexed) {
			// All chromosomes share track.dat; the index is needed alongside it
			// so that the staged copy resolves the same offsets.
			stage_file("track.idx", false);
			stage_file("track.dat", true);
		} else
			stage_file(chrom_filename, true);

		return m_stage_dir + "/" + chrom_filename;
	}

	// Renames the staged data files over the originals. Nothing else may be
	// interposed here: this is the commit.
	void commit()
	{
		for (vector<string>::const_iterator i = m_committable.begin(); i != m_committable.end(); ++i) {
			string src = m_stage_dir + "/" + *i;
			string dst = m_track_dir + "/" + *i;
			if (rename(src.c_str(), dst.c_str()))
				verror("Failed to commit %s to %s: %s", src.c_str(), dst.c_str(), strerror(errno));
		}
		m_committable.clear();
	}

private:
	string         m_track_dir;
	string         m_stage_dir;
	bool           m_indexed;
	set<string>    m_staged;
	vector<string> m_committable;

	void stage_file(const string &name, bool committable)
	{
		if (!m_staged.insert(name).second)
			return;
		copy_file(m_track_dir + "/" + name, m_stage_dir + "/" + name);
		if (committable)
			m_committable.push_back(name);
	}

	static void copy_file(const string &src, const string &dst)
	{
		int in = ::open(src.c_str(), O_RDONLY);
		if (in < 0)
			verror("Cannot open %s: %s", src.c_str(), strerror(errno));

		int out = ::open(dst.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
		if (out < 0) {
			int err = errno;
			::close(in);
			verror("Cannot create %s: %s", dst.c_str(), strerror(err));
		}

		// The staged file is renamed over the original, so it must carry the
		// original's permissions rather than whatever the umask happens to be.
		struct stat st;
		if (!::fstat(in, &st))
			(void)::fchmod(out, st.st_mode & 07777);

		vector<char> buf(1 << 20);
		while (1) {
			ssize_t nread = ::read(in, &buf[0], buf.size());
			if (nread < 0) {
				if (errno == EINTR)
					continue;
				int err = errno;
				::close(in);
				::close(out);
				verror("Reading %s: %s", src.c_str(), strerror(err));
			}
			if (!nread)
				break;

			ssize_t written = 0;
			while (written < nread) {
				ssize_t n = ::write(out, &buf[written], nread - written);
				if (n < 0) {
					if (errno == EINTR)
						continue;
					int err = errno;
					::close(in);
					::close(out);
					verror("Writing %s: %s", dst.c_str(), strerror(err));
				}
				written += n;
			}

			// A multi-GB track.dat must stay interruptible. check_interrupt()
			// throws TGLException, which unwinds through the guards below and
			// leaves the staging directory for R to remove.
			try {
				check_interrupt();
			} catch (...) {
				::close(in);
				::close(out);
				throw;
			}
		}

		::close(in);
		// close() is where a deferred ENOSPC/EDQUOT surfaces on some filesystems.
		if (::close(out))
			verror("Writing %s: %s", dst.c_str(), strerror(errno));
	}
};

}

extern "C" {

SEXP gtrack_modify(SEXP _track, SEXP _track_expr, SEXP _intervals, SEXP _iterator_policy, SEXP _workdir, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		// check the arguments
		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(_track_expr) || Rf_length(_track_expr) != 1)
			verror("Track expression argument is not a string vector");

		string stage_dir;
		if (!Rf_isNull(_workdir)) {
			if (!Rf_isString(_workdir) || Rf_length(_workdir) != 1)
				verror("Staging directory argument is not a string");
			stage_dir = CHAR(STRING_ELT(_workdir, 0));
		}

		IntervUtils iu(_envir);
		GIntervalsFetcher1D *intervals = NULL;
		iu.convert_rintervs(_intervals, &intervals, NULL);
		unique_ptr<GIntervalsFetcher1D> intervals1d_guard(intervals);
		intervals->sort();
		intervals->unify_overlaps(false);
		const char *track = CHAR(STRING_ELT(_track, 0));
		string trackpath = track2path(_envir, track);
		char filename[FILENAME_MAX];
		TrackExprScanner scanner(iu);
		GenomeTrackFixedBin gtrack;
		GInterval last_interval(-1, -1, -1, -1);

		if (GenomeTrack::get_type(trackpath.c_str(), iu.get_chromkey()) != GenomeTrack::FIXED_BIN)
			verror("Cannot modify track %s: modification is supported only for dense tracks", track);

		Stager stager(trackpath, stage_dir);

		scanner.begin(_track_expr, intervals, NULL, _iterator_policy);

		if (scanner.get_iterator()->get_type() != TrackExpressionIteratorBase::FIXED_BIN)
			verror("gtrack.modify() requires the iterator policy to be fixed bin.\n");

		for (; !scanner.isend(); scanner.next()) {
			const GInterval &cur_interval = scanner.last_interval1d();

			if (last_interval.chromid != cur_interval.chromid) {
				string resolved = GenomeTrack::find_existing_1d_filename(iu.get_chromkey(), trackpath, cur_interval.chromid);
				snprintf(filename, sizeof(filename), "%s/%s", trackpath.c_str(), resolved.c_str());
				gtrack.init_read(filename, cur_interval.chromid);
				if (gtrack.get_bin_size() != ((TrackExpressionFixedBinIterator *)scanner.get_iterator())->get_bin_size())
					verror("Cannot modify track %s: iterator policy must be set to the bin size of the track (%d).\n", track, gtrack.get_bin_size());

				// The values go to the staging copy; the live track keeps the
				// old ones until commit() below.
				if (stager.enabled()) {
					string staged = stager.stage(resolved);
					gtrack.init_update(staged.c_str(), cur_interval.chromid);
				} else
					gtrack.init_update(filename, cur_interval.chromid);
			}

			if (last_interval.chromid != cur_interval.chromid || last_interval.end != cur_interval.start)
				gtrack.goto_bin((uint64_t)(cur_interval.start / gtrack.get_bin_size()));

			double val = scanner.last_real(0);
			gtrack.write_next_bin(val);
			last_interval = cur_interval;
		}

		// Everything is written and reported before anything is published: a
		// short write must fail here, not silently at fclose() after commit.
		gtrack.flush_writes();
		stager.commit();
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}

}
