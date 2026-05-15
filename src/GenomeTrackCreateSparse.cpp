#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <atomic>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"

#include "GenomeTrackSparse.h"

using namespace std;
using namespace rdb;

namespace {

// Write a 4-byte sparse-format signature to `path` via raw POSIX syscalls.
// Used by the parallel empty-chrom path in gtrack_create_sparse; avoids the
// libc stdio FILE* wrapping and per-call umask twiddling that the serial
// GenomeTrackSparse::init_write path goes through. Best-effort: a non-zero
// errno is reported back via `first_errno` so the caller can verror after
// all workers finish.
void write_empty_sparse_signature(const std::string &path, int32_t sig,
                                  std::atomic<int> &first_errno)
{
    int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) {
        int e = errno ? errno : EIO;
        int expected = 0;
        first_errno.compare_exchange_strong(expected, e);
        return;
    }
    ssize_t n = ::write(fd, &sig, sizeof(sig));
    int werr = (n != (ssize_t)sizeof(sig)) ? (errno ? errno : EIO) : 0;
    ::close(fd);
    if (werr) {
        int expected = 0;
        first_errno.compare_exchange_strong(expected, werr);
    }
}

// Read the worker count from R option `gtrack.create.threads`. Default 4
// (empirically saturates NFSv3 CREATE concurrency); cap at 32 to keep
// worst-case thread count sane on misconfigured systems. Must be called on
// the main thread; R API is not thread-safe.
unsigned get_empty_chrom_workers()
{
    SEXP opt = Rf_GetOption1(Rf_install("gtrack.create.threads"));
    int v = -1;
    if (Rf_isInteger(opt) && Rf_length(opt) >= 1)
        v = INTEGER(opt)[0];
    else if (Rf_isReal(opt) && Rf_length(opt) >= 1)
        v = (int)REAL(opt)[0];
    if (v < 0)
        v = 4;
    if (v > 32)
        v = 32;
    return (unsigned)v;
}

} // namespace

extern "C" {

SEXP gtrack_create_sparse(SEXP _track, SEXP _intervs, SEXP _values, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		IntervUtils iu(_envir);
		GIntervals intervs;
		iu.convert_rintervs(_intervs, &intervs, NULL);
		intervs.sort();
		intervs.verify_no_overlaps(iu.get_chromkey());

		if (!Rf_isReal(_values) && !Rf_isInteger(_values))
			verror("Values argument is not numeric");

		if (Rf_length(_values) != (int)intervs.size())
			verror("Number of intervals (%ld) does not match the number of values (%d)", intervs.size(), (int)Rf_length(_values));

		const char *track = CHAR(STRING_ELT(_track, 0));

		string dirname = create_track_dir(_envir, track);
		int cur_chromid = -1;
		set<int> created_chromids;
		GenomeTrackSparse gtrack;
		char filename[FILENAME_MAX];
		GIntervals all_genome_intervs;
		iu.get_all_genome_intervs(all_genome_intervs);

		Progress_reporter progress;
		progress.init(intervs.size(), 100000);

		for (GIntervals::const_iterator iinterv = intervs.begin(); iinterv != intervs.end(); ++iinterv) {
			if (cur_chromid != iinterv->chromid) {
				cur_chromid = iinterv->chromid;
				snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), cur_chromid).c_str());
				gtrack.init_write(filename, cur_chromid);
				created_chromids.insert(cur_chromid);
			}

			if (Rf_isReal(_values))
				gtrack.write_next_interval(*iinterv, REAL(_values)[iu.get_orig_interv_idx(*iinterv)]);
			else
				gtrack.write_next_interval(*iinterv, INTEGER(_values)[iu.get_orig_interv_idx(*iinterv)]);

			progress.report(1);
			check_interrupt();
		}

		// For indexed databases, empty chromosome files are not needed since the track
		// will be converted to indexed format (track.idx + track.dat) automatically.
		// This avoids creating thousands of empty files on network filesystems.
		// For non-indexed databases, we still create empty chromosome files for backward compatibility.
		if (!is_db_indexed(_envir)) {
			vector<string> paths;
			paths.reserve(all_genome_intervs.size());
			for (GIntervals::const_iterator iinterv = all_genome_intervs.begin(); iinterv != all_genome_intervs.end(); ++iinterv) {
				if (created_chromids.find(iinterv->chromid) == created_chromids.end()) {
					snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), iinterv->chromid).c_str());
					paths.emplace_back(filename);
				}
			}

			if (!paths.empty()) {
				const int32_t sig = GenomeTrack::FORMAT_SIGNATURES[GenomeTrack::SPARSE];
				unsigned nworkers = get_empty_chrom_workers();
				if (nworkers == 0)
					nworkers = 1;
				if (nworkers > (unsigned)paths.size())
					nworkers = (unsigned)paths.size();

				std::atomic<int> first_errno{0};

				if (nworkers <= 1) {
					for (const auto &p : paths)
						write_empty_sparse_signature(p, sig, first_errno);
				} else {
					std::atomic<size_t> next{0};
					auto worker = [&]() {
						for (;;) {
							size_t i = next.fetch_add(1, std::memory_order_relaxed);
							if (i >= paths.size())
								return;
							write_empty_sparse_signature(paths[i], sig, first_errno);
						}
					};
					std::vector<std::thread> threads;
					threads.reserve(nworkers - 1);
					for (unsigned t = 0; t < nworkers - 1; ++t)
						threads.emplace_back(worker);
					worker();
					for (auto &th : threads)
						th.join();
				}

				int e = first_errno.load();
				if (e != 0)
					TGLError("Failed to create empty chromosome track file: %s", strerror(e));
			}
		}

		progress.report_last();
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
