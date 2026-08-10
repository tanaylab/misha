#include <cstdint>
#include <cstring>
#include <ctime>
#include "GIntervalsBigSet1D.h"
#include "GenomeSeqFetch.h"
#include "GenomeIndex.h"
#include "rdbinterval.h"
#include "rdbutils.h"

using namespace std;
using namespace rdb;

// Cost of one interval expressed in "range" units so it can be summed with interval lengths:
// a seek plus BufferedFile's 128 KB read-ahead window (see BufferedFile::read). On a cold page
// cache that dominates no matter how few bases were asked for, so balancing kids by bases alone
// would badly skew a many-small-intervals workload.
static const uint64_t SEQ_INTERVAL_OVERHEAD = 131072;

// Whether distributing pays depends almost entirely on whether the pages are already cached: a
// cold read costs milliseconds, a cached one well under a microsecond, and forking a process
// costs ~10 ms. That is a four-orders-of-magnitude spread, so no static interval count can
// decide it -- instead the first few intervals are read serially and timed.
static const uint64_t SEQ_PROBE_INTERVALS = 16;
static const uint64_t SEQ_MIN_INTERVALS4MT = 500;    // below this even a cold read cannot repay a fork
static const double   SEQ_PROBE_MIN_USEC = 100.;     // per interval, above which we are I/O bound
static const double   SEQ_MIN_USEC4PROCESS = 100000.;  // give each kid ~10x the cost of forking it

// One shared-memory record: [uint64 original row index][uint64 sequence length][length bytes]
static const uint64_t SEQ_REC_HEADER = 2 * sizeof(uint64_t);

// Reads are served out of BufferedFile's read-ahead window, so the unit of I/O cost is a window
// miss, not an interval: 5000 tiled 200 bp intervals over 1 Mb touch 8 windows, not 5000. Counting
// the windows a (sorted) interval range spans lets the probe below extrapolate by window instead of
// by interval. Timing per interval and scaling by interval count overestimates a clustered set by
// however many intervals share a window - measured 4.9x SLOWER than sequential on a cold tiled scan
// before this was accounted for.
static uint64_t count_readahead_windows(GIntervals::const_iterator ibegin, GIntervals::const_iterator iend)
{
	uint64_t windows = 0;
	int last_chromid = -1;
	int64_t last_window = -1;

	for (GIntervals::const_iterator iinterv = ibegin; iinterv != iend; ++iinterv) {
		int64_t first = iinterv->start / (int64_t)SEQ_INTERVAL_OVERHEAD;
		int64_t last = max(iinterv->end - 1, iinterv->start) / (int64_t)SEQ_INTERVAL_OVERHEAD;

		if (iinterv->chromid != last_chromid) {
			windows += last - first + 1;
			last_chromid = iinterv->chromid;
			last_window = last;
		} else {
			// intervals are sorted, so only the part past the previous interval is a fresh miss
			int64_t from = max(first, last_window + 1);

			if (last >= from)
				windows += last - from + 1;
			if (last > last_window)
				last_window = last;
		}
	}

	return windows;
}

// gseq.extract.probe.usec overrides SEQ_PROBE_MIN_USEC: 0 always distributes (which is how the
// tests reach the parallel path on a cached example genome), a huge value never does.
static double get_probe_min_usec()
{
	SEXP opt = Rf_GetOption1(Rf_install("gseq.extract.probe.usec"));

	if (Rf_isReal(opt) && Rf_length(opt) == 1 && !ISNAN(REAL(opt)[0]))
		return REAL(opt)[0];
	if (Rf_isInteger(opt) && Rf_length(opt) == 1 && INTEGER(opt)[0] != NA_INTEGER)
		return (double)INTEGER(opt)[0];
	return SEQ_PROBE_MIN_USEC;
}

extern "C" {

SEXP gseqread(SEXP _intervals, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		IntervUtils iu(_envir);
		GIntervalsFetcher1D *intervals = NULL;
		iu.convert_rintervs(_intervals, &intervals, NULL);
		unique_ptr<GIntervalsFetcher1D> intervals_guard(intervals);
		intervals->sort();

		uint64_t num_intervs = intervals->size();
		if (!num_intervs)
			return R_NilValue;

		SEXP answer;
		rprotect(answer = RSaneAllocVector(STRSXP, num_intervs));

		string seqdir = string(rdb::get_groot(_envir)) + "/seq";
		GenomeSeqFetch seqfetch;
		seqfetch.set_seqdir(seqdir);

		vector<char> buf;
		// Reserve based on a heuristic: start with 256KB and grow as needed
		buf.reserve(262144);

		uint64_t seqlen = 0;
		uint64_t num_probed = 0;
		int num_kids = 0;

		// Only in-memory interval sets are distributed; big sets stay on the serial path.
		GIntervals *plain_intervals = dynamic_cast<GIntervals *>(intervals);

		// A kid must not fork again: the nested distribute_task would find the outer job's shared
		// memory already mapped and reuse its record size and offsets.
		if (iu.get_multitasking() && !RdbInitializer::is_kid() &&
			plain_intervals && num_intervs >= SEQ_MIN_INTERVALS4MT) {
			struct timespec probe_start, probe_end;

			// The probed sequences are kept, so a probe that decides against forking costs nothing.
			num_probed = min(num_intervs, SEQ_PROBE_INTERVALS);
			clock_gettime(CLOCK_REALTIME, &probe_start);

			for (uint64_t i = 0; i < num_probed; ++i) {
				const GInterval &interval = (*plain_intervals)[i];

				seqfetch.read_interval(interval, iu.get_chromkey(), buf);
				seqlen += buf.size();
				SET_STRING_ELT(answer, iu.get_orig_interv_idx(interval),
							   Rf_mkCharLenCE(buf.empty() ? "" : &buf[0], (int)buf.size(), CE_NATIVE));
			}

			clock_gettime(CLOCK_REALTIME, &probe_end);
			iu.verify_max_data_size(seqlen, "Result sequence");

			// CLOCK_REALTIME to match the rest of misha (macOS builds shim it); a backwards NTP
			// step would only cost us one missed chance to distribute, so clamp and move on.
			double probe_usec = (probe_end.tv_sec - probe_start.tv_sec) * 1000000. +
								(probe_end.tv_nsec - probe_start.tv_nsec) / 1000.;

			if (probe_usec < 0.)
				probe_usec = 0.;

			// Cost per read-ahead window, extrapolated over the windows the whole set spans, then
			// amortised back over the intervals. A clustered set has far fewer windows than
			// intervals and so lands below the threshold, which is what keeps it sequential.
			uint64_t probe_windows = count_readahead_windows(plain_intervals->begin(),
															 plain_intervals->begin() + num_probed);
			uint64_t total_windows = count_readahead_windows(plain_intervals->begin(),
															 plain_intervals->end());
			double usec_per_interv = probe_usec / (double)max<uint64_t>(1, probe_windows) *
									 (double)total_windows / (double)num_intervs;

			double probe_min_usec = get_probe_min_usec();

			if (usec_per_interv >= probe_min_usec) {
				// A zero threshold means "always distribute"; otherwise scale the number of kids to
				// the measured cost so each one gets meaningfully more work than forking it costs.
				double kids_by_cost = usec_per_interv * num_intervs / SEQ_MIN_USEC4PROCESS;
				uint64_t desired_kids = probe_min_usec > 0 ?
					(uint64_t)min(kids_by_cost, (double)MAX_KIDS) : (uint64_t)MAX_KIDS;

				num_kids = iu.prepare4multitasking_whole_intervals(plain_intervals->begin() + num_probed,
																   plain_intervals->end(),
																   SEQ_INTERVAL_OVERHEAD, desired_kids);
			}
		}

		if (num_kids > 1) {
			uint64_t num_distributed = num_intervs - num_probed;

			// A sequence is never longer than its interval, so this bounds the shared memory.
			uint64_t max_res_bytes = num_distributed * SEQ_REC_HEADER;
			uint64_t max_seqlen = seqlen;
			for (GIntervals::const_iterator iinterv = plain_intervals->begin() + num_probed; iinterv != plain_intervals->end(); ++iinterv) {
				max_res_bytes += (uint64_t)iinterv->range();
				max_seqlen += (uint64_t)iinterv->range();
			}

			// Check the same quantity the sequential path checks, up front. Otherwise the per-record
			// headers push an over-budget call past the limit only once the kids have already done
			// all the I/O, and it reports the generic shared-memory message instead of the one
			// naming the sizes and gmax.data.size.
			iu.verify_max_data_size(max_seqlen, "Result sequence");

			// Shared memory big enough for the whole result may not be available even though the
			// sequential path would have coped; fall through to it rather than failing outright,
			// the same way gscreen/gextract retry without multitasking.
			bool do_child;
			try {
				do_child = iu.distribute_task(0, 1, rdb::MT_MODE_MMAP, max_res_bytes);
			} catch (TGLException &e) {
				// Only safe while nothing has been forked; the mmap happens before the first
				// launch_process(), so get_num_kids() is still 0. If that ever stops holding,
				// re-reading everything serially would race the orphaned children.
				if (!strstr(e.msg(), "Failed to allocate shared memory") || get_num_kids())
					throw;
				num_kids = 0;
				do_child = false;
			}

			if (num_kids > 1 && do_child) { // child
				GIntervalsFetcher1D *kid_intervals = iu.get_kid_intervals1d();
				uint64_t kid_num_intervs = kid_intervals->size();
				uint64_t num_read = 0;
				vector<char> packed;

				// The seqfetch above was opened before the fork, and forked processes SHARE the
				// file offset of an inherited descriptor -- concurrent seeks would hand a kid
				// another kid's bytes. Open our own.
				GenomeSeqFetch kid_seqfetch;
				kid_seqfetch.set_seqdir(seqdir);

				for (kid_intervals->begin_iter(); !kid_intervals->isend(); kid_intervals->next()) {
					const GInterval &interval = kid_intervals->cur_interval();

					kid_seqfetch.read_interval(interval, iu.get_chromkey(), buf);

					uint64_t idx = iu.get_orig_interv_idx(interval);
					uint64_t len = buf.size();
					size_t offset = packed.size();

					packed.resize(offset + SEQ_REC_HEADER + len);
					memcpy(&packed[offset], &idx, sizeof(idx));
					memcpy(&packed[offset + sizeof(idx)], &len, sizeof(len));
					if (len)
						memcpy(&packed[offset + SEQ_REC_HEADER], &buf[0], len);

					if (++num_read % 128 == 0)
						update_progress((unsigned char)(100 * num_read / max<uint64_t>(1, kid_num_intervs)));

					check_interrupt();
				}

				void *res = allocate_res(packed.size());
				if (!packed.empty())
					memcpy(res, &packed[0], packed.size());

				rreturn(R_NilValue);
			}
		}

		if (num_kids > 1) {
			// parent: reassemble the kids' records into the original row order
			uint64_t num_distributed = num_intervs - num_probed;
			uint64_t num_unpacked = 0;

			for (int ikid = 0; ikid < get_num_kids(); ++ikid) {
				const char *ptr = (const char *)get_kid_res(ikid);
				const char *end = ptr + get_kid_res_size(ikid);

				while (ptr < end) {
					uint64_t idx, len;

					if ((uint64_t)(end - ptr) < SEQ_REC_HEADER)
						verror("Truncated sequence record from process %d", ikid);
					memcpy(&idx, ptr, sizeof(idx));
					memcpy(&len, ptr + sizeof(idx), sizeof(len));
					ptr += SEQ_REC_HEADER;

					if (idx >= num_intervs || (uint64_t)(end - ptr) < len)
						verror("Corrupted sequence record from process %d", ikid);

					SET_STRING_ELT(answer, idx, Rf_mkCharLenCE(ptr, (int)len, CE_NATIVE));
					ptr += len;
					++num_unpacked;
				}
			}

			if (num_unpacked != num_distributed)
				verror("Got sequences for %llu intervals out of %llu",
					   (unsigned long long)num_unpacked, (unsigned long long)num_distributed);

			return answer;
		}

		uint64_t iinterv_idx = 0;

		for (intervals->begin_iter(); !intervals->isend(); intervals->next()) {
			if (iinterv_idx++ < num_probed)   // already read (and kept) by the probe above
				continue;

			seqfetch.read_interval(intervals->cur_interval(), iu.get_chromkey(), buf);
			seqlen += buf.size();
			iu.verify_max_data_size(seqlen, "Result sequence");
			// Avoid strlen by constructing string with known length
			SET_STRING_ELT(answer,
						   iu.get_orig_interv_idx(intervals->cur_interval()),
						   Rf_mkCharLenCE(buf.empty() ? "" : &buf[0], (int)buf.size(), CE_NATIVE));
			check_interrupt();
		}
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

// Validate genome index file (called during gdb.init)
// Throws error if index is corrupt or has checksum mismatch
SEXP gseq_validate_index(SEXP _seqdir, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_seqdir) || Rf_length(_seqdir) != 1)
			verror("seqdir argument must be a string");

		const char *seqdir = CHAR(STRING_ELT(_seqdir, 0));
		string idx_path = string(seqdir) + "/genome.idx";

		// Try to load index - will throw TGLException if corrupt
		GenomeIndex index;
		index.load(idx_path);

		return R_NilValue;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
