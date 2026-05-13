/*
 * GenomeTrackCreate.cpp
 *
 *  Created on: May 2, 2010
 *      Author: hoichman
 */

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "rdbinterval.h"
#include "rdbutils.h"

#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "GenomeTrackIndexedWriter.h"
#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "TrackExpressionScanner.h"
#include "TrackExpressionIterator.h"
#include "TrackExpressionFixedBinIterator.h"
#include "TrackExpressionSparseIterator.h"

namespace {

// Mirror of R's .gdb.is_indexed: an indexed DB has both seq/genome.idx
// and seq/genome.seq under GROOT. We re-check from C++ rather than
// passing a flag from R so the writer's behaviour stays tied to actual
// on-disk state.
bool is_indexed_db(SEXP envir) {
    SEXP groot = find_in_misha(envir, "GROOT");
    if (groot == R_NilValue || groot == R_UnboundValue ||
        TYPEOF(groot) != STRSXP || Rf_length(groot) < 1) {
        return false;
    }
    const std::string root = CHAR(STRING_ELT(groot, 0));
    struct stat st;
    return ::stat((root + "/seq/genome.idx").c_str(), &st) == 0 &&
           ::stat((root + "/seq/genome.seq").c_str(), &st) == 0;
}

} // namespace

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrackcreate(SEXP track, SEXP expr, SEXP _iterator_policy, SEXP _band, SEXP envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(track) || Rf_length(track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(expr) || Rf_length(expr) != 1)
			verror("Track expression argument is not a string");

		const char *track_str = CHAR(STRING_ELT(track, 0));

		string dirname = create_track_dir(envir, track_str);
		IntervUtils iu(envir);
		TrackExprScanner scanner(iu);
		char filename[FILENAME_MAX];
		GIntervals all_genome_intervs1d;
		GIntervals2D all_genome_intervs2d;
		iu.get_all_genome_intervs(all_genome_intervs1d);
		iu.get_all_genome_intervs(all_genome_intervs2d);

		scanner.begin(expr, &all_genome_intervs1d, &all_genome_intervs2d, _iterator_policy, _band);

		TrackExpressionIteratorBase::Type itr_type = scanner.get_iterator()->get_type();

		if (itr_type == TrackExpressionIteratorBase::FIXED_BIN) {
			const unsigned bin_size = ((TrackExpressionFixedBinIterator *)scanner.get_iterator())->get_bin_size();

			if (is_indexed_db(envir)) {
				// Streaming indexed-direct path: write track.dat + track.idx
				// in one pass instead of N per-chrom files that the R-level
				// gtrack.convert_to_indexed step would unlink immediately
				// afterward. Output is bit-for-bit identical to the
				// per-chrom + pack flow (tests/testthat/test-track-indexed-direct.R).
				GenomeTrackIndexedWriter iw;
				iw.init(dirname, GenomeTrack::FIXED_BIN, (uint32_t)iu.get_chromkey().get_num_chroms());

				vector<char> header_bytes;
				GenomeTrackFixedBin::pack_header(header_bytes, bin_size);

				int cur_chromid = -1;
				for (; !scanner.isend(); scanner.next()) {
					int chromid = scanner.last_interval1d().chromid;
					if (cur_chromid != chromid) {
						if (cur_chromid != -1)
							iw.end_chrom();
						cur_chromid = chromid;
						iw.begin_chrom(chromid);
						iw.write_bytes(header_bytes.data(), header_bytes.size());
					}
					float val = scanner.last_real(0);
					iw.write_bytes(&val, sizeof(val));
				}
				if (cur_chromid != -1)
					iw.end_chrom();
				iw.finalize();
			} else {
				int cur_chromid = -1;
				GenomeTrackFixedBin gtrack;

				for (; !scanner.isend(); scanner.next()) {
					if (cur_chromid != scanner.last_interval1d().chromid) {
						cur_chromid = scanner.last_interval1d().chromid;
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), cur_chromid).c_str());
						gtrack.init_write(filename, bin_size, cur_chromid);
					}
					gtrack.write_next_bin(scanner.last_real(0));
				}
			}
		} else if (itr_type == TrackExpressionIteratorBase::INTERVALS1D) {
			if (is_indexed_db(envir)) {
				// Streaming indexed-direct sparse path: write track.dat +
				// track.idx in one pass. Mirrors the FIXED_BIN branch
				// above. Per-chrom byte layout = 4-byte format signature
				// header + 20-byte records (int64 start, int64 end,
				// float val).
				//
				// To preserve bit-for-bit equivalence with the legacy
				// per-chrom + pack flow we emit a header-only block for
				// every genome chrom that the scanner does NOT cover.
				// The legacy INTERVALS1D branch calls init_write() for
				// each chrom in all_genome_intervs1d (writes the 4-byte
				// signature even for chroms with no intervals), and
				// pack then concatenates those bytes into track.dat
				// with length=4 entries.
				//
				// Drive in chromid order: build the ordered chromid
				// list from all_genome_intervs1d, then for each
				// chromid C emit a header and drain any scanner
				// records belonging to C. This keeps track.dat in
				// strict genome order and satisfies begin_chrom's
				// monotonicity check.
				GenomeTrackIndexedWriter iw;
				iw.init(dirname, GenomeTrack::SPARSE, (uint32_t)iu.get_chromkey().get_num_chroms());

				vector<char> header_bytes;
				GenomeTrackSparse::pack_header(header_bytes);

				vector<int> all_chromids;
				all_chromids.reserve(all_genome_intervs1d.size());
				for (GIntervals::const_iterator it = all_genome_intervs1d.begin(); it != all_genome_intervs1d.end(); ++it)
					all_chromids.push_back(it->chromid);
				std::sort(all_chromids.begin(), all_chromids.end());
				all_chromids.erase(std::unique(all_chromids.begin(), all_chromids.end()), all_chromids.end());

				vector<char> rec_buf;
				for (size_t i = 0; i < all_chromids.size(); ++i) {
					int chromid = all_chromids[i];
					iw.begin_chrom(chromid);
					iw.write_bytes(header_bytes.data(), header_bytes.size());

					// Drain records belonging to this chrom. The scanner
					// advances chroms in genome order, so once we see a
					// record with chromid != C we leave the rest for the
					// later iteration.
					while (!scanner.isend() && scanner.last_interval1d().chromid == chromid) {
						const GInterval &iv = scanner.last_interval1d();
						rec_buf.clear();
						GenomeTrackSparse::pack_record(rec_buf, iv.start, iv.end, scanner.last_real(0));
						iw.write_bytes(rec_buf.data(), rec_buf.size());
						scanner.next();
					}
					iw.end_chrom();
				}
				iw.finalize();
			} else {
				int cur_chromid = -1;
				GenomeTrackSparse gtrack;
				set<int> created_chromids;

				for (; !scanner.isend(); scanner.next()) {
					if (cur_chromid != scanner.last_interval1d().chromid) {
						cur_chromid = scanner.last_interval1d().chromid;
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), cur_chromid).c_str());
						gtrack.init_write(filename, cur_chromid);
						created_chromids.insert(cur_chromid);
					}
					gtrack.write_next_interval(scanner.last_interval1d(), scanner.last_real(0));
				}

				// some of the chromosome could be previously skipped; we still must create them even if they are empty
				for (GIntervals::const_iterator iinterv = all_genome_intervs1d.begin(); iinterv != all_genome_intervs1d.end(); ++iinterv) {
					if (created_chromids.find(iinterv->chromid) == created_chromids.end()) {
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), iinterv->chromid).c_str());
						gtrack.init_write(filename, iinterv->chromid);
					}
				}
			}
		} else if (itr_type == TrackExpressionIteratorBase::INTERVALS2D) {
			int cur_chromid1 = -1;
			int cur_chromid2 = -1;
			GenomeTrackRectsRects gtrack(iu.get_track_chunk_size(), iu.get_track_num_chunks());
			RectsQuadTree qtree;

			for (; !scanner.isend(); scanner.next()) {
				const GInterval2D &interv = scanner.last_interval2d();

				if (cur_chromid1 != interv.chromid1() || cur_chromid2 != interv.chromid2()) {
					if (gtrack.opened())
						gtrack.write(qtree);

					cur_chromid1 = interv.chromid1();
					cur_chromid2 = interv.chromid2();
					snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_2d_filename(iu.get_chromkey(), cur_chromid1, cur_chromid2).c_str());

					qtree.reset(0, 0, iu.get_chromkey().get_chrom_size(cur_chromid1), iu.get_chromkey().get_chrom_size(cur_chromid2));
					gtrack.init_write(filename, cur_chromid1, cur_chromid2);
				}

				qtree.insert(RectsQuadTree::ValueType(interv, scanner.last_real(0)));
			}

			if (gtrack.opened())
				gtrack.write(qtree);
		} else {
			if (itr_type >= 0 && itr_type < sizeof(TrackExpressionIteratorBase::TYPE_NAMES) / sizeof(TrackExpressionIteratorBase::TYPE_NAMES[0])) {
    			verror("Iterator type %s is not supported by the function", TrackExpressionIteratorBase::TYPE_NAMES[itr_type]);
			} else {
    			verror("Invalid iterator type encountered");
			}
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	return R_NilValue;
}


SEXP gtrackcreate_multitask(SEXP track, SEXP expr, SEXP _iterator_policy, SEXP _band, SEXP envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(track) || Rf_length(track) != 1)
			verror("Track argument is not a string");

		if (!Rf_isString(expr) || Rf_length(expr) != 1)
			verror("Track expression argument is not a string");

		const char *track_str = CHAR(STRING_ELT(track, 0));

		// On indexed DBs, fixed-bin tracks must be written in a single
		// process: GenomeTrackIndexedWriter appends to a single track.dat,
		// so concurrent children would corrupt offsets. The R-level
		// caller routes through gtrackcreate when this returns; here we
		// short-circuit by falling back to the single-process path
		// inline. We still honour the multitasking entrypoint contract
		// (must rreturn before falling off the end).
		//
		// We can't easily distinguish "fixed-bin" before the scanner
		// starts (the iterator type is derived from the expression),
		// so this branch handles fixed-bin specifically inside the
		// distribute_task body, mirroring gtrackcreate. The parallelism
		// loss is small relative to the FS-op savings of avoiding the
		// per-chrom + pack roundtrip.
		const bool indexed = is_indexed_db(envir);

		string dirname = create_track_dir(envir, track_str);
		IntervUtils iu(envir);
		char filename[FILENAME_MAX];
		GIntervals all_genome_intervs1d;
		GIntervals2D all_genome_intervs2d;
		iu.get_all_genome_intervs(all_genome_intervs1d);
		iu.get_all_genome_intervs(all_genome_intervs2d);

		if (indexed) {
			// Single-process path. Same body as the indexed branch in
			// gtrackcreate. We can't dispatch on iterator type until the
			// scanner is constructed, so this path handles all iterator
			// types (FIXED_BIN streams direct; others fall through to
			// the per-chrom path that pack_per_chrom_to_indexed converts
			// at the R level - unchanged behaviour for those types in
			// Phase 4b).
			TrackExprScanner scanner(iu);
			scanner.begin(expr, &all_genome_intervs1d, &all_genome_intervs2d, _iterator_policy, _band);
			TrackExpressionIteratorBase::Type itr_type = scanner.get_iterator()->get_type();

			if (itr_type == TrackExpressionIteratorBase::FIXED_BIN) {
				const unsigned bin_size = ((TrackExpressionFixedBinIterator *)scanner.get_iterator())->get_bin_size();
				GenomeTrackIndexedWriter iw;
				iw.init(dirname, GenomeTrack::FIXED_BIN, (uint32_t)iu.get_chromkey().get_num_chroms());

				vector<char> header_bytes;
				GenomeTrackFixedBin::pack_header(header_bytes, bin_size);

				int cur_chromid = -1;
				for (; !scanner.isend(); scanner.next()) {
					int chromid = scanner.last_interval1d().chromid;
					if (cur_chromid != chromid) {
						if (cur_chromid != -1)
							iw.end_chrom();
						cur_chromid = chromid;
						iw.begin_chrom(chromid);
						iw.write_bytes(header_bytes.data(), header_bytes.size());
					}
					float val = scanner.last_real(0);
					iw.write_bytes(&val, sizeof(val));
				}
				if (cur_chromid != -1)
					iw.end_chrom();
				iw.finalize();
				rreturn(R_NilValue);
			}

			// INTERVALS1D on indexed DBs: stream the sparse data
			// directly into track.dat + track.idx in the parent
			// process. Same write pattern as gtrackcreate's indexed
			// INTERVALS1D branch (see above for the byte-layout and
			// header-for-skipped-chroms rationale). The single-file
			// append model of GenomeTrackIndexedWriter precludes
			// parallel children, so we forgo the multitasking
			// parallelism here, same trade-off as fixed-bin.
			if (itr_type == TrackExpressionIteratorBase::INTERVALS1D) {
				GenomeTrackIndexedWriter iw;
				iw.init(dirname, GenomeTrack::SPARSE, (uint32_t)iu.get_chromkey().get_num_chroms());

				vector<char> header_bytes;
				GenomeTrackSparse::pack_header(header_bytes);

				vector<int> all_chromids;
				all_chromids.reserve(all_genome_intervs1d.size());
				for (GIntervals::const_iterator it = all_genome_intervs1d.begin(); it != all_genome_intervs1d.end(); ++it)
					all_chromids.push_back(it->chromid);
				std::sort(all_chromids.begin(), all_chromids.end());
				all_chromids.erase(std::unique(all_chromids.begin(), all_chromids.end()), all_chromids.end());

				vector<char> rec_buf;
				for (size_t i = 0; i < all_chromids.size(); ++i) {
					int chromid = all_chromids[i];
					iw.begin_chrom(chromid);
					iw.write_bytes(header_bytes.data(), header_bytes.size());

					while (!scanner.isend() && scanner.last_interval1d().chromid == chromid) {
						const GInterval &iv = scanner.last_interval1d();
						rec_buf.clear();
						GenomeTrackSparse::pack_record(rec_buf, iv.start, iv.end, scanner.last_real(0));
						iw.write_bytes(rec_buf.data(), rec_buf.size());
						scanner.next();
					}
					iw.end_chrom();
				}
				iw.finalize();
				rreturn(R_NilValue);
			}
			if (itr_type == TrackExpressionIteratorBase::INTERVALS2D) {
				int cur_chromid1 = -1;
				int cur_chromid2 = -1;
				GenomeTrackRectsRects gtrack(iu.get_track_chunk_size(), iu.get_track_num_chunks());
				RectsQuadTree qtree;

				for (; !scanner.isend(); scanner.next()) {
					const GInterval2D &interv = scanner.last_interval2d();

					if (cur_chromid1 != interv.chromid1() || cur_chromid2 != interv.chromid2()) {
						if (gtrack.opened())
							gtrack.write(qtree);

						cur_chromid1 = interv.chromid1();
						cur_chromid2 = interv.chromid2();
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_2d_filename(iu.get_chromkey(), cur_chromid1, cur_chromid2).c_str());

						qtree.reset(0, 0, iu.get_chromkey().get_chrom_size(cur_chromid1), iu.get_chromkey().get_chrom_size(cur_chromid2));
						gtrack.init_write(filename, cur_chromid1, cur_chromid2);
					}

					qtree.insert(RectsQuadTree::ValueType(interv, scanner.last_real(0)));
				}

				if (gtrack.opened())
					gtrack.write(qtree);
				rreturn(R_NilValue);
			}
			if (itr_type >= 0 && itr_type < sizeof(TrackExpressionIteratorBase::TYPE_NAMES) / sizeof(TrackExpressionIteratorBase::TYPE_NAMES[0])) {
				verror("Iterator type %s is not supported by the function", TrackExpressionIteratorBase::TYPE_NAMES[itr_type]);
			} else {
				verror("Invalid iterator type encountered");
			}
		}

		if (!iu.prepare4multitasking(expr, &all_genome_intervs1d, &all_genome_intervs2d, _iterator_policy, _band))
			rreturn(R_NilValue);

		if (iu.distribute_task(0, 0)) {  // child process
			TrackExprScanner scanner(iu);

			scanner.begin(expr, iu.get_kid_intervals1d(), iu.get_kid_intervals2d(), _iterator_policy, _band);

			TrackExpressionIteratorBase::Type itr_type = scanner.get_iterator()->get_type();

			if (itr_type == TrackExpressionIteratorBase::FIXED_BIN) {
				int cur_chromid = -1;
				GenomeTrackFixedBin gtrack;

				for (; !scanner.isend(); scanner.next()) {
                    if (cur_chromid != scanner.last_interval1d().chromid) {
                        cur_chromid = scanner.last_interval1d().chromid;
                        snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), cur_chromid).c_str());
                        gtrack.init_write(filename, ((TrackExpressionFixedBinIterator *)scanner.get_iterator())->get_bin_size(), cur_chromid);
                    }
                    gtrack.write_next_bin(scanner.last_real(0));
				}
			} else if (itr_type == TrackExpressionIteratorBase::INTERVALS1D) {
				int cur_chromid = -1;
				GenomeTrackSparse gtrack;
				set<int> created_chromids;

				for (; !scanner.isend(); scanner.next()) {
					if (cur_chromid != scanner.last_interval1d().chromid) {
						cur_chromid = scanner.last_interval1d().chromid;
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), cur_chromid).c_str());
						gtrack.init_write(filename, cur_chromid);
						created_chromids.insert(cur_chromid);
					}
					gtrack.write_next_interval(scanner.last_interval1d(), scanner.last_real(0));
				}

				// some of the chromosome could be previously skipped; we still must create them even if they are empty
				GIntervals *kid_intervals1d = (GIntervals *)iu.get_kid_intervals1d();
				for (GIntervals::const_iterator iinterv = kid_intervals1d->begin(); iinterv != kid_intervals1d->end(); ++iinterv) {
					if (created_chromids.find(iinterv->chromid) == created_chromids.end()) {
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), iinterv->chromid).c_str());
						gtrack.init_write(filename, iinterv->chromid);
					}
				}
			} else if (itr_type == TrackExpressionIteratorBase::INTERVALS2D) {
				int cur_chromid1 = -1;
				int cur_chromid2 = -1;
				GenomeTrackRectsRects gtrack(iu.get_track_chunk_size(), iu.get_track_num_chunks());
				RectsQuadTree qtree;

				for (; !scanner.isend(); scanner.next()) {
					const GInterval2D &interv = scanner.last_interval2d();

					if (cur_chromid1 != interv.chromid1() || cur_chromid2 != interv.chromid2()) {
						if (gtrack.opened())
							gtrack.write(qtree);

						cur_chromid1 = interv.chromid1();
						cur_chromid2 = interv.chromid2();
						snprintf(filename, sizeof(filename), "%s/%s", dirname.c_str(), GenomeTrack::get_2d_filename(iu.get_chromkey(), cur_chromid1, cur_chromid2).c_str());

						qtree.reset(0, 0, iu.get_chromkey().get_chrom_size(cur_chromid1), iu.get_chromkey().get_chrom_size(cur_chromid2));
						gtrack.init_write(filename, cur_chromid1, cur_chromid2);
					}

					qtree.insert(RectsQuadTree::ValueType(interv, scanner.last_real(0)));
				}

				if (gtrack.opened())
					gtrack.write(qtree);
			} else {
				if (itr_type >= 0 && itr_type < sizeof(TrackExpressionIteratorBase::TYPE_NAMES) / sizeof(TrackExpressionIteratorBase::TYPE_NAMES[0])) {
    				verror("Iterator type %s is not supported by the function", TrackExpressionIteratorBase::TYPE_NAMES[itr_type]);
				} else {
    				verror("Invalid iterator type encountered");
				}
			}				
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

	rreturn(R_NilValue);
}

}
