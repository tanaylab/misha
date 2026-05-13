#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <dirent.h>
#include <cstring>
#include <set>

#include "GenomeTrackArrays.h"
#include "GenomeTrackComputed.h"
#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "GIntervalsMeta2D.h"
#include "GTrackIntervalsFetcher.h"
#include "GIntervalsBigSet1D.h"
#include "GIntervalsBigSet2D.h"
#include "TrackIndex.h"
#include "TrackIndex2D.h"
#include "rdbprogress.h"
#include "rdbutils.h"

bool GTrackIntervalsFetcher::isbig(const char *track_name, const IntervUtils &iu)
{
	string path = interv2path(iu.get_env(), track_name);
	SEXP gtracks;

    gtracks = rprotect_ptr(find_in_misha(iu.get_env(), "GTRACKS"));
	for (int itrack = 0; itrack < Rf_length(gtracks); ++itrack) {
		const char *track = CHAR(STRING_ELT(gtracks, itrack));
		if (!strcmp(track_name, track))
			return true;
	}
	return false;
}

void GTrackIntervalsFetcher::init(const char *track_name, const IntervUtils &iu)
{
	m_track_name = track_name;
	m_iu = (IntervUtils *)&iu;
}

void GTrackIntervalsFetcher::create_track_meta(const char *track_name, const IntervUtils &iu)
{
	string trackpath(track2path(iu.get_env(), track_name));
	GenomeTrack::Type track_type = GenomeTrack::get_type(trackpath.c_str(), iu.get_chromkey(), true);

	if (track_type == GenomeTrack::FIXED_BIN)
		verror("Track %s is of type %s which cannot be used as a substitute for an intervals set", track_name, GenomeTrack::TYPE_NAMES[track_type]);

	REprintf("Preparing the track %s to be used as an intervals set...\n", track_name);

	if (track_type == GenomeTrack::SPARSE || track_type == GenomeTrack::ARRAYS) {
		const uint64_t num_chroms = iu.get_chromkey().get_num_chroms();
		vector<GIntervalsMeta1D::ChromStat> chromstats(num_chroms);

		// For indexed tracks, the on-disk track.idx tells us exactly which
		// chromids have data (length > 0). This lets us skip O(N) per-chrom
		// find_existing_1d_filename + access() syscalls on indexed-only tracks
		// (where the per-chrom files don't exist anyway). On a 1M-contig DB
		// this drops ~3-5M syscalls per meta creation to ~N.
		const string idx_path = trackpath + "/track.idx";
		struct stat idx_st;
		const bool indexed = (::stat(idx_path.c_str(), &idx_st) == 0);

		if (indexed) {
			std::shared_ptr<TrackIndex> idx = GenomeTrack::get_track_index(trackpath);
			if (!idx)
				verror("Failed to load track index for %s", trackpath.c_str());

			// Count non-empty contigs for progress reporting.
			const vector<TrackContigEntry> &entries = idx->get_all_entries();
			uint64_t num_nonempty = 0;
			for (const TrackContigEntry &e : entries)
				if (e.length > 0)
					++num_nonempty;

			Progress_reporter progress;
			progress.init(num_nonempty, 1);

			// init_read on an indexed track ignores the per-chrom filename
			// component once it has located track.idx; pass the trackdir
			// prefix so get_track_dir() recovers the right directory.
			const string dummy_filename = trackpath + "/.idx_dummy";

			for (const TrackContigEntry &entry : entries) {
				if (entry.length == 0)
					continue;
				const uint64_t chromid = entry.chrom_id;
				if (chromid >= num_chroms)
					continue;

				if (track_type == GenomeTrack::SPARSE) {
					GenomeTrackSparse track;
					track.init_read(dummy_filename.c_str(), chromid);
					GIntervals intervals(track.get_intervals());
					chromstats[chromid] = GIntervalsBigSet1D::get_chrom_stat(&intervals).second;
				} else { // ARRAYS
					GenomeTrackArrays track;
					track.init_read(dummy_filename.c_str(), chromid);
					GIntervals intervals(track.get_intervals());
					chromstats[chromid] = GIntervalsBigSet1D::get_chrom_stat(&intervals).second;
				}

				progress.report(1);
				check_interrupt();
			}
			progress.report_last();
		} else {
			Progress_reporter progress;
			progress.init(num_chroms, 1);

			for (uint64_t chromid = 0; chromid < num_chroms; chromid++) {
				string resolved = GenomeTrack::find_existing_1d_filename(iu.get_chromkey(), trackpath, chromid);
				string filename(trackpath + "/" + resolved);

				if (access(filename.c_str(), R_OK) && errno == ENOENT) {
					progress.report(1);
					continue;
				}

				if (track_type == GenomeTrack::SPARSE) {
					GenomeTrackSparse track;
					track.init_read(filename.c_str(), chromid);
					GIntervals intervals(track.get_intervals());
					chromstats[chromid] = GIntervalsBigSet1D::get_chrom_stat(&intervals).second;
				} else if (track_type == GenomeTrack::ARRAYS) {
					GenomeTrackArrays track;
					track.init_read(filename.c_str(), chromid);
					GIntervals intervals(track.get_intervals());
					chromstats[chromid] = GIntervalsBigSet1D::get_chrom_stat(&intervals).second;
				}

				progress.report(1);
				check_interrupt();
			}
			progress.report_last();
		}
		GIntervalsMeta1D::save_plain_intervals_meta(trackpath.c_str(), chromstats, iu);
	} else if (track_type == GenomeTrack::RECTS || track_type == GenomeTrack::POINTS || track_type == GenomeTrack::COMPUTED) {
		// Phase 7b: enumerate populated chrom-pairs sparsely.
		//
		// The previous implementation walked all N*N chrom-pairs and called
		// init_read + read_interval for each, which loaded a full qtree per
		// pair. On a 1M-contig DB this was 10^12 attempts (multi-day OOM).
		//
		// Two fast paths now:
		//   1. Indexed 2D track (track.idx present): iterate index entries.
		//   2. Legacy per-pair files: readdir the track directory and parse
		//      "chrom1-chrom2" file names.
		GIntervalsMeta2D::ChromStats2D chromstats;

		// Collect (chromid1, chromid2) pairs to process.
		vector<pair<uint64_t, uint64_t>> pairs;

		const string idx_path = trackpath + "/track.idx";
		struct stat idx_st;
		bool indexed_2d = false;
		if (::stat(idx_path.c_str(), &idx_st) == 0) {
			try {
				auto idx2d = TrackIndex2D::get_track_index_2d(trackpath);
				if (idx2d && idx2d->is_loaded()) {
					indexed_2d = true;
					for (const Track2DPairEntry &e : idx2d->entries()) {
						if (e.length == 0)
							continue;
						if (e.chrom1_id >= iu.get_chromkey().get_num_chroms() ||
						    e.chrom2_id >= iu.get_chromkey().get_num_chroms())
							continue;
						pairs.emplace_back((uint64_t)e.chrom1_id, (uint64_t)e.chrom2_id);
					}
				}
			} catch (...) {
				// fall through to legacy per-pair enumeration
				indexed_2d = false;
			}
		}

		if (!indexed_2d) {
			// Legacy per-pair files: scan the track directory and match
			// "chrom1-chrom2" names against the chromkey. This is O(num files)
			// rather than O(N*N), which is the only tractable choice when N is
			// large.
			DIR *d = opendir(trackpath.c_str());
			if (d) {
				struct dirent *de;
				while ((de = readdir(d)) != NULL) {
					if (de->d_name[0] == '.')
						continue;
					const char *dash = strrchr(de->d_name, '-');
					if (!dash || dash == de->d_name)
						continue;
					string chrom1_str(de->d_name, dash - de->d_name);
					string chrom2_str(dash + 1);
					try {
						int c1 = iu.get_chromkey().chrom2id(chrom1_str);
						int c2 = iu.get_chromkey().chrom2id(chrom2_str);
						pairs.emplace_back((uint64_t)c1, (uint64_t)c2);
					} catch (TGLException &) {
						// not a chrom-pair file (attributes, etc.)
						continue;
					}
				}
				closedir(d);
			}
		}

		Progress_reporter progress;
		progress.init(pairs.size(), 1);

		for (const auto &p : pairs) {
			uint64_t chromid1 = p.first;
			uint64_t chromid2 = p.second;
			string filename(trackpath + "/" + GenomeTrack::get_2d_filename(iu.get_chromkey(), chromid1, chromid2));
			unique_ptr<GenomeTrack2D> track;

			if (track_type == GenomeTrack::RECTS) {
				track = unique_ptr<GenomeTrack2D>(new GenomeTrackRectsRects(iu.get_track_chunk_size(), iu.get_track_num_chunks()));
				((GenomeTrackRectsRects *)track.get())->init_read(filename.c_str(), chromid1, chromid2);
			} else if (track_type == GenomeTrack::POINTS) {
				track = unique_ptr<GenomeTrack2D>(new GenomeTrackRectsPoints(iu.get_track_chunk_size(), iu.get_track_num_chunks()));
				((GenomeTrackRectsPoints *)track.get())->init_read(filename.c_str(), chromid1, chromid2);
			} else if (track_type == GenomeTrack::COMPUTED) {
				track = unique_ptr<GenomeTrack2D>(new GenomeTrackComputed(get_groot(iu.get_env()), iu.get_track_chunk_size(), iu.get_track_num_chunks()));
				((GenomeTrackComputed *)track.get())->init_read(filename.c_str(), chromid1, chromid2);
			}

			track->read_interval(GInterval2D(chromid1, 0, iu.get_chromkey().get_chrom_size(chromid1),
											 chromid2, 0, iu.get_chromkey().get_chrom_size(chromid2)),
								 DiagonalBand());

			GIntervalsMeta2D::ChromStat stat;
			stat.contains_overlaps = false;
			if (track_type == GenomeTrack::RECTS)
				stat.size = ((GenomeTrackRectsRects *)track.get())->get_qtree().get_num_objs();
			else if (track_type == GenomeTrack::POINTS)
				stat.size = ((GenomeTrackRectsPoints *)track.get())->get_qtree().get_num_objs();
			else if (track_type == GenomeTrack::COMPUTED)
				stat.size = ((GenomeTrackComputed *)track.get())->get_qtree().get_num_objs();
			stat.surface = track->last_occupied_area();

			chromstats.set((int)chromid1, (int)chromid2, stat);
			progress.report(1);
			check_interrupt();
		}

		GIntervalsMeta2D::save_plain_intervals_meta(trackpath.c_str(), chromstats, iu);
		progress.report_last();
	}
	REprintf("Track %s is modified and ready to be used as an intervals set.\n", track_name);
}
