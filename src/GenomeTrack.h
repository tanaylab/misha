/*
 * GenomeTrack.h
 *
 *  Created on: Mar 10, 2010
 *      Author: hoichman
 */

#ifndef GENOMETRACK_H_
#define GENOMETRACK_H_

#include <map>
#include <memory>
#include <mutex>
#include <string>

#include "BufferedFile.h"
#include "GenomeChromKey.h"

using namespace std;

// Forward declaration
class TrackIndex;

// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeTrack {
public:
	typedef map<string, string> TrackAttrs;

	enum Type { FIXED_BIN, SPARSE, ARRAYS, RECTS, POINTS, COMPUTED, OLD_RECTS1, OLD_RECTS2, OLD_COMPUTED1, OLD_COMPUTED2, OLD_COMPUTED3, NUM_TYPES };
	enum Error { NOT_2D, BAD_FORMAT, OBSOLETE_FORMAT, MISMATCH_FORMAT, FILE_ERROR, BAD_ATTRS };

	static const char *TYPE_NAMES[NUM_TYPES];
	static const int   FORMAT_SIGNATURES[NUM_TYPES];

	virtual ~GenomeTrack() {}

	Type get_type() const { return m_type; }

	bool is_1d() const { return is_1d(m_type); }
	bool is_2d() const { return is_2d(m_type); }

	static bool is_1d(Type type) { return IS_1D_TRACK[type]; }
	static bool is_2d(Type type) { return !IS_1D_TRACK[type]; }

	// returns true if the track file is opened
	bool opened() const { return m_bfile.opened(); }

	const string &file_name() const { return m_bfile.file_name(); }

    static void set_rnd_func(double (*rnd_func)()) { s_rnd_func = rnd_func; }

	static Type get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types = false);

	static void load_attrs(const char *track, const char *filename, TrackAttrs &attrs);

	static void save_attrs(const char *track, const char *filename, const TrackAttrs &attrs);

	static const string &get_1d_filename(const GenomeChromKey &chromkey, int chromid) { return chromkey.id2chrom(chromid); }
	static string find_existing_1d_filename(const GenomeChromKey &chromkey, const string &track_dir, int chromid);

	static const string get_2d_filename(const GenomeChromKey &chromkey, int chromid1, int chromid2) {
		return chromkey.id2chrom(chromid1) + "-" + chromkey.id2chrom(chromid2);
	}

	// Per-pair files of a 2D track named by aliases of their chromosomes (e.g. "1-2" or "chr1-2"
	// for chr1 and chr2), keyed by chrom pair
	typedef map<pair<int, int>, string> Pair2Filename;

	// The per-pair files in track_dir whose names resolve (get_chromid_2d) to a chrom pair that
	// has no file named get_2d_filename; when several do, the first in sorted order. This is
	// the 2D counterpart of find_existing_1d_filename. Empty for an indexed track.
	static void get_2d_alias_filenames(const GenomeChromKey &chromkey, const string &track_dir, Pair2Filename &alias_filenames);

	// The file of (chromid1, chromid2) in alias_filenames if it is there, otherwise get_2d_filename
	static string get_2d_filename(const GenomeChromKey &chromkey, int chromid1, int chromid2, const Pair2Filename &alias_filenames) {
		if (!alias_filenames.empty()) {
			Pair2Filename::const_iterator ifilename = alias_filenames.find(pair<int, int>(chromid1, chromid2));
			if (ifilename != alias_filenames.end())
				return ifilename->second;
		}
		return get_2d_filename(chromkey, chromid1, chromid2);
	}

	// The same for a caller that resolves one pair per call: the alias filenames of track_dir are
	// listed once and cached per directory, and listed again when the directory's modification
	// time changes or the cached file of the pair does not exist
	static string get_2d_filename(const GenomeChromKey &chromkey, const string &track_dir, int chromid1, int chromid2);

	static const int get_chromid_1d(const GenomeChromKey &chromkey, const string &filename) { return chromkey.chrom2id(filename); }

	static const pair<int, int> get_chromid_2d(const GenomeChromKey &chromkey, const string &filename);

//protected:
public:
	static const bool IS_1D_TRACK[NUM_TYPES];

    static double (*s_rnd_func)();

	BufferedFile m_bfile;
	Type         m_type;

	GenomeTrack(Type type) : m_type(type) {}

	void read_type(const char *filename, const char *mode = "rb");

	void write_type(const char *filename, const char *mode = "wb");

	static Type s_read_type(const char *filename, const char *mode = "rb");

	static Type s_read_type(BufferedFile &bfile, const char *filename, const char *mode = "rb");

	// Helper to get-or-load the track index (cached, thread-safe).
	// Returns nullptr if track.idx is not present in track_dir.
	static std::shared_ptr<TrackIndex> get_track_index(const std::string &track_dir);

	// Drop the cached TrackIndex for `track_dir`. Call this whenever the
	// on-disk track contents at `track_dir` have changed (rm, create,
	// convert) so subsequent get_track_index() calls re-read track.idx.
	// Safe to call when no entry is cached.
	static void invalidate_index_cache(const std::string &track_dir);

	// Wipe every cached TrackIndex entry. Used by tests and as a coarse
	// reset when callers cannot enumerate the affected dirs.
	static void clear_index_cache();

	// Helper to extract track directory from a per-chrom filename.
	static std::string get_track_dir(const std::string &filename);

protected:
	// Track index cache (thread-safe)
	static std::map<std::string, std::shared_ptr<TrackIndex>> s_index_cache;
	// get_2d_alias_filenames per directory, with the directory's mtime in ns; guarded by s_cache_mutex too
	static std::map<std::string, std::pair<int64_t, Pair2Filename>> s_alias_filenames_cache;
	static std::mutex s_cache_mutex;
};

#endif /* GENOMETRACK_H_ */
