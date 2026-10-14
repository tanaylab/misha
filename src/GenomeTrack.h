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
#include <unordered_map>

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

	// The per-chromosome files of a track or intervals set directory, from one listing of it. A file
	// is named by a chromosome (1D) or by two chromosomes joined by "-" (2D), each by its canonical
	// name or by an alias (e.g. "1" for chr1; "1-2", "1-chr2" or "chr1-2" for chr1 and chr2). The
	// file of a chromosome or pair is the one with the canonical name; without one, the one named by
	// aliases. Two different files of a chromosome or pair, neither with the canonical name, are an
	// error when the files are read as 1D (2D). A symlink counts only when it leads to a regular
	// file and is named by a chromosome or pair. An indexed track (track.idx) has no files listed,
	// and the index files of an intervals set are not listed.
	class ChromFiles {
	public:
		struct File {
			string name;
			int    chromid;             // the chromosome it is named by, -1 if none
			int    chromid1, chromid2;  // the chromosome pair it is named by, -1 if none
		};

		ChromFiles(const GenomeChromKey &chromkey, const string &dir);

		const string &dir() const { return m_dir; }

		// sorted by name
		const vector<File> &files() const { return m_files; }

		// the file of a chromosome or pair, or the canonical name when there is none
		const string &filename(int chromid) const;
		string filename(int chromid1, int chromid2) const;
		string path(int chromid) const { return m_dir + "/" + filename(chromid); }
		string path(int chromid1, int chromid2) const { return m_dir + "/" + filename(chromid1, chromid2); }

		bool exists(int chromid) const;
		bool exists(int chromid1, int chromid2) const;

		// the chromosomes and pairs that have a file, sorted
		vector<int> chromids() const;
		vector<pair<int, int>> chrom_pairs() const;

	private:
		// indices to m_files; conflict: another file with no canonical name, NONE if there is none
		struct Choice {
			size_t file;
			bool   canonical;
			size_t conflict;
		};
		static const size_t NONE = (size_t)-1;

		const GenomeChromKey             *m_chromkey;
		string                            m_dir;
		vector<File>                      m_files;
		// hashed: a reader looks up every chromosome pair of its scope, 207,025 for 455 contigs
		unordered_map<int, Choice>        m_chrom2file;
		unordered_map<uint64_t, Choice>   m_pair2file;  // by pair_key()
		string                            m_conflict1d;
		string                            m_conflict2d;

		static uint64_t pair_key(int chromid1, int chromid2) { return (uint64_t)(uint32_t)chromid1 << 32 | (uint32_t)chromid2; }
		template <class Key> void choose(unordered_map<Key, Choice> &choices, const Key &key, size_t ifile, bool canonical);
		template <class Key> const Choice *find(const unordered_map<Key, Choice> &choices, const string &conflict, const Key &key) const;
	};

	static Type get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types = false);
	static Type get_type(const ChromFiles &files, bool return_obsolete_types = false);

	static void load_attrs(const char *track, const char *filename, TrackAttrs &attrs);

	static void save_attrs(const char *track, const char *filename, const TrackAttrs &attrs);

	static const string &get_1d_filename(const GenomeChromKey &chromkey, int chromid) { return chromkey.id2chrom(chromid); }

	static const string get_2d_filename(const GenomeChromKey &chromkey, int chromid1, int chromid2) {
		return chromkey.id2chrom(chromid1) + "-" + chromkey.id2chrom(chromid2);
	}

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
	static std::mutex s_cache_mutex;
};

#endif /* GENOMETRACK_H_ */
