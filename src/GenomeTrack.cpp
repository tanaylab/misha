/*
 * GenomeTrack.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: hoichman
 */

#include <algorithm>
#include <cstdint>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_set>

#include "GenomeTrack.h"
#include "TrackIndex.h"
#include "TrackIndex2D.h"
#include "TGLException.h"

const char *GenomeTrack::TYPE_NAMES[GenomeTrack::NUM_TYPES] = { "dense", "sparse", "array", "rectangles", "points", "computed", "obsolete rectangles", "obsolete rectangles", "obsolete computed", "obsolete computed", "obsolete computed" };
const bool  GenomeTrack::IS_1D_TRACK[NUM_TYPES] =             { true,    true,     true,     false,        false,    false,      false,                 false,                 false,               false,               false };
const int   GenomeTrack::FORMAT_SIGNATURES[NUM_TYPES] =       { 0,       -1,       -8,       -9,           -10,      -11,        -4,                    -6,                    -3,                  -5,                  -7 };

double (*GenomeTrack::s_rnd_func)() = drand48;

// Static members for track index cache
std::map<std::string, std::shared_ptr<TrackIndex>> GenomeTrack::s_index_cache;
std::mutex GenomeTrack::s_cache_mutex;

std::shared_ptr<TrackIndex> GenomeTrack::get_track_index(const std::string &track_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);

	// Check if already in cache
	auto it = s_index_cache.find(track_dir);
	if (it != s_index_cache.end()) {
		return it->second;
	}

	// Create new index and try to load
	auto idx = std::make_shared<TrackIndex>();
	std::string idx_path = track_dir + "/track.idx";

	if (!idx->load(idx_path)) {
		// Index file doesn't exist - return nullptr
		return nullptr;
	}

	// Cache and return
	s_index_cache[track_dir] = idx;
	return idx;
}

void GenomeTrack::invalidate_index_cache(const std::string &track_dir) {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	s_index_cache.erase(track_dir);
}

void GenomeTrack::clear_index_cache() {
	std::lock_guard<std::mutex> lock(s_cache_mutex);
	s_index_cache.clear();
}

std::string GenomeTrack::get_track_dir(const std::string &filename) {
	// Extract directory from filename
	// filename is typically like "/path/to/trackdir/chrN"
	size_t pos = filename.find_last_of("/");
	if (pos == std::string::npos) {
		// No path separator - assume current directory
		return ".";
	}
	return filename.substr(0, pos);
}

GenomeTrack::ChromFiles::ChromFiles(const GenomeChromKey &chromkey, const string &dir) : m_chromkey(&chromkey), m_dir(dir)
{
	struct stat st;

	// an indexed track keeps all the chromosomes in track.dat
	if (!stat((dir + "/track.idx").c_str(), &st))
		return;

	DIR *d = opendir(dir.c_str());
	if (!d)
		TGLError<GenomeTrack>(FILE_ERROR, "Accessing directory %s: %s\n", dir.c_str(), strerror(errno));

	static const unordered_set<string> index_files = {
		"track.idx", "track.dat", "intervals.idx", "intervals.dat", "intervals2d.idx", "intervals2d.dat",
		"track.idx.tmp", "track.dat.tmp", "intervals.idx.tmp", "intervals.dat.tmp", "intervals2d.idx.tmp", "intervals2d.dat.tmp"
	};
	struct dirent *dirp;

	// also by the "chr" toggle, which the aliases of over 1000 contigs without "chr" leave out (.compute_chrom_aliases)
	auto find_chromid = [&chromkey](const string &name) {
		int chromid = chromkey.find_chromid(name);
		return chromid >= 0 ? chromid : chromkey.find_chromid(name.compare(0, 3, "chr") ? "chr" + name : name.substr(3));
	};

	while ((dirp = readdir(d))) {
		const string name(dirp->d_name);

		if (name[0] == '.' || index_files.count(name))
			continue;

		bool is_symlink = dirp->d_type == DT_LNK;
		if (is_symlink || dirp->d_type == DT_UNKNOWN) {
			string path = dir + "/" + name;
			if (!lstat(path.c_str(), &st) && S_ISLNK(st.st_mode))
				is_symlink = true;
			if (stat(path.c_str(), &st) || !S_ISREG(st.st_mode))  // a dangling symlink, a directory, ...
				continue;
		} else if (dirp->d_type != DT_REG)
			continue;

		File file{name, find_chromid(name), -1, -1};
		for (uint64_t pos = name.find('-'); pos != string::npos; pos = name.find('-', pos + 1)) {
			int chromid1 = find_chromid(name.substr(0, pos));
			int chromid2 = chromid1 < 0 ? -1 : find_chromid(name.substr(pos + 1));
			if (chromid2 >= 0) {
				file.chromid1 = chromid1;
				file.chromid2 = chromid2;
				break;
			}
		}

		if (!is_symlink || file.chromid >= 0 || file.chromid1 >= 0)
			m_files.push_back(file);
	}
	closedir(d);

	sort(m_files.begin(), m_files.end(), [](const File &f1, const File &f2) { return f1.name < f2.name; });

	for (size_t i = 0; i < m_files.size(); ++i) {
		const File &file = m_files[i];
		if (file.chromid >= 0)
			choose(m_chrom2file, file.chromid, i, file.name == chromkey.id2chrom(file.chromid));
		if (file.chromid1 >= 0)
			choose(m_pair2file, pair_key(file.chromid1, file.chromid2), i, file.name == get_2d_filename(chromkey, file.chromid1, file.chromid2));
	}

	// the error names the conflict of the first chromosome (pair) by chrom id
	const Choice *conflict = NULL;
	int conflict_chromid = 0;
	for (const auto &ichoice : m_chrom2file) {
		if (ichoice.second.conflict != NONE && (!conflict || ichoice.first < conflict_chromid)) {
			conflict = &ichoice.second;
			conflict_chromid = ichoice.first;
		}
	}
	if (conflict)
		m_conflict1d = "Directory " + dir + " has two files of chromosome " + chromkey.id2chrom(conflict_chromid) + ", " +
			m_files[conflict->file].name + " and " + m_files[conflict->conflict].name + ", and none named " +
			chromkey.id2chrom(conflict_chromid) + "; remove one of them";

	conflict = NULL;
	uint64_t conflict_key = 0;
	for (const auto &ichoice : m_pair2file) {
		if (ichoice.second.conflict != NONE && (!conflict || ichoice.first < conflict_key)) {
			conflict = &ichoice.second;
			conflict_key = ichoice.first;
		}
	}
	if (conflict) {
		const File &file = m_files[conflict->file];
		m_conflict2d = "Directory " + dir + " has two files of chromosome pair (" + chromkey.id2chrom(file.chromid1) + ", " +
			chromkey.id2chrom(file.chromid2) + "), " + file.name + " and " + m_files[conflict->conflict].name + ", and none named " +
			get_2d_filename(chromkey, file.chromid1, file.chromid2) + "; remove one of them";
	}
}

template <class Key>
void GenomeTrack::ChromFiles::choose(unordered_map<Key, Choice> &choices, const Key &key, size_t ifile, bool canonical)
{
	auto res = choices.emplace(key, Choice{ifile, canonical, NONE});
	Choice &choice = res.first->second;

	if (res.second || choice.canonical)
		return;

	if (canonical) {
		choice = Choice{ifile, true, NONE};
		return;
	}

	// two names of one file (a symlink, a hard link) are not a conflict
	struct stat st1, st2;
	if (choice.conflict == NONE &&
		(stat((m_dir + "/" + m_files[choice.file].name).c_str(), &st1) || stat((m_dir + "/" + m_files[ifile].name).c_str(), &st2) ||
		 st1.st_dev != st2.st_dev || st1.st_ino != st2.st_ino))
		choice.conflict = ifile;
}

// two files of a chromosome or pair, neither with the canonical name
static void check_conflict(const string &conflict)
{
	if (!conflict.empty())
		TGLError<GenomeTrack>(GenomeTrack::BAD_FORMAT, "%s", conflict.c_str());
}

template <class Key>
const GenomeTrack::ChromFiles::Choice *GenomeTrack::ChromFiles::find(const unordered_map<Key, Choice> &choices, const string &conflict, const Key &key) const
{
	check_conflict(conflict);
	auto ichoice = choices.find(key);
	return ichoice == choices.end() ? NULL : &ichoice->second;
}

const string &GenomeTrack::ChromFiles::filename(int chromid) const
{
	const Choice *choice = find(m_chrom2file, m_conflict1d, chromid);
	return choice ? m_files[choice->file].name : get_1d_filename(*m_chromkey, chromid);
}

string GenomeTrack::ChromFiles::filename(int chromid1, int chromid2) const
{
	const Choice *choice = find(m_pair2file, m_conflict2d, pair_key(chromid1, chromid2));
	return choice ? m_files[choice->file].name : get_2d_filename(*m_chromkey, chromid1, chromid2);
}

bool GenomeTrack::ChromFiles::exists(int chromid) const
{
	return find(m_chrom2file, m_conflict1d, chromid) != NULL;
}

bool GenomeTrack::ChromFiles::exists(int chromid1, int chromid2) const
{
	return find(m_pair2file, m_conflict2d, pair_key(chromid1, chromid2)) != NULL;
}

vector<int> GenomeTrack::ChromFiles::chromids() const
{
	check_conflict(m_conflict1d);
	vector<int> res;
	for (const auto &ichoice : m_chrom2file)
		res.push_back(ichoice.first);
	sort(res.begin(), res.end());
	return res;
}

vector<pair<int, int>> GenomeTrack::ChromFiles::chrom_pairs() const
{
	check_conflict(m_conflict2d);
	vector<uint64_t> keys;
	for (const auto &ichoice : m_pair2file)
		keys.push_back(ichoice.first);
	sort(keys.begin(), keys.end());
	vector<pair<int, int>> res;
	for (uint64_t key : keys)
		res.push_back(make_pair((int)(key >> 32), (int)(uint32_t)key));
	return res;
}

GenomeTrack::Type GenomeTrack::get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types)
{
	if (access(track_dir, F_OK))
		TGLError<GenomeTrack>(FILE_ERROR, "Accessing directory %s: %s\n", track_dir, strerror(errno));

	return get_type(ChromFiles(chromkey, track_dir), return_obsolete_types);
}

GenomeTrack::Type GenomeTrack::get_type(const ChromFiles &files, bool return_obsolete_types)
{
	const char *track_dir = files.dir().c_str();

	// First, try to read from track.idx if it exists
	std::string idx_path = std::string(track_dir) + "/track.idx";
	struct stat st;
	if (stat(idx_path.c_str(), &st) == 0) {
		// Try 1D index first
		try {
			auto idx = get_track_index(track_dir);
			if (idx) {
				// Map MishaTrackType to GenomeTrack::Type
				switch (idx->get_track_type()) {
					case MishaTrackType::DENSE:  return FIXED_BIN;
					case MishaTrackType::SPARSE: return SPARSE;
					case MishaTrackType::ARRAY:  return ARRAYS;
					default:
						// Unknown type, fall through
						break;
				}
			}
		} catch (...) {
			// 1D index failed (likely wrong magic header) -- try 2D index
		}

		// Try 2D index
		try {
			auto idx2d = TrackIndex2D::get_track_index_2d(track_dir);
			if (idx2d && idx2d->is_loaded()) {
				switch (idx2d->get_track_type()) {
					case MishaTrack2DType::RECTS:  return RECTS;
					case MishaTrack2DType::POINTS: return POINTS;
				}
			}
		} catch (...) {
			// Fall through to per-chromosome probing on any index error
		}
	}

	// Fall back to per-chromosome probing (per-chrom files)
	for (const ChromFiles::File &file : files.files()) {
		// not named by a chromosome or a pair (attributes, notes, ...)
		if (file.chromid < 0 && file.chromid1 < 0)
			continue;

		string fullpath = files.dir() + "/" + file.name;
		bool is_2d = file.chromid1 >= 0;

		try {
			Type type = s_read_type(fullpath.c_str());

			// a name can be both a chromosome's (with a dash) and a pair's
			if (file.chromid >= 0 && (!is_2d || type == FIXED_BIN || type == SPARSE || type == ARRAYS)) {
				if (type != FIXED_BIN && type != SPARSE && type != ARRAYS)
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file at %s", track_dir);
				return type;
			}

			// 2D track file
			if (type == OLD_RECTS1 || type == OLD_RECTS2 || type == OLD_COMPUTED1 || type == OLD_COMPUTED2 || type == OLD_COMPUTED3) {
				if (return_obsolete_types)
					return type;
				TGLError<GenomeTrack>(OBSOLETE_FORMAT, "Track file at %s is in obsolete format and requires conversion", track_dir);
			}

			if (type != RECTS && type != POINTS && type != COMPUTED)
				TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file at %s", track_dir);
			return type;
		} catch (TGLException &) {
			// ignore files that cannot be read as track data (attributes, etc.)
			continue;
		}
	}

	// If no valid track files found, treat as empty sparse track
	// This allows for empty tracks created by liftover or other operations
	if (files.files().empty()) {
		return SPARSE;
	}

	TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track at %s", track_dir);
	return NUM_TYPES;
}

void GenomeTrack::read_type(const char *filename, const char *mode)
{
	Type type = s_read_type(m_bfile, filename, mode);

	if (type != m_type)
		TGLError<GenomeTrack>(MISMATCH_FORMAT, "Track file %s is in %s format while expected to be in %s format", filename, TYPE_NAMES[type], TYPE_NAMES[m_type]);
}

GenomeTrack::Type GenomeTrack::s_read_type(const char *filename, const char *mode)
{
	BufferedFile bfile;
	return s_read_type(bfile, filename, mode);
}

GenomeTrack::Type GenomeTrack::s_read_type(BufferedFile &bfile, const char *filename, const char *mode)
{
	if (bfile.open(filename, mode))
		TGLError<GenomeTrack>(FILE_ERROR, "Opening a track file %s: %s", filename, strerror(errno));

	int format_signature;

	if (bfile.read(&format_signature, sizeof(format_signature)) != sizeof(format_signature)) {
		if (bfile.error())
			TGLError<GenomeTrack>(FILE_ERROR, "Reading a track file %s: %s", filename, strerror(errno));
		TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of track file %s", filename);
	}

	if (format_signature > 0)
		return FIXED_BIN;

	for (int type = FIXED_BIN + 1; type < NUM_TYPES; ++type) {
		if (format_signature == FORMAT_SIGNATURES[type])
			return (Type)type;
	}

	TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of genome track file %s", filename);
	return NUM_TYPES;
}

void GenomeTrack::write_type(const char *filename, const char *mode)
{
	mode_t old_umask = umask(07);

	if (m_bfile.open(filename, mode)) {
		umask(old_umask);
		TGLError<GenomeTrack>(FILE_ERROR, "Opening a track file %s: %s", filename, strerror(errno));
	}

	umask(old_umask);

	if (m_bfile.write(&FORMAT_SIGNATURES[m_type], sizeof(FORMAT_SIGNATURES[m_type])) != sizeof(FORMAT_SIGNATURES[m_type])) {
		if (m_bfile.error())
			TGLError<GenomeTrack>(FILE_ERROR, "Failed to write a %s track file %s: %s", TYPE_NAMES[m_type], filename, strerror(errno));
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write a %s track file %s", TYPE_NAMES[m_type], filename);
	}
}

void GenomeTrack::load_attrs(const char *, const char *filename, TrackAttrs &attrs)
{
	BufferedFile bfile;
	int c;
	int idx = 0;
	string name;
	string val;

	attrs.clear();

	if (bfile.open(filename, "rb")) {
		if (errno == ENOENT)   // no file = no attributes
			return; 
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to read attributes file %s: %s", filename, strerror(errno));
	}

	while ((c = bfile.getc()) >= 0) {
		if (c) {
			if (idx) 
				val.push_back((char)c);
			else
				name.push_back((char)c);
		} else {
			if (idx) {
				if (name.empty() || val.empty())
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 

				if (attrs.find(name) != attrs.end()) // duplicated attributes
					TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 

				attrs[name] = val;
				name.clear();
				val.clear();
			}
			idx = 1 - idx;
		}
	}

	if (bfile.error()) 
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to read attributes file %s: %s", filename, strerror(errno));

	if (idx) 
		TGLError<GenomeTrack>(BAD_FORMAT, "Invalid format of attributes file %s", filename); 
}

void GenomeTrack::save_attrs(const char *track, const char *filename, const TrackAttrs &attrs)
{
	bool empty_attrs = true;

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) { 
		if (!iattr->second.empty()) {
			empty_attrs = false;
			break;
		}
	}

	if (empty_attrs) {
		if (unlink(filename) && errno != ENOENT)
			TGLError<GenomeTrack>(FILE_ERROR, "Failed accessing attributes file %s: %s", filename, strerror(errno));
		return;
	}

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) {
		if (iattr->first.empty())
			TGLError<GenomeTrack>(BAD_ATTRS, "Track %s: attribute name is an empty string", track); 
	}

	BufferedFile bfile;

	if (bfile.open(filename, "wb"))
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write attributes file %s: %s", filename, strerror(errno));

	for (TrackAttrs::const_iterator iattr = attrs.begin(); iattr != attrs.end(); ++iattr) {
		if (!iattr->second.empty())  {
			bfile.write(iattr->first.c_str(), iattr->first.length() + 1);
			bfile.write(iattr->second.c_str(), iattr->second.length() + 1);
		}
	}

	if (bfile.error())
		TGLError<GenomeTrack>(FILE_ERROR, "Failed to write attributes file %s: %s", filename, strerror(errno));
}
