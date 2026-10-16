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
#include <fstream>
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

GenomeTrack::ChromFiles::ChromFiles(const GenomeChromKey &chromkey, const string &dir, bool converting) : m_chromkey(&chromkey), m_dir(dir)
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

	// the name with "chr" added or removed
	auto toggled = [](const string &name) { return name.compare(0, 3, "chr") ? "chr" + name : name.substr(3); };

	// also by the "chr" toggle, which the aliases of over 1000 contigs without "chr" leave out (.compute_chrom_aliases)
	auto find_chromid = [&chromkey, &toggled](const string &name) {
		int chromid = chromkey.find_chromid(name);
		return chromid >= 0 ? chromid : chromkey.find_chromid(toggled(name));
	};

	unordered_map<string, string> unreachable;

	// a listing cut short by an I/O error (EIO, ESTALE on NFS) would silently miss files
	for (;;) {
		errno = 0;
		if (!(dirp = readdir(d))) {
			if (errno) {
				int err = errno;
				closedir(d);
				TGLError<GenomeTrack>(FILE_ERROR, "Reading directory %s: %s\n", dir.c_str(), strerror(err));
			}
			break;
		}
		const string name(dirp->d_name);

		if (name[0] == '.' || index_files.count(name))
			continue;

		if (dirp->d_type != DT_LNK && dirp->d_type != DT_UNKNOWN && dirp->d_type != DT_REG)
			continue;

		File file{name, find_chromid(name), -1, -1, dirp->d_type == DT_LNK};
		for (uint64_t pos = name.find('-'); pos != string::npos; pos = name.find('-', pos + 1)) {
			int chromid1 = find_chromid(name.substr(0, pos));
			int chromid2 = chromid1 < 0 ? -1 : find_chromid(name.substr(pos + 1));
			if (chromid2 >= 0) {
				file.chromid1 = chromid1;
				file.chromid2 = chromid2;
				break;
			}
		}

		// leftovers of interrupted writes (*.tmp, which master's 2D converter skipped), unless named by a chromosome
		if (file.chromid < 0 && file.chromid1 < 0 && name.size() > 4 && !name.compare(name.size() - 4, 4, ".tmp"))
			continue;

		bool &is_symlink = file.link;
		if (is_symlink || dirp->d_type == DT_UNKNOWN) {
			string path = dir + "/" + name;
			if (!lstat(path.c_str(), &st) && S_ISLNK(st.st_mode))
				is_symlink = true;
			if (stat(path.c_str(), &st)) {
				// a dangling symlink is skipped, as is any file not named by a chromosome; one named by a chromosome
				// or pair that cannot be reached (EACCES, ELOOP, an I/O error) is an error where it is the file
				// read, below
				if (errno != ENOENT && (file.chromid >= 0 || file.chromid1 >= 0)) {
					unreachable.emplace(name, path + ": " + strerror(errno));
					m_files.push_back(file);
				}
				continue;
			}
			if (!S_ISREG(st.st_mode))  // a directory, ...
				continue;
		}

		if (!is_symlink || file.chromid >= 0 || file.chromid1 >= 0)
			m_files.push_back(file);
	}
	closedir(d);

	sort(m_files.begin(), m_files.end(), [](const File &f1, const File &f2) { return f1.name < f2.name; });

	// Next to an intervals set's index only a canonically named file holds data newer than it: no misha writes a
	// file named by an alias there, and one left by an earlier conversion is older than the index. A conversion
	// packs them all, as master did.
	// (an index that cannot be stat'ed for a reason other than its absence counts as there, so readers report it)
	auto has_index = [&st](const string &path) { return !stat(path.c_str(), &st) || errno != ENOENT; };
	const bool indexed1d = !converting && has_index(dir + "/intervals.idx");
	const bool indexed2d = !converting && has_index(dir + "/intervals2d.idx");

	for (size_t i = 0; i < m_files.size(); ++i) {
		const File &file = m_files[i];
		if (file.chromid >= 0) {
			const string &name = chromkey.id2chrom(file.chromid);
			int rank = file.name == name ? 0 : file.name == toggled(name) ? 1 : 2;
			if (!rank || !indexed1d)
				choose(m_chrom2file, file.chromid, i, rank);
		}
		if (file.chromid1 >= 0) {
			const string &name1 = chromkey.id2chrom(file.chromid1);
			const string &name2 = chromkey.id2chrom(file.chromid2);
			int rank = file.name == name1 + "-" + name2 ? 0 : file.name == toggled(name1) + "-" + toggled(name2) ? 1 : 2;
			if (!rank || !indexed2d)
				choose(m_pair2file, pair_key(file.chromid1, file.chromid2), i, rank);
		}
	}

	// a file read that cannot be reached is an error when the files are read as 1D (2D), as a conflict
	for (const auto &ichoice : m_chrom2file) {
		if (unreachable.empty())
			break;
		auto iunreachable = unreachable.find(m_files[ichoice.second.file].name);
		if (iunreachable != unreachable.end() && m_conflict1d.empty())
			m_conflict1d = "Accessing " + iunreachable->second;
	}
	for (const auto &ichoice : m_pair2file) {
		if (unreachable.empty())
			break;
		auto iunreachable = unreachable.find(m_files[ichoice.second.file].name);
		if (iunreachable != unreachable.end() && m_conflict2d.empty())
			m_conflict2d = "Accessing " + iunreachable->second;
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
			chromkey.id2chrom(conflict_chromid) + " or " + toggled(chromkey.id2chrom(conflict_chromid)) + "; remove one of them";

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
			get_2d_filename(chromkey, file.chromid1, file.chromid2) + " or " + toggled(chromkey.id2chrom(file.chromid1)) + "-" +
			toggled(chromkey.id2chrom(file.chromid2)) + "; remove one of them";
	}
}

// whether two files hold the same bytes
static bool same_contents(const string &path1, const string &path2, off_t size1, off_t size2)
{
	if (size1 != size2)
		return false;
	ifstream f1(path1, ios::binary), f2(path2, ios::binary);
	vector<char> buf1(1 << 20), buf2(1 << 20);
	while (f1 && f2) {
		f1.read(buf1.data(), buf1.size());
		f2.read(buf2.data(), buf2.size());
		if (f1.gcount() != f2.gcount() || memcmp(buf1.data(), buf2.data(), f1.gcount()))
			return false;
	}
	return f1.eof() && f2.eof();
}

template <class Key>
void GenomeTrack::ChromFiles::choose(unordered_map<Key, Choice> &choices, const Key &key, size_t ifile, int rank)
{
	auto res = choices.emplace(key, Choice{ifile, rank, NONE});
	Choice &choice = res.first->second;

	if (res.second || rank > choice.rank)
		return;

	if (rank < choice.rank) {
		choice = Choice{ifile, rank, NONE};
		return;
	}

	// two files of one rank (not the canonical one, a single name): two names of one file (a symlink, a hard link) are not a conflict, nor two files with the same bytes (a copy
	// of such names, which makes a file of each)
	const string path1 = m_dir + "/" + m_files[choice.file].name;
	const string path2 = m_dir + "/" + m_files[ifile].name;
	struct stat st1, st2;
	if (choice.conflict == NONE &&
		(stat(path1.c_str(), &st1) || stat(path2.c_str(), &st2) ||
		 ((st1.st_dev != st2.st_dev || st1.st_ino != st2.st_ino) && !same_contents(path1, path2, st1.st_size, st2.st_size))))
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

	// Fall back to per-chromosome probing (per-chrom files), regular files before symlinks (master listed no symlinks)
	vector<const ChromFiles::File *> probed;
	for (bool link : {false, true}) {
		for (const ChromFiles::File &file : files.files()) {
			if (file.link == link)
				probed.push_back(&file);
		}
	}
	for (const ChromFiles::File *pfile : probed) {
		const ChromFiles::File &file = *pfile;
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
