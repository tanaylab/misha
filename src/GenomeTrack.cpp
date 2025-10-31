/*
 * GenomeTrack.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: hoichman
 */

#include <cstdint>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_set>

#include "GenomeTrack.h"
#include "TrackIndex.h"
#include "TGLException.h"
#include "rdbutils.h"

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

const pair<int, int> GenomeTrack::get_chromid_2d(const GenomeChromKey &chromkey, const string &filename)
{
	uint64_t pos = filename.find_first_of("-");

	if (pos == string::npos)
		TGLError<GenomeTrack>(NOT_2D, "File %s does not belong to 2D track", filename.c_str());

	string chrom1(filename, 0, pos);
	string chrom2(filename, pos + 1);

	return pair<int, int>(chromkey.chrom2id(chrom1), chromkey.chrom2id(chrom2));
}

string GenomeTrack::find_existing_1d_filename(const GenomeChromKey &chromkey, const string &track_dir, int chromid)
{
	const string &base = get_1d_filename(chromkey, chromid);
	vector<string> candidates;
	candidates.push_back(base);

	if (base.compare(0, 3, "chr") == 0 && base.size() > 3)
		candidates.push_back(base.substr(3));
	else
		candidates.push_back(string("chr") + base);

	vector<string> aliases;
	chromkey.get_aliases(chromid, aliases);
	candidates.insert(candidates.end(), aliases.begin(), aliases.end());

	unordered_set<string> seen;
	for (const string &candidate : candidates) {
		if (candidate.empty())
			continue;
		if (!seen.insert(candidate).second)
			continue;
		string full = track_dir + "/" + candidate;
		if (access(full.c_str(), F_OK) == 0)
			return candidate;
	}

	return base;
}


GenomeTrack::Type GenomeTrack::get_type(const char *track_dir, const GenomeChromKey &chromkey, bool return_obsolete_types)
{
	if (access(track_dir, F_OK))
		TGLError<GenomeTrack>(FILE_ERROR, "Accessing directory %s: %s\n", track_dir, strerror(errno));

	// First, try to read from track.idx if it exists
	std::string idx_path = std::string(track_dir) + "/track.idx";
	struct stat st;
	if (stat(idx_path.c_str(), &st) == 0) {
		try {
			auto idx = get_track_index(track_dir);
			if (idx) {
				// Map MishaTrackType to GenomeTrack::Type
				switch (idx->get_track_type()) {
					case MishaTrackType::DENSE:  return FIXED_BIN;
					case MishaTrackType::SPARSE: return SPARSE;
					case MishaTrackType::ARRAY:  return ARRAYS;
					default:
						// Unknown type, fall through to per-chromosome probing
						break;
				}
			}
		} catch (...) {
			// Fall through to per-chromosome probing on any index error
		}
	}

	// Fall back to per-chromosome probing (per-chrom files)
	vector<string> filenames;
	rdb::get_chrom_files(track_dir, filenames);

	sort(filenames.begin(), filenames.end());

	for (const string &fname : filenames) {
		string fullpath = string(track_dir) + "/" + fname;
		bool is_2d = fname.find('-') != string::npos;

		try {
			Type type = s_read_type(fullpath.c_str());

			if (!is_2d) {
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
	if (filenames.empty()) {
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
	umask(07);

	if (m_bfile.open(filename, mode))
		TGLError<GenomeTrack>(FILE_ERROR, "Opening a track file %s: %s", filename, strerror(errno));

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
