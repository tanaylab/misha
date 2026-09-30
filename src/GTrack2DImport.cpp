#include <cstdint>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <unordered_map>

#include "ConfigurationDefaults.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "strutil.h"

#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "GTrack2DImport.h"

using namespace std;
using namespace rdb;

typedef pair<uint32_t, uint32_t> Chrom_pair;

class BinIntervalsFiles : public unordered_map<Chrom_pair, BufferedFile *> {
public:
	BinIntervalsFiles() : unordered_map<Chrom_pair, BufferedFile *>() {}

	~BinIntervalsFiles() {
		for (iterator ifile = begin(); ifile != end(); ++ifile)
			delete ifile->second;
	}
};

static const int64_t INTERVAL_RECORD_SIZE = 4 * sizeof(int64_t) + sizeof(float);

static bool read_interval(BufferedFile &f, int64_t &start1, int64_t &end1, int64_t &start2, int64_t &end2, float &val)
{
	f.read(&start1, sizeof(start1));
	f.read(&end1, sizeof(end1));
	f.read(&start2, sizeof(start2));
	f.read(&end2, sizeof(end2));
	f.read(&val, sizeof(val));

	if (f.eof())
		return false;

	if (f.error())
		verror("Reading file %s: %s\n", f.file_name().c_str(), strerror(errno));

	return true;
}

static void write_interval(BufferedFile &f, int64_t start1, int64_t end1, int64_t start2, int64_t end2, float val)
{
	f.write(&start1, sizeof(start1));
	f.write(&end1, sizeof(end1));
	f.write(&start2, sizeof(start2));
	f.write(&end2, sizeof(end2));
	f.write(&val, sizeof(val));
	if (f.error())
		verror("Writing file %s: %s\n", f.file_name().c_str(), strerror(errno));
}

// Calls f on every interval of the intermediate files, in order, and deletes each file once it is read
template <typename F>
static void for_each_interval(const vector<string> &files, F f)
{
	int64_t start1, end1, start2, end2;
	float val;

	for (const string &fname : files) {
		BufferedFile file;

		if (file.open(fname.c_str(), "r"))
			verror("Opening an intermediate file %s: %s\n", fname.c_str(), strerror(errno));

		while (read_interval(file, start1, end1, start2, end2, val)) {
			f(start1, end1, start2, end2, val);
			check_interrupt();
		}
		file.close();
		unlink(fname.c_str());
	}
}

//------------------------------- SHARED WITH GTrack2DImportContacts.cpp ---------------------------

int import_num_kids(const IntervUtils &iu, uint64_t num_tasks)
{
	if (!iu.get_multitasking())
		return 1;

	uint64_t num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
	uint64_t max_kids = min(min(iu.get_max_processes2core() * num_cores, iu.get_max_processes()), (uint64_t)MAX_KIDS);

	return (int)max(min(max_kids, num_tasks), (uint64_t)1);
}

vector<int> split_files(const vector<int64_t> &sizes, int num_kids)
{
	int num_files = sizes.size();
	uint64_t total = 0;
	uint64_t cum = 0;
	vector<int> bounds(1, 0);

	// +1 per file: an empty file still costs an open
	for (int64_t size : sizes)
		total += size + 1;

	// A block ends once it reaches its share of the total, or when the files left are just enough
	// to give each remaining kid one.
	for (int i = 0; i < num_files - 1; ++i) {
		int num_blocks = bounds.size();

		cum += sizes[i] + 1;
		if (num_blocks < num_kids && (cum * num_kids >= total * num_blocks || num_files - 1 - i == num_kids - num_blocks))
			bounds.push_back(i + 1);
	}
	bounds.push_back(num_files);
	return bounds;
}

string pair_file_name(const string &dirname, int chromid1, int chromid2, int kid)
{
	char buf[100];

	snprintf(buf, sizeof(buf), "/.%d-%d.%d", chromid1, chromid2, kid);
	return dirname + buf;
}

vector<PairFiles> list_pair_files(const string &dirname)
{
	map<Chrom_pair, map<int, string>> pair_kid_files;
	DIR *dir = opendir(dirname.c_str());
	struct dirent *dirp;

	if (!dir)
		verror("Failed to read directory %s: %s\n", dirname.c_str(), strerror(errno));

	while ((dirp = readdir(dir))) {
		int chromid1, chromid2, kid, len = 0;

		if (sscanf(dirp->d_name, ".%d-%d.%d%n", &chromid1, &chromid2, &kid, &len) == 3 && !dirp->d_name[len])
			pair_kid_files[Chrom_pair(chromid1, chromid2)][kid] = dirname + "/" + dirp->d_name;
	}
	closedir(dir);

	vector<PairFiles> pairs;

	for (const auto &ipair : pair_kid_files) {
		PairFiles pair_files;

		pair_files.chromid1 = ipair.first.first;
		pair_files.chromid2 = ipair.first.second;
		pair_files.size = 0;
		for (const auto &ifile : ipair.second) {
			pair_files.size += BufferedFile::file_size(ifile.second.c_str());
			pair_files.files.push_back(ifile.second);
		}
		pairs.push_back(pair_files);
	}
	return pairs;
}

vector<vector<int>> assign_pairs(const vector<PairFiles> &pairs, int num_kids)
{
	vector<int> order(pairs.size());
	vector<uint64_t> loads(num_kids, 0);
	vector<vector<int>> kid_pairs(num_kids);

	iota(order.begin(), order.end(), 0);
	stable_sort(order.begin(), order.end(), [&](int a, int b) { return pairs[a].size > pairs[b].size; });

	for (int ipair : order) {
		int kid = min_element(loads.begin(), loads.end()) - loads.begin();

		kid_pairs[kid].push_back(ipair);
		loads[kid] += pairs[ipair].size + 1;
	}
	return kid_pairs;
}

vector<char> run_kids(IntervUtils &iu, int num_kids, const function<char(int)> &work)
{
	vector<char> res(num_kids);

	if (num_kids == 1) {
		res[0] = work(0);
		return res;
	}

	prepare4multitasking(sizeof(char), 0, num_kids * sizeof(char), iu.get_max_mem_usage(), num_kids);

	for (int kid = 0; kid < num_kids; ++kid) {
		if (!launch_process()) {
			// an error thrown by work() reaches the entry point's catch, whose rerror() ends the kid
			char kid_res = work(kid);

			*(char *)allocate_res(0) = kid_res;
			rexit();
		}
	}

	wait_for_kids(iu);

	for (int kid = 0; kid < num_kids; ++kid)
		res[kid] = *(char *)get_kid_res(kid);
	return res;
}

//--------------------------------------------------------------------------------------------------

// STAGE 1 of one kid: reads files [fbegin, fend) into the kid's intermediate files, one per chromosome pair.
// Returns true if all the intervals are points.
static bool read_input_files(IntervUtils &iu, SEXP _files, const vector<int64_t> &sizes, int fbegin, int fend, const string &dirname, int kid)
{
	Progress_reporter progress;
	BinIntervalsFiles bin_intervals_files;
	bool are_all_points = true;
	int64_t start1, start2, end1, end2;
	float val;

	progress.init(accumulate(sizes.begin() + fbegin, sizes.begin() + fend, (int64_t)0), 10000000);
	for (int ifile = fbegin; ifile < fend; ++ifile) {
		vector<string> fields;
		long lineno = 0;
		BufferedFile infile;
		infile.open(CHAR(STRING_ELT(_files, ifile)), "r");

		lineno += split_line(infile, fields, '\t');

		if (fields.empty())
			continue;

		if (fields.size() < GInterval2D::NUM_COLS + 1)
			verror("File %s, line %ld: invalid format", infile.file_name().c_str(), lineno);

		for (int i = 0; i < GInterval2D::NUM_COLS; ++i) {
			if (fields[i] != GInterval2D::COL_NAMES[i])
				verror("File %s, line %ld: invalid format", infile.file_name().c_str(), lineno);
		}

		while (1) {
			int64_t fpos = infile.tell();
			int chromid1, chromid2;
			char *endptr;

			// read the interval from tab-delimited file
			lineno += split_line(infile, fields, '\t');

			progress.report(infile.tell() - fpos);
			check_interrupt();

			if (fields.empty())
				break;

			if (fields.size() < GInterval2D::NUM_COLS + 1)
				verror("File %s, line %ld: invalid format", infile.file_name().c_str(), lineno);

			try {
				chromid1 = iu.chrom2id(fields[GInterval2D::CHROM1]);
				chromid2 = iu.chrom2id(fields[GInterval2D::CHROM2]);
			} catch (TGLException &) {
				// there might be unrecognized chromosomes, ignore them
				continue;
			}

			start1 = strtoll(fields[GInterval2D::START1].c_str(), &endptr, 10);
			if (*endptr || start1 < 0)
				verror("File %s, line %ld: invalid format of start1 coordinate", infile.file_name().c_str(), lineno);

			end1 = strtoll(fields[GInterval2D::END1].c_str(), &endptr, 10);
			if (*endptr || end1 < 0)
				verror("File %s, line %ld: invalid format of end1 coordinate", infile.file_name().c_str(), lineno);

			if (start1 >= end1)
				verror("File %s, line %ld: start1 coordinate exceeds or equals the end1 coordinate", infile.file_name().c_str(), lineno);

			if ((uint64_t)end1 > iu.get_chromkey().get_chrom_size(chromid1))
				verror("File %s, line %ld: end1 coordinate exceeds chromosome's size", infile.file_name().c_str(), lineno);

			start2 = strtoll(fields[GInterval2D::START2].c_str(), &endptr, 10);
			if (*endptr || start1 < 0)
				verror("File %s, line %ld: invalid format of start2 coordinate", infile.file_name().c_str(), lineno);

			end2 = strtoll(fields[GInterval2D::END2].c_str(), &endptr, 10);
			if (*endptr || end2 < 0)
				verror("File %s, line %ld: invalid format of end2 coordinate", infile.file_name().c_str(), lineno);

			if (start2 >= end2)
				verror("File %s, line %ld: start2 coordinate exceeds or equals the end1 coordinate", infile.file_name().c_str(), lineno);

			if ((uint64_t)end2 > iu.get_chromkey().get_chrom_size(chromid2))
				verror("File %s, line %ld: end2 coordinate exceeds chromosome's size", infile.file_name().c_str(), lineno);

			val = strtod(fields[GInterval2D::NUM_COLS].c_str(), &endptr);
			if (*endptr)
				verror("File %s, line %ld: invalid value", infile.file_name().c_str(), lineno);

			are_all_points &= start1 == end1 - 1 && start2 == end2 - 1;

			// Write down the interval + value in binary format
			BinIntervalsFiles::iterator ibin_intervals_file = bin_intervals_files.find(Chrom_pair(chromid1, chromid2));
			if (ibin_intervals_file == bin_intervals_files.end()) {
				string filename = pair_file_name(dirname, chromid1, chromid2, kid);

				ibin_intervals_file = bin_intervals_files.insert(make_pair(Chrom_pair(chromid1, chromid2), new BufferedFile())).first;
				if (ibin_intervals_file->second->open(filename.c_str(), "wb"))
					verror("Writing an intermediate file %s: %s\n", filename.c_str(), strerror(errno));
			}

			BufferedFile *out_file = ibin_intervals_file->second;
			write_interval(*out_file, start1, end1, start2, end2, val);
		}
	}
	progress.report_last();

	// Close the files here and check: a kid ends without unwinding, and the last buffered writes can still fail.
	// It also keeps the number of open files low for the next stage.
	for (BinIntervalsFiles::iterator ifile = bin_intervals_files.begin(); ifile != bin_intervals_files.end(); ++ifile) {
		if (ifile->second->close())
			verror("Writing an intermediate file %s: %s\n", ifile->second->file_name().c_str(), strerror(errno));
	}

	return are_all_points;
}

// STAGES 2 - 4 for one chromosome pair
static void write_pair(IntervUtils &iu, const string &dirname, const PairFiles &pair_files, bool are_all_points)
{
	// STAGE 2. Check what is the number of contacts in a chromosome pair and based on that number decide how many subtrees would be needed in StatQuadTreeCachedSerializer.
	int chromid1 = pair_files.chromid1;
	int chromid2 = pair_files.chromid2;
	int64_t num_intervals = pair_files.size / INTERVAL_RECORD_SIZE;
	int64_t num_subtrees = max(num_intervals / (int64_t)iu.get_max_data_size(), (int64_t)1);
	num_subtrees = 1 << (2 * (int)(log2(num_subtrees) / 2));  // round the number of subtrees to the lowest power of 4

	GenomeTrackRectsPoints gtrack_points(iu.get_track_chunk_size(), iu.get_track_num_chunks());
	GenomeTrackRectsRects gtrack_rects(iu.get_track_chunk_size(), iu.get_track_num_chunks());
	PointsQuadTreeCachedSerializer points_serializer;
	RectsQuadTreeCachedSerializer rects_serializer;
	string filename = dirname + "/" + GenomeTrack::get_2d_filename(iu.get_chromkey(), chromid1, chromid2);

	if (are_all_points) {
		gtrack_points.init_write(filename.c_str(), chromid1, chromid2);
		gtrack_points.init_serializer(points_serializer, 0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2), num_subtrees, true);
	} else {
		gtrack_rects.init_write(filename.c_str(), chromid1, chromid2);
		gtrack_rects.init_serializer(rects_serializer, 0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2), num_subtrees, true);
	}

	// the files to read for each subtree
	vector<vector<string>> subtree_srcs(1, pair_files.files);

	// STAGE 3: Split the intervals of a chromosome pair to binary files - each holding contacts of a subtree.
	if (num_subtrees > 1) {
		const Rectangles &subarenas = are_all_points ? points_serializer.get_subarenas() : rects_serializer.get_subarenas();
		BufferedFiles subtrees_files(num_subtrees);

		subtree_srcs.resize(num_subtrees);
		for (uint64_t i = 0; i < subtrees_files.size(); ++i) {
			char buf[FILENAME_MAX];

			snprintf(buf, sizeof(buf), "%s/.%d-%d.s%ld", dirname.c_str(), chromid1, chromid2, (long)i);
			subtree_srcs[i].assign(1, buf);
			subtrees_files[i] = new BufferedFile();
			if (subtrees_files[i]->open(buf, "w"))
				verror("Opening an intermediate file %s: %s\n", buf, strerror(errno));
		}

		for_each_interval(pair_files.files, [&](int64_t start1, int64_t end1, int64_t start2, int64_t end2, float val) {
			Rectangle rect(start1, start2, end1, end2);
			for (Rectangles::const_iterator isubarena = subarenas.begin(); isubarena != subarenas.end(); ++isubarena) {
				if (rect.do_intersect(*isubarena)) {
					write_interval(*subtrees_files[isubarena - subarenas.begin()], start1, end1, start2, end2, val);
					break;
				}
			}
		});

		for (uint64_t i = 0; i < subtrees_files.size(); ++i) {
			if (subtrees_files[i]->close())
				verror("Writing an intermediate file %s: %s\n", subtrees_files[i]->file_name().c_str(), strerror(errno));
		}
	}

	// Stage 4: Read the contacts of a subtree and insert them to StatQuadTreeCachedSerializer.
	for (const vector<string> &srcs : subtree_srcs) {
		for_each_interval(srcs, [&](int64_t start1, int64_t end1, int64_t start2, int64_t end2, float val) {
			if (are_all_points)
				points_serializer.insert(PointsQuadTree::ValueType(start1, start2, val));
			else
				rects_serializer.insert(RectsQuadTree::ValueType(start1, start2, end1, end2, val));
		});
	}

	if (are_all_points)
		points_serializer.end();
	else
		rects_serializer.end();
}

extern "C" {

SEXP gtrack_2d_import(SEXP _track, SEXP _files, SEXP _envir)
{
	try {
		string dirname;
		vector<PairFiles> pairs;
		bool are_all_points = true;

		// Each stage has its own RdbInitializer: the multitasking state (shared memory, kid bookkeeping)
		// is set up once per RdbInitializer, so the kids of the second stage need a fresh one.
		{
			RdbInitializer rdb_init;

			if (!Rf_isString(_track) || Rf_length(_track) != 1)
				verror("Track argument is not a string");

			if (!Rf_isString(_files) || Rf_length(_files) < 1)
				verror("Files argument is not a vector of strings");

			IntervUtils iu(_envir);
			const char *track = CHAR(STRING_ELT(_track, 0));
			dirname = create_track_dir(_envir, track);

			// The number of 2D intervals might be huge. We might not be even able to hold all the intervals of one chromosome pair in memory and to build the quad tree with it.
			// Our strategy is therefore that:
			// 1. Read the input files and split them into binary files each holding the intervals of specific chromosomes pair.
			// 2. Check what is the number of contacts in a chromosome pair and based on that number decide how many subtrees would be needed in StatQuadTreeCachedSerializer.
			//    (Creating a quad tree through seriaizer would be the only feasible method when only part of the intervals of a chromosome pair can be loaded into RAM.)
			// 3. Split the intervals of a chromosome pair to binary files - each holding contacts of a subtree. We assume that we will be able to hold in the memory the contacts of one subtree.
			// 4. Read the contacts of a subtree and insert them to StatQuadTreeCachedSerializer.
			// Step 1 runs in parallel over the input files, steps 2 - 4 over the chromosome pairs (see GTrack2DImport.h).

			REprintf("Reading input file(s)...\n");

			vector<int64_t> sizes;
			for (int ifile = 0; ifile < Rf_length(_files); ++ifile) {
				const char *fname = CHAR(STRING_ELT(_files, ifile));
				struct stat st;

				if (stat(fname, &st))
					verror("Accessing file %s: %s", fname, strerror(errno));
				sizes.push_back(st.st_size);
			}

			// STEP 1: Read the input files and split them into binary files each holding the intervals of specific chromosomes pair.
			vector<int> bounds = split_files(sizes, import_num_kids(iu, sizes.size()));
			vector<char> kids_all_points = run_kids(iu, bounds.size() - 1, [&](int kid) {
					return (char)read_input_files(iu, _files, sizes, bounds[kid], bounds[kid + 1], dirname, kid);
				});

			for (char kid_all_points : kids_all_points)
				are_all_points &= (bool)kid_all_points;

			pairs = list_pair_files(dirname);
		}

		{
			RdbInitializer rdb_init;
			IntervUtils iu(_envir);
			uint64_t num_records = 0;

			for (const PairFiles &pair_files : pairs)
				num_records += pair_files.size / INTERVAL_RECORD_SIZE;

			REprintf("Writing the track...\n");

			int num_kids = import_num_kids(iu, min((uint64_t)pairs.size(), num_records / misha::config::MIN_RECORDS_PER_PROCESS));
			vector<vector<int>> kid_pairs = assign_pairs(pairs, num_kids);

			run_kids(iu, num_kids, [&](int kid) {
					Progress_reporter progress;

					progress.init(kid_pairs[kid].size(), 1);
					for (int ipair : kid_pairs[kid]) {
						write_pair(iu, dirname, pairs[ipair], are_all_points);
						progress.report(1);
					}
					progress.report_last();
					return (char)0;
				});
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	}

	return R_NilValue;
}

}
