#include <cstdint>
#include <algorithm>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <numeric>
#include <unordered_map>

#include "ConfigurationDefaults.h"
#include "HashFunc.h"
#include "rdbinterval.h"
#include "rdbprogress.h"
#include "rdbutils.h"
#include "strutil.h"

#include "GenomeTrackRects.h"
#include "GenomeTrackSparse.h"
#include "GTrack2DImport.h"

using namespace std;
using namespace rdb;

// The checks of the 1D importers (GenomeTrackImportWig.cpp)
void verify_chroms_matched(const GenomeChromKey &chromkey, const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);
void report_skipped_chroms(const char *fname, const UnknownChroms &unknown, uint64_t num_matched, SEXP envir);

typedef pair<uint32_t, uint32_t> Chrom_pair;

class Contact_chrom_files : public unordered_map<Chrom_pair, BufferedFile *> {
public:
	Contact_chrom_files() : unordered_map<Chrom_pair, BufferedFile *>() {}

	~Contact_chrom_files() {
		for (iterator ifile = begin(); ifile != end(); ++ifile)
			delete ifile->second;
	}

	// Closes the files and checks: a kid ends without unwinding, and the last buffered writes can still fail.
	// It also keeps the number of open files low for the next stage.
	void close() {
		for (iterator ifile = begin(); ifile != end(); ++ifile) {
			if (ifile->second->close())
				verror("Writing an intermediate file %s: %s\n", ifile->second->file_name().c_str(), strerror(errno));
		}
	}
};

// Fend ids index arrays (see process_contacts_as_fends), so they are bounded
static const int64_t MAX_FEND_ID = INT32_MAX;

static const int64_t CONTACT_RECORD_SIZE = 2 * sizeof(int64_t) + sizeof(float);

// Peak memory of building a chromosome pair, in bytes per contact (see pair_mem): the sorted vector of the
// records that sums the duplicates plus the quad trees, which hold every contact twice (mirrored in a cis pair,
// one tree per orientation in a trans pair). Measured as the peak RSS of a serial import of one pair of 2M - 40M
// Hi-C-like contacts: at most 167.3 bytes a contact (cis, no duplicates; trans and duplicates need less).
// Rounded up to a multiple of 8.
static const uint64_t CONTACT_MEM = 168;

static unsigned read_header(BufferedFile &file, const char *fname, const char *ftype, const char *colnames[], int num_cols, vector<int> &fcol_idx)
{
	vector<string> fields;

	if (file.open(fname, "r"))
		verror("Failed to open %s file %s: %s", ftype, fname, strerror(errno));

	fcol_idx.resize(num_cols, -1);

	split_line(file, fields, '\t');
	if (fields.empty())
		verror("Invalid format of %s file %s", ftype, fname);

	if (file.error())
		verror("Error while reading %s file %s: %s", ftype, fname, strerror(errno));

	for (int i = 0; i < num_cols; ++i) {
		for (vector<string>::const_iterator ifield = fields.begin(); ifield != fields.end(); ++ifield) {
			if (*ifield == colnames[i]) {
				if (fcol_idx[i] == -1)
					fcol_idx[i] = ifield - fields.begin();
				else
					verror("Invalid format of %s file %s: column %s appears more than once", ftype, fname, colnames[i]);
			}
		}

		if (fcol_idx[i] < 0)
			verror("Invalid format of %s file %s: column %s was not found", ftype, fname, colnames[i]);
	}

	return fields.size();
}

static void write_contact(BufferedFile &file, int64_t coord1, int64_t coord2, float value)
{
	file.write(&coord1, sizeof(coord1));
	file.write(&coord2, sizeof(coord2));
	file.write(&value, sizeof(value));
	if (file.error())
		verror("Writing file %s: %s\n", file.file_name().c_str(), strerror(errno));
}

// Calls f on every contact of the intermediate files, in order, and deletes each file once it is read
template <typename F>
static void for_each_contact(const vector<string> &files, F f)
{
	int64_t coord1, coord2;
	float value;

	for (const string &fname : files) {
		BufferedFile file;

		if (file.open(fname.c_str(), "r"))
			verror("Opening an intermediate file %s: %s\n", fname.c_str(), strerror(errno));

		while (1) {
			file.read(&coord1, sizeof(coord1));
			file.read(&coord2, sizeof(coord2));
			file.read(&value, sizeof(value));

			if (file.eof())
				break;

			if (file.error())
				verror("Reading file %s: %s\n", file.file_name().c_str(), strerror(errno));

			f(coord1, coord2, value);
			check_interrupt();
		}
		file.close();
		unlink(fname.c_str());
	}
}

static vector<int64_t> get_files_sizes(SEXP _files)
{
	vector<int64_t> sizes;

	for (int ifile = 0; ifile < Rf_length(_files); ++ifile) {
		const char *fname = CHAR(STRING_ELT(_files, ifile));
		struct stat st;

		if (stat(fname, &st))
			verror("Accessing file %s: %s", fname, strerror(errno));
		sizes.push_back(st.st_size);
	}
	return sizes;
}

// Returns the number of chromosome names of the input that the database has
static uint64_t process_contacts_as_intervals(IntervUtils &iu, SEXP _files, const string &dirname, UnknownChroms &unknown)
{
	REprintf("Reading input file(s)...\n");

	vector<int64_t> sizes = get_files_sizes(_files);

	// STEP 1: Read the input files and split them into binary files each holding the intervals of specific chromosomes pair.
	// Each kid reads a contiguous block of the files.
	vector<int> bounds = split_files(sizes, import_num_kids(iu, sizes.size()));

	run_kids(iu, bounds.size() - 1, [&](int kid) {
		Progress_reporter progress;
		Contact_chrom_files contact_chrom_files;
		KidChroms kid_chroms;
		int64_t start1, start2, end1, end2;
		float val;

		progress.init(accumulate(sizes.begin() + bounds[kid], sizes.begin() + bounds[kid + 1], (int64_t)0), 10000000);
		for (int ifile = bounds[kid]; ifile < bounds[kid + 1]; ++ifile) {
			vector<string> fields;
			long lineno = 0;
			BufferedFile infile;
			if (infile.open(CHAR(STRING_ELT(_files, ifile)), "r"))
				verror("Failed to open file %s: %s", CHAR(STRING_ELT(_files, ifile)), strerror(errno));

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

				// there might be unrecognized chromosomes: skip them, and report them once the track is written
				chromid1 = iu.get_chromkey().find_chromid(fields[GInterval2D::CHROM1]);
				chromid2 = iu.get_chromkey().find_chromid(fields[GInterval2D::CHROM2]);
				kid_chroms.note(chromid1, fields[GInterval2D::CHROM1]);
				kid_chroms.note(chromid2, fields[GInterval2D::CHROM2]);
				if (chromid1 < 0 || chromid2 < 0)
					continue;

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
				if (*endptr || start2 < 0)
					verror("File %s, line %ld: invalid format of start2 coordinate", infile.file_name().c_str(), lineno);

				end2 = strtoll(fields[GInterval2D::END2].c_str(), &endptr, 10);
				if (*endptr || end2 < 0)
					verror("File %s, line %ld: invalid format of end2 coordinate", infile.file_name().c_str(), lineno);

				if (start2 >= end2)
					verror("File %s, line %ld: start2 coordinate exceeds or equals the end2 coordinate", infile.file_name().c_str(), lineno);

				if ((uint64_t)end2 > iu.get_chromkey().get_chrom_size(chromid2))
					verror("File %s, line %ld: end2 coordinate exceeds chromosome's size", infile.file_name().c_str(), lineno);

				val = strtod(fields[GInterval2D::NUM_COLS].c_str(), &endptr);
				if (*endptr)
					verror("File %s, line %ld: invalid value", infile.file_name().c_str(), lineno);

				int64_t coord1 = (start1 + end1) / 2;
				int64_t coord2 = (start2 + end2) / 2;

				if (chromid1 > chromid2 || (chromid1 == chromid2 && coord1 > coord2)) {
					swap(chromid1, chromid2);
					swap(coord1, coord2);
				}


				// Write down the interval + value in binary format
				Chrom_pair chrom_pair(chromid1, chromid2);
				Contact_chrom_files::iterator icontact_chrom_file = contact_chrom_files.find(chrom_pair);

				if (icontact_chrom_file == contact_chrom_files.end()) {
					string filename = pair_file_name(dirname, chromid1, chromid2, kid);

					icontact_chrom_file = contact_chrom_files.insert(make_pair(chrom_pair, new BufferedFile())).first;
					if (icontact_chrom_file->second->open(filename.c_str(), "wb"))
						verror("Writing an intermediate file %s: %s\n", filename.c_str(), strerror(errno));
				}

				write_contact(*icontact_chrom_file->second, coord1, coord2, val);
			}
		}
		progress.report_last();
		contact_chrom_files.close();
		kid_chroms.save(dirname, kid);
		return (char)0;
	});

	set<int> matched;
	load_kid_chroms(dirname, bounds.size() - 1, matched, unknown);
	return matched.size();
}

// unknown receives the chromosome names of the fends file that the database does not have; the contacts
// of their fends are skipped. Returns the number of chromosomes of the fends file that the database has.
static uint64_t process_contacts_as_fends(IntervUtils &iu, SEXP _contacts, SEXP _fends, const string &dirname, UnknownChroms &unknown)
{
	enum { COL_FEND, COL_CHR, COL_COORD, NUM_FENDS_COLS };
	const char *FENDS_COLS[] = { "fend", "chr", "coord" };
	vector<int> fends_col_idx;
	unsigned fends_cols_total;
	const char *fends_fname = CHAR(STRING_ELT(_fends, 0));

	enum { COL_FEND1, COL_FEND2, COL_COUNT, NUM_CONTACTS_COLS };
	const char *CONTACTS_COLS[] = { "fend1", "fend2", "count" };

	vector<string> fields;
	vector<int64_t> sizes = get_files_sizes(_contacts);

	Progress_reporter progress;
	BufferedFile fends_file;
	fends_cols_total = read_header(fends_file, fends_fname, "fragment ends", FENDS_COLS, NUM_FENDS_COLS, fends_col_idx);

	// read fends file
	REprintf("Reading fends file...\n");

	// Indexed by fend id, so their size follows the largest id (12 bytes an id): real fends files number their
	// fends 1..N. chromids: FEND_NOT_LISTED, FEND_UNKNOWN_CHROM, or the chromosome of the fend.
	enum { FEND_NOT_LISTED = -2, FEND_UNKNOWN_CHROM = -1 };
	vector<int> chromids;
	vector<int64_t> coords;
	set<int> matched_chromids;
	int64_t lineno;
	char *endptr;
	int64_t fendid;
	int64_t coord;

	progress.init(fends_file.file_size(), 10000000);
	lineno = 1;
	while (1) {
		int64_t fpos = fends_file.tell();

		lineno++;
		split_line(fends_file, fields, '\t');
		if (fields.empty())
			break;

		if (fields.size() != fends_cols_total)
			verror("Error reading fragment ends file %s, line %ld: invalid file format", fends_fname, lineno);

		const char *fendid_str = fields[fends_col_idx[COL_FEND]].c_str();
		const char *chr_str = fields[fends_col_idx[COL_CHR]].c_str();
		const char *coord_str = fields[fends_col_idx[COL_COORD]].c_str();

		errno = 0;
		fendid = strtoll(fendid_str, &endptr, 10);
		if (*endptr || endptr == fendid_str)
			verror("Error reading fragment ends file %s, line %ld: invalid %s", fends_fname, lineno, FENDS_COLS[COL_FEND]);

		if (fendid < 0 || errno == ERANGE || fendid > MAX_FEND_ID)
			verror("Error reading fragment ends file %s, line %ld: %s %s is out of range (0 - %lld)", fends_fname, lineno,
				   FENDS_COLS[COL_FEND], fendid_str, (long long)MAX_FEND_ID);

		// grow only: the ids need not be sorted, and shrinking would drop the fends read so far
		if ((uint64_t)fendid >= chromids.size()) {
			chromids.resize(fendid + 1, FEND_NOT_LISTED);
			coords.resize(fendid + 1);
		}

		if (chromids[fendid] != FEND_NOT_LISTED)
			verror("Error reading fragment ends file %s, line %ld: %s %lld appears more than once", fends_fname, lineno,
				   FENDS_COLS[COL_FEND], (long long)fendid);

		// Use chromosome name as-is (no forced chr prefix)
		chromids[fendid] = iu.get_chromkey().find_chromid(chr_str);
		if (chromids[fendid] < 0)
			chromids[fendid] = FEND_UNKNOWN_CHROM;

		coord = strtoll(coord_str, &endptr, 10);
		if (*endptr)
			verror("Error reading fragment ends file %s, line %ld: invalid %s", fends_fname, lineno, FENDS_COLS[COL_COORD]);

		if (chromids[fendid] == FEND_UNKNOWN_CHROM)
			unknown.record(chr_str);
		else if (coord < 0 || coord >= (int64_t)iu.get_chromkey().get_chrom_size(chromids[fendid]))
			verror("Error reading fragment ends file %s, line %ld: %s is out of range", fends_fname, lineno, FENDS_COLS[COL_COORD]);
		else
			matched_chromids.insert(chromids[fendid]);

		coords[fendid] = coord;
		progress.report(fends_file.tell() - fpos);
		check_interrupt();
	}
	progress.report_last();
	fends_file.close();

	// STEP 1: Read the contacts and split them into binary files each holding the contacts of specific pair of chromosomes.
	// Each kid reads a contiguous block of the contacts files.
	REprintf("Reading contacts...\n");

	vector<int> bounds = split_files(sizes, import_num_kids(iu, sizes.size()));

	run_kids(iu, bounds.size() - 1, [&](int kid) {
		Progress_reporter progress;
		Contact_chrom_files contact_chrom_files;
		vector<string> fields;
		int64_t lineno;
		char *endptr;
		float value;

		progress.init(accumulate(sizes.begin() + bounds[kid], sizes.begin() + bounds[kid + 1], (int64_t)0), 10000000);

		for (int icontact_file = bounds[kid]; icontact_file < bounds[kid + 1]; ++icontact_file) {
			int64_t fend1, fend2;
			const char *contacts_fname = CHAR(STRING_ELT(_contacts, icontact_file));
			vector<int> contacts_col_idx;
			unsigned contacts_cols_total;
			BufferedFile contacts_file;

			contacts_cols_total = read_header(contacts_file, contacts_fname, "contacts", CONTACTS_COLS, NUM_CONTACTS_COLS, contacts_col_idx);

			lineno = 1;
			while (1) {
				int64_t fpos = contacts_file.tell();

				lineno++;
				split_line(contacts_file, fields, '\t');
				if (fields.empty())
					break;

				if (fields.size() != contacts_cols_total)
					verror("Error reading contacts ends file %s, line %ld: invalid file format", contacts_fname, lineno);

				const char *fend1_str = fields[contacts_col_idx[COL_FEND1]].c_str();
				const char *fend2_str = fields[contacts_col_idx[COL_FEND2]].c_str();
				const char *value_str = fields[contacts_col_idx[COL_COUNT]].c_str();

				fend1 = strtoll(fend1_str, &endptr, 10);
				if (*endptr)
					verror("Error reading contacts file %s, line %ld: invalid %s", contacts_fname, lineno, CONTACTS_COLS[COL_FEND1]);
				if (fend1 < 0 || fend1 >= (int64_t)chromids.size() || chromids[fend1] < 0)
					continue;

				fend2 = strtoll(fend2_str, &endptr, 10);
				if (*endptr)
					verror("Error reading contacts file %s, line %ld: invalid %s", contacts_fname, lineno, CONTACTS_COLS[COL_FEND2]);
				if (fend2 < 0 || fend2 >= (int64_t)chromids.size() || chromids[fend2] < 0)
					continue;

				value = strtod(value_str, &endptr);
				if (*endptr)
					verror("Error reading contacts file %s, line %ld: invalid %s", contacts_fname, lineno, CONTACTS_COLS[COL_COUNT]);

				if (chromids[fend1] > chromids[fend2] || (chromids[fend1] == chromids[fend2] && coords[fend1] > coords[fend2]))
					swap(fend1, fend2);

				Chrom_pair chrom_pair(chromids[fend1], chromids[fend2]);
				Contact_chrom_files::iterator icontact_chrom_file = contact_chrom_files.find(chrom_pair);

				if (icontact_chrom_file == contact_chrom_files.end()) {
					string filename = pair_file_name(dirname, chromids[fend1], chromids[fend2], kid);

					icontact_chrom_file = contact_chrom_files.insert(make_pair(chrom_pair, new BufferedFile())).first;
					if (icontact_chrom_file->second->open(filename.c_str(), "wb"))
						verror("Writing an intermediate file %s: %s\n", filename.c_str(), strerror(errno));
				}

				write_contact(*icontact_chrom_file->second, coords[fend1], coords[fend2], value);
				progress.report(contacts_file.tell() - fpos);
				check_interrupt();
			}
		}
		progress.report_last();
		contact_chrom_files.close();
		return (char)0;
	});
	return matched_chromids.size();
}

// The number of subtrees of a pair's quad trees: a cis pair holds every contact twice, mirrored around the diagonal
static int64_t contacts_num_subtrees(const IntervUtils &iu, const PairFiles &pair_files)
{
	int64_t num_contacts = pair_files.size / CONTACT_RECORD_SIZE;

	return pair_num_subtrees(iu, pair_files.chromid1 == pair_files.chromid2 ? 2 * num_contacts : num_contacts);
}

// STAGE 2 for a chromosome pair whose quad trees have several subtrees: splits its contacts into one file per
// subtree, mirrored around the diagonal in a cis pair. Returns the number of contacts of the largest subtree,
// which the memory of the pair is estimated by.
static uint64_t split_pair(IntervUtils &iu, const string &dirname, const PairFiles &pair_files)
{
	int chromid1 = pair_files.chromid1;
	int chromid2 = pair_files.chromid2;
	int64_t num_subtrees = contacts_num_subtrees(iu, pair_files);
	Rectangles subarenas = PointsQuadTreeCachedSerializer::split_arena(0, 0, iu.get_chromkey().get_chrom_size(chromid1),
																		iu.get_chromkey().get_chrom_size(chromid2), num_subtrees);
	BufferedFiles subtrees_files(num_subtrees);
	vector<uint64_t> num_contacts(num_subtrees, 0);

	for (uint64_t i = 0; i < subtrees_files.size(); ++i) {
		string fname = subtree_file_name(dirname, chromid1, chromid2, i);

		subtrees_files[i] = new BufferedFile();
		if (subtrees_files[i]->open(fname.c_str(), "w"))
			verror("Opening an intermediate file %s: %s\n", fname.c_str(), strerror(errno));
	}

	auto write_to_subtree = [&](int64_t coord1, int64_t coord2, float value) {
		Point point(coord1, coord2);
		for (Rectangles::const_iterator isubarena = subarenas.begin(); isubarena != subarenas.end(); ++isubarena) {
			if (point.do_intersect(*isubarena)) {
				write_contact(*subtrees_files[isubarena - subarenas.begin()], coord1, coord2, value);
				num_contacts[isubarena - subarenas.begin()]++;
				break;
			}
		}
	};

	for_each_contact(pair_files.files, [&](int64_t coord1, int64_t coord2, float value) {
		write_to_subtree(coord1, coord2, value);

		// if it's a pair of identical chromosomes, mirror the contact around the diagonal
		if (chromid1 == chromid2 && coord1 != coord2)
			write_to_subtree(coord2, coord1, value);
	});

	for (uint64_t i = 0; i < subtrees_files.size(); ++i) {
		if (subtrees_files[i]->close())
			verror("Writing an intermediate file %s: %s\n", subtrees_files[i]->file_name().c_str(), strerror(errno));
	}

	return *max_element(num_contacts.begin(), num_contacts.end());
}

// STAGE 3 for one chromosome pair: builds its quad trees from its intermediate files, or from its subtree files
// if split_pair split it.
static void write_pair(IntervUtils &iu, const string &dirname, const PairFiles &pair_files, bool allow_duplicates)
{
	uint64_t chromid1 = pair_files.chromid1;
	uint64_t chromid2 = pair_files.chromid2;
	int64_t num_subtrees = contacts_num_subtrees(iu, pair_files);

	GenomeTrackRectsPoints gtrack1(iu.get_track_chunk_size(), iu.get_track_num_chunks());
	GenomeTrackRectsPoints gtrack2(iu.get_track_chunk_size(), iu.get_track_num_chunks());
	PointsQuadTreeCachedSerializer qtree_serializer1, qtree_serializer2;
	string filename = dirname + "/" + GenomeTrack::get_2d_filename(iu.get_chromkey(), chromid1, chromid2);

	gtrack1.init_write(filename.c_str(), chromid1, chromid2);
	gtrack1.init_serializer(qtree_serializer1, 0, 0, iu.get_chromkey().get_chrom_size(chromid1), iu.get_chromkey().get_chrom_size(chromid2), num_subtrees, false);

	if (chromid1 != chromid2) {
		filename = dirname + "/" + GenomeTrack::get_2d_filename(iu.get_chromkey(), chromid2, chromid1);
		gtrack2.init_write(filename.c_str(), chromid2, chromid1);
		gtrack2.init_serializer(qtree_serializer2, 0, 0, iu.get_chromkey().get_chrom_size(chromid2), iu.get_chromkey().get_chrom_size(chromid1), num_subtrees, false);
	}

	// the files to read for each subtree
	vector<vector<string>> subtree_srcs(1, pair_files.files);

	if (num_subtrees > 1) {
		subtree_srcs.clear();
		for (int64_t i = 0; i < num_subtrees; ++i)
			subtree_srcs.push_back(vector<string>(1, subtree_file_name(dirname, chromid1, chromid2, i)));
	}

	struct Contact {
		uint64_t coord1;
		uint64_t coord2;
		float value;
	};

	// Read the contacts of a subtree, sum up the duplicates and insert them to StatQuadTreeCachedSerializer
	// in coordinate order. A hash map used to sum them, and its iteration order, which depends on the C++
	// library, decided the bytes of the track files. Sorting needs less memory than the map did.
	for (const vector<string> &srcs : subtree_srcs) {
		uint64_t num_records = 0;

		for (const string &fname : srcs) {
			struct stat st;

			if (stat(fname.c_str(), &st))
				verror("Accessing an intermediate file %s: %s\n", fname.c_str(), strerror(errno));
			num_records += st.st_size / CONTACT_RECORD_SIZE;
		}

		vector<Contact> contacts;
		contacts.reserve(num_records);
		for_each_contact(srcs, [&](int64_t coord1, int64_t coord2, float value) {
			contacts.push_back(Contact{(uint64_t)coord1, (uint64_t)coord2, value});
		});

		// stable: the duplicates of a contact are summed in the order they were read
		stable_sort(contacts.begin(), contacts.end(), [](const Contact &a, const Contact &b) {
			return a.coord1 < b.coord1 || (a.coord1 == b.coord1 && a.coord2 < b.coord2);
		});

		for (size_t i = 0; i < contacts.size(); ) {
			int64_t coord1 = contacts[i].coord1;
			int64_t coord2 = contacts[i].coord2;
			double sum = contacts[i].value;
			size_t j = i + 1;

			for (; j < contacts.size() && contacts[j].coord1 == contacts[i].coord1 && contacts[j].coord2 == contacts[i].coord2; ++j)
				sum += contacts[j].value;

			if (j - i > 1 && !allow_duplicates)
				verror("Duplicated contact (%s, %ld)-(%s, %ld)", iu.id2chrom(chromid1).c_str(), coord1, iu.id2chrom(chromid2).c_str(), coord2);

			float value = sum;

			qtree_serializer1.insert(PointsQuadTree::ValueType(coord1, coord2, value));

			if (chromid1 != chromid2)
				qtree_serializer2.insert(PointsQuadTree::ValueType(coord2, coord1, value));
			else if (num_subtrees == 1 && coord1 != coord2)
				qtree_serializer1.insert(PointsQuadTree::ValueType(coord2, coord1, value));

			i = j;
			check_interrupt();
		}

		check_interrupt();
	}

	qtree_serializer1.end();

	if (chromid1 != chromid2)
		qtree_serializer2.end();
}

extern "C" {

SEXP gtrack_import_contacts(SEXP _track, SEXP _contacts, SEXP _fends, SEXP _allow_duplicates, SEXP _envir)
{
	try {
		string dirname;
		string srcname;
		vector<PairFiles> pairs;
		UnknownChroms unknown;
		uint64_t num_matched = 0;

		// Each stage has its own RdbInitializer: the multitasking state (shared memory, kid bookkeeping)
		// is set up once per RdbInitializer, so the kids of each stage need a fresh one. A new RdbInitializer
		// also clears a pending Ctrl-C, so stages 1 and 2 check for one before they end.
		{
			RdbInitializer rdb_init;

			if (!Rf_isString(_track) || Rf_length(_track) != 1)
				verror("Track argument is not a string");

			if (!Rf_isNull(_fends) && (!Rf_isString(_fends) || Rf_length(_fends) != 1))
				verror("Fends argument is not a string");

			if (!Rf_isString(_contacts) || Rf_length(_contacts) < 1)
				verror("Contacts argument is not a string");

			if (!Rf_isLogical(_allow_duplicates) || Rf_length(_allow_duplicates) != 1)
				verror("Allow duplicates argument is not a boolean");

			IntervUtils iu(_envir);
			const char *track = CHAR(STRING_ELT(_track, 0));
			dirname = create_track_dir(_envir, track);

			// The number of contacts might be huge. We might not be even able to hold the contacts of one chromosome in the memory and to build the quad tree with it.
			// Our strategy is therefore that:
			// 1. Read the contacts and split them into binary files each holding the contacts of specific pair of chromosomes.
			// 2. Check what is the number of contacts in a chromosome pair and based on that number decide how many subtrees would be needed in StatQuadTreeCachedSerializer.
			//    (Creating a quad tree through seriaizer would be the only feasible method considering the enormous number of contacts.)
			//    Split contacts of a chromosome pair of several subtrees to binary files - each holding contacts of a subtree. We assume that we will be able to hold in the memory the contacts of one subtree.
			// 3. Read the contacts of a subtree, sum up the duplicates and insert them to StatQuadTreeCachedSerializer.
			// Step 1 runs in parallel over the input files, steps 2 and 3 over the chromosome pairs (see GTrack2DImport.h).

			// STAGE 1: Read the contacts and split them into binary files each holding the contacts of specific pair of chromosomes.
			vector<string> fnames;

			for (int ifile = 0; ifile < Rf_length(_contacts); ++ifile)
				fnames.push_back(CHAR(STRING_ELT(_contacts, ifile)));

			if (Rf_isNull(_fends)) {
				srcname = import_src_name(fnames);
				num_matched = process_contacts_as_intervals(iu, _contacts, dirname, unknown);
			} else {
				srcname = CHAR(STRING_ELT(_fends, 0));
				num_matched = process_contacts_as_fends(iu, _contacts, _fends, dirname, unknown);
			}

			pairs = list_pair_files(dirname);
			verify_chroms_matched(iu.get_chromkey(), srcname.c_str(), unknown, num_matched, _envir);
			if (pairs.empty()) {
				if (!Rf_isNull(_fends))
					verror("No contact of %s could be placed: none was read, or their fragment ends are missing from %s or lie on chromosomes "
						   "that the genome database does not have. A 2D track needs at least one contact.",
						   import_src_name(fnames).c_str(), srcname.c_str());
				if (num_matched)
					verror("No contact of %s has both of its chromosomes in the genome database, so nothing is left to import: a 2D track needs at least one contact",
						   srcname.c_str());
				verror("No contacts were read from %s: a 2D track needs at least one", srcname.c_str());
			}
			check_interrupt();
		}

		// the contacts that the memory estimate of each pair counts: all of them, or those of its largest subtree
		vector<uint64_t> part_records;

		// STAGE 2: split the pairs of several subtrees into one file per subtree
		{
			RdbInitializer rdb_init;
			IntervUtils iu(_envir);
			vector<PairFiles> split;
			vector<int> split_idx;
			uint64_t num_split_records = 0;

			REprintf("Writing the track...\n");

			for (size_t ipair = 0; ipair < pairs.size(); ++ipair) {
				int64_t num_records = pairs[ipair].size / CONTACT_RECORD_SIZE;

				part_records.push_back(num_records);
				if (contacts_num_subtrees(iu, pairs[ipair]) > 1) {
					split.push_back(pairs[ipair]);
					split_idx.push_back(ipair);
					num_split_records += num_records;
				}
			}

			if (!split.empty()) {
				int num_kids = import_num_kids(iu, min((uint64_t)split.size(), num_split_records / misha::config::MIN_RECORDS_PER_PROCESS));
				vector<uint64_t> res = split_pairs(iu, split, num_kids, [&](const PairFiles &pair_files) {
						return split_pair(iu, dirname, pair_files);
					});

				for (size_t i = 0; i < split.size(); ++i)
					part_records[split_idx[i]] = res[i];
			}
			check_interrupt();
		}

		// STAGE 3: build the quad trees
		{
			RdbInitializer rdb_init;
			IntervUtils iu(_envir);
			bool allow_duplicates = LOGICAL(_allow_duplicates)[0];
			uint64_t num_records = 0;

			for (const PairFiles &pair_files : pairs)
				num_records += pair_files.size / CONTACT_RECORD_SIZE;

			int num_kids = import_num_kids(iu, min((uint64_t)pairs.size(), num_records / misha::config::MIN_RECORDS_PER_PROCESS));
			vector<uint64_t> mem;

			for (uint64_t num_part_records : part_records)
				mem.push_back(pair_mem(num_part_records, CONTACT_MEM));

			build_pairs(iu, pairs, mem, num_kids, [&](const PairFiles &pair_files) {
					write_pair(iu, dirname, pair_files, allow_duplicates);
				});

			// the chromosomes that got data, not all the names that matched
			report_skipped_chroms(srcname.c_str(), unknown, pairs_num_chroms(pairs), _envir);
		}
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
		rerror("Out of memory");
	} catch (const std::exception &e) {
		rerror("%s", e.what());
	}

	return R_NilValue;
}

}
