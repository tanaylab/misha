#ifndef GTRACK2DIMPORT_H_
#define GTRACK2DIMPORT_H_

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// Shared by gtrack_2d_import and gtrack_import_contacts. Both run in two stages:
//   1. read the input files and split the records into intermediate binary files, one per
//      chromosome pair and kid;
//   2. build the quad tree of each chromosome pair from its intermediate files.
// With multitasking on, stage 1 runs in kids that each read a contiguous block of the input files,
// and stage 2 in kids that each build the quad trees of a subset of the chromosome pairs, one pair
// at a time. A pair's intermediate files are read back in kid order, which is input file order,
// so every quad tree gets its records in the same order as a serial run and the track files come
// out byte-identical.

namespace rdb { class IntervUtils; }

struct PairFiles {
	int                      chromid1;
	int                      chromid2;
	int64_t                  size;   // total size of the intermediate files in bytes
	std::vector<std::string> files;  // intermediate files ordered by kid
};

// Number of kids to split num_tasks tasks between: 1 when multitasking is off
int import_num_kids(const rdb::IntervUtils &iu, uint64_t num_tasks);

// Splits the input files into at most num_kids contiguous blocks balanced by size.
// Block k is files [bounds[k], bounds[k + 1]); the number of blocks is bounds.size() - 1.
std::vector<int> split_files(const std::vector<int64_t> &sizes, int num_kids);

// Intermediate file of a chromosome pair written by a kid in stage 1
std::string pair_file_name(const std::string &dirname, int chromid1, int chromid2, int kid);

// Collects the intermediate files that stage 1 left in dirname, ordered by chromosome pair
std::vector<PairFiles> list_pair_files(const std::string &dirname);

// Assigns the pairs to num_kids kids: largest first, each to the least loaded kid.
// Returns the indices of the pairs of each kid.
std::vector<std::vector<int>> assign_pairs(const std::vector<PairFiles> &pairs, int num_kids);

// Runs work(kid) for each kid - in forked kids when num_kids > 1, in this process otherwise -
// and returns the byte each work() returned.
std::vector<char> run_kids(rdb::IntervUtils &iu, int num_kids, const std::function<char(int)> &work);

#endif
