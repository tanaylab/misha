#ifndef GTRACK2DIMPORT_H_
#define GTRACK2DIMPORT_H_

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// Shared by gtrack_2d_import and gtrack_import_contacts. Both run in three stages:
//   1. read the input files and split the records into intermediate binary files, one per
//      chromosome pair and kid;
//   2. split each chromosome pair whose quad tree has several subtrees (see pair_num_subtrees)
//      into one file per subtree, so that its memory can be estimated by its largest subtree;
//   3. build the quad tree of each chromosome pair from its intermediate or subtree files.
// With multitasking on, stage 1 runs in kids that each read a contiguous block of the input files,
// and stages 2 and 3 in kids that each take the next chromosome pair, one pair at a time. A pair's
// intermediate files are read back in kid order, which is input file order, so every quad tree gets
// its records in the same order as a serial run and the track files come out byte-identical,
// whichever kid splits or builds the pair.

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

// Runs work(kid) for each kid - in forked kids when num_kids > 1, in this process otherwise -
// and returns the byte each work() returned. With throttle, a kid is suspended while the kids together
// use more than gmax.mem.usage (see RdbInitializer::report_alloc).
std::vector<char> run_kids(rdb::IntervUtils &iu, int num_kids, const std::function<char(int)> &work, bool throttle = true);

// Number of subtrees of the quad tree of a chromosome pair of num_records records: num_records / gmax.data.size
// rounded down to a power of 4, at least 1
int64_t pair_num_subtrees(const rdb::IntervUtils &iu, int64_t num_records);

// Name of the file that holds the records of subtree isubtree of a chromosome pair (stage 2)
std::string subtree_file_name(const std::string &dirname, int chromid1, int chromid2, int64_t isubtree);

// Stage 2: runs split(pair) for every pair in num_kids kids (see run_kids), largest first, and returns what
// each call returned, in the order of pairs
std::vector<uint64_t> split_pairs(rdb::IntervUtils &iu, const std::vector<PairFiles> &pairs, int num_kids,
								  const std::function<uint64_t(const PairFiles &)> &split);

// Estimated peak memory in bytes of building a chromosome pair whose largest subtree holds num_records records,
// bytes_per_record each. The serializer holds one subtree at a time.
uint64_t pair_mem(uint64_t num_records, uint64_t bytes_per_record);

// Stage 3: runs build(pair) for every pair in num_kids kids (see run_kids). Each kid takes the next pair,
// largest first, and a pair starts only when its estimated memory mem[pair] fits into gmax.mem.usage
// together with the estimates of the pairs being built. A smaller pair that fits may start while a larger
// one waits. A pair whose estimate exceeds gmax.mem.usage is built alone.
// The throttle of run_kids is off here: it keeps one kid running and suspends the others, and if the kid
// it keeps is one waiting for memory, the suspended kids never resume.
void build_pairs(rdb::IntervUtils &iu, const std::vector<PairFiles> &pairs, const std::vector<uint64_t> &mem, int num_kids,
				 const std::function<void(const PairFiles &)> &build);

#endif
