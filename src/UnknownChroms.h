#ifndef UNKNOWNCHROMS_H_
#define UNKNOWNCHROMS_H_

#include <cctype>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------
// Chromosome names that a file mentioned and the genome database does not have.
//
// Aliases are resolved by GenomeChromKey::chrom2id before a name can get here,
// so a name the alias chain bridges (chr1 <-> 1, M <-> MT, chrom_aliases.tsv)
// never lands in this collection.
//
// Two kinds of unmatched name mean very different things to a caller:
//   - "chr1_gl000191_random", "GL000220.1", "chrUn_KI270302v1": a contig the
//     database deliberately does not have. Everyday, and no alias can fix it.
//   - "chr7", "X": a primary chromosome. The database should have had it, so the
//     naming is probably wrong.
// They are collected separately so each can be reported in its own channel.
//
// Mitochondria are the first kind, not the second, even though "MT" looks like a
// chromosome name: .compute_chrom_aliases auto-generates every mito spelling
// (M / MT / chrM) for whatever mito contig a database has, so a mito name that
// reached here proves the database has no mito contig at all. That can never be
// fixed by an alias, and hg19 / mm10 as built here have no chrM - so treating it
// as primary would put every whole-genome bigWig back in the warning channel.

class UnknownChroms {
public:
	// MAX_TRACKED: distinct names are tracked up to this many; past it the count is
	//   reported as "N+". Bounds the memory a mis-columned file (first column = read
	//   id, say) can make us allocate.
	// MAX_REPORTED: how many names are kept for the message text.
	// (An enum rather than static const members so no out-of-line definition is
	// needed when they are passed by reference, e.g. to std::min.)
	enum { MAX_TRACKED = 100, MAX_REPORTED = 5 };

	UnknownChroms() : m_truncated(false) {}

	void clear() {
		m_names.clear();
		m_primary_names.clear();
		m_seen.clear();
		m_primary_seen.clear();
		m_truncated = false;
	}

	void record(const string &chrom) {
		// Primary-shaped names are tracked without a cap, deliberately: the cap exists to
		// bound a mis-columned file's junk names, and there are only so many names that can
		// look like a chromosome. Capping them would let a chr7 that appears after 100
		// scaffolds go unnoticed and silently downgrade the report to a message.
		if (is_primary_chrom_name(chrom) && m_primary_seen.insert(chrom).second &&
		    m_primary_names.size() < (size_t)MAX_REPORTED)
			m_primary_names.push_back(chrom);

		if (m_seen.size() >= (size_t)MAX_TRACKED) {
			// only a name that is not already tracked is actually being dropped
			if (m_seen.find(chrom) == m_seen.end())
				m_truncated = true;
			return;
		}

		if (!m_seen.insert(chrom).second)
			return;

		if (m_names.size() < (size_t)MAX_REPORTED)
			m_names.push_back(chrom);
	}

	bool empty() const { return m_seen.empty() && m_primary_seen.empty(); }

	// Number of distinct names; capped at MAX_TRACKED, in which case truncated() is true.
	uint64_t num() const { return (uint64_t)m_seen.size(); }
	bool truncated() const { return m_truncated; }

	// Up to MAX_REPORTED names, in order of appearance.
	const vector<string> &names() const { return m_names; }

	// Up to MAX_REPORTED of the names that look like a primary chromosome, and the true
	// number of distinct such names.
	const vector<string> &primary_names() const { return m_primary_names; }
	uint64_t num_primary() const { return (uint64_t)m_primary_seen.size(); }

	// "chr7", "7", "X" (any case, with or without the chr prefix) - the names a genome
	// database is expected to have. Scaffolds, patches and unplaced contigs carry extra
	// tokens and do not match; neither does a mitochondrial name (see above).
	//
	// Allocation-free: this runs once per unmatched record, not once per distinct name.
	static bool is_primary_chrom_name(const string &name) {
		size_t off = 0;

		if (name.size() > 3 && tolower((unsigned char)name[0]) == 'c' &&
		    tolower((unsigned char)name[1]) == 'h' && tolower((unsigned char)name[2]) == 'r')
			off = 3;

		size_t len = name.size() - off;

		if (!len || len > 2)
			return false;

		bool all_digits = true;
		for (size_t i = off; i < name.size(); ++i) {
			if (!isdigit((unsigned char)name[i])) {
				all_digits = false;
				break;
			}
		}
		if (all_digits)
			return true;

		if (len > 1)
			return false;

		char c = (char)toupper((unsigned char)name[off]);
		return c == 'X' || c == 'Y' || c == 'Z' || c == 'W';
	}

private:
	vector<string> m_names;
	vector<string> m_primary_names;
	set<string>    m_seen;
	set<string>    m_primary_seen;
	bool           m_truncated;
};

#endif /* UNKNOWNCHROMS_H_ */
