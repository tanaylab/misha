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
//   - "chr7", "X", "MT": a primary chromosome. The database should have had it,
//     so the naming is probably wrong.
// They are collected separately so each can be reported in its own channel.

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
		m_truncated = false;
	}

	void record(const string &chrom) {
		if (m_seen.size() >= (size_t)MAX_TRACKED) {
			m_truncated = true;
			return;
		}

		if (!m_seen.insert(chrom).second)
			return;

		if (m_names.size() < (size_t)MAX_REPORTED)
			m_names.push_back(chrom);

		if (is_primary_chrom_name(chrom) && m_primary_names.size() < (size_t)MAX_REPORTED)
			m_primary_names.push_back(chrom);
	}

	bool empty() const { return m_seen.empty(); }

	// Number of distinct names; equals MAX_TRACKED when truncated() is true.
	uint64_t num() const { return (uint64_t)m_seen.size(); }
	bool truncated() const { return m_truncated; }

	// Up to MAX_REPORTED names, in order of appearance.
	const vector<string> &names() const { return m_names; }

	// Up to MAX_REPORTED of those names that look like a primary chromosome.
	const vector<string> &primary_names() const { return m_primary_names; }

	// "chr7", "7", "X", "MT" (any case, with or without the chr prefix) - the names a
	// genome database is expected to have. Scaffolds, patches and unplaced contigs
	// carry extra tokens and do not match.
	static bool is_primary_chrom_name(const string &name) {
		string s = name;

		if (s.size() > 3) {
			string prefix = s.substr(0, 3);
			for (size_t i = 0; i < prefix.size(); ++i)
				prefix[i] = (char)tolower((unsigned char)prefix[i]);
			if (prefix == "chr")
				s = s.substr(3);
		}

		if (s.empty() || s.size() > 2)
			return false;

		bool all_digits = true;
		for (size_t i = 0; i < s.size(); ++i) {
			if (!isdigit((unsigned char)s[i])) {
				all_digits = false;
				break;
			}
		}
		if (all_digits)
			return true;

		for (size_t i = 0; i < s.size(); ++i)
			s[i] = (char)toupper((unsigned char)s[i]);
		return s == "X" || s == "Y" || s == "Z" || s == "W" || s == "M" || s == "MT";
	}

private:
	vector<string> m_names;
	vector<string> m_primary_names;
	set<string>    m_seen;
	bool           m_truncated;
};

#endif /* UNKNOWNCHROMS_H_ */
