#ifndef TRACKNAMEINDEX_H_
#define TRACKNAMEINDEX_H_

#include <string>
#include <unordered_map>
#include <vector>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

// The names of a character vector (GTRACKS, or the virtual track names), by name. A database can
// hold ~10^5 tracks and every track expression is searched for the names it uses: with the index
// that is one lookup per identifier of the expression (find_names) instead of a substring search
// of the expression per name.
class TrackNameIndex {
public:
	// Index of `names` (anything but a character vector has no names)
	explicit TrackNameIndex(SEXP names);

	// The index of GTRACKS `gtracks`, kept between calls: rebuilt only when GTRACKS is another
	// vector than the one indexed (it is replaced, not modified, when the track list changes).
	static const TrackNameIndex &gtracks(SEXP gtracks);

	// Positions (0-based) of `name` in the vector, or NULL
	const std::vector<int> *find(const std::string &name) const;

	// Positions of the names that occur in `expr`, each once, ascending: a name occurs where it
	// stands with no R variable character (rdb::is_R_var_char) right before or after it, inside
	// a string or not. The positions of the vector's elements are `names`'s.
	void find_names(const std::string &expr, SEXP names, std::vector<int> &pos) const;

private:
	std::unordered_map<std::string, std::vector<int>> m_pos;
	// names with a character other than an R variable character, which cannot be found as a
	// whole identifier: the expression is searched for them as substrings (not for an empty name)
	std::vector<int> m_other;
};

#endif /* TRACKNAMEINDEX_H_ */
