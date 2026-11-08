#include <algorithm>
#include <memory>

#include "rdbutils.h"
#include "TrackNameIndex.h"

using namespace std;

TrackNameIndex::TrackNameIndex(SEXP names)
{
	if (!Rf_isString(names))
		return;

	int n = Rf_length(names);
	m_pos.reserve(n);
	for (int i = 0; i < n; ++i) {
		const char *name = CHAR(STRING_ELT(names, i));
		m_pos[name].push_back(i);
		bool identifier = true;
		for (const char *p = name; *p && identifier; ++p)
			identifier = rdb::is_R_var_char(*p);
		if (*name && !identifier)
			m_other.push_back(i);
	}
}

const TrackNameIndex &TrackNameIndex::gtracks(SEXP gtracks)
{
	// The indexed vector is kept from the GC while it is cached, so that another vector cannot be
	// allocated at its address. R code does not modify it in place: find_in_misha() marks what it
	// returns as not mutable.
	static unique_ptr<TrackNameIndex> s_index;
	static SEXP s_indexed = NULL;
	static const TrackNameIndex s_no_names(R_NilValue);

	if (!Rf_isString(gtracks))
		return s_no_names;
	if (!s_index || s_indexed != gtracks) {
		s_index.reset(new TrackNameIndex(gtracks));
		if (s_indexed)
			R_ReleaseObject(s_indexed);
		R_PreserveObject(gtracks);
		s_indexed = gtracks;
	}
	return *s_index;
}

const vector<int> *TrackNameIndex::find(const string &name) const
{
	auto it = m_pos.find(name);
	return it == m_pos.end() ? NULL : &it->second;
}

void TrackNameIndex::find_names(const string &expr, SEXP names, vector<int> &pos) const
{
	pos.clear();

	// A name made of R variable characters stands with no such character right before or after
	// it: it is a whole run of them.
	string run;
	for (size_t i = 0; i < expr.size(); ) {
		if (!rdb::is_R_var_char(expr[i])) {
			++i;
			continue;
		}
		size_t start = i;
		while (i < expr.size() && rdb::is_R_var_char(expr[i]))
			++i;
		run.assign(expr, start, i - start);
		const vector<int> *p = find(run);
		if (p)
			pos.insert(pos.end(), p->begin(), p->end());
	}

	for (int iname : m_other) {
		string name = CHAR(STRING_ELT(names, iname));
		uint64_t at = 0;

		while ((at = expr.find(name, at)) != string::npos) {
			uint64_t end = at + name.size();
			if ((!at || !rdb::is_R_var_char(expr[at - 1])) && (end == expr.size() || !rdb::is_R_var_char(expr[end]))) {
				pos.push_back(iname);
				break;
			}
			at += name.size();
		}
	}

	sort(pos.begin(), pos.end());
	pos.erase(unique(pos.begin(), pos.end()), pos.end());
}
