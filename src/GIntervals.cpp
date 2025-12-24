#include <cstdint>
#include <algorithm>

#include "GIntervals.h"

void GIntervals::sort(bool (*cmp_function)(const GInterval &, const GInterval &))
{
	// if the intervals are empty do not sort
	if (empty()){
		return;
	}

	// Only sort if not already sorted (std::is_sorted may have better optimizations)
	if (!std::is_sorted(begin(), end(), cmp_function)) {
		std::sort(begin(), end(), cmp_function);
	}
}

// intervs are expected to be already sorted
void GIntervals::verify_no_overlaps(const GenomeChromKey &chromkey, const char *error_prefix) const
{
	for (const_iterator iinterv = begin() + 1; iinterv < end(); ++iinterv) {
		if (*iinterv < *(iinterv - 1))
			TGLError<GIntervalsFetcher1D>(UNSORTED_INTERVALS, "%sTo verify overlaps intervals must be sorted", error_prefix);

		if (iinterv->chromid == (iinterv - 1)->chromid && ((iinterv - 1)->end > iinterv->start))
			TGLError<GIntervalsFetcher1D>(OVERLAPPING_INTERVAL, "%sIntervals (%s, %ld, %ld) and (%s, %ld, %ld) overlap",
										  error_prefix,
										  chromkey.id2chrom((iinterv - 1)->chromid).c_str(), (iinterv - 1)->start, (iinterv - 1)->end,
										  chromkey.id2chrom(iinterv->chromid).c_str(), iinterv->start, iinterv->end);
	}
}

int64_t GIntervals::range() const
{
	int64_t range = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv)
		range += iinterv->range();
	return range;
}

int64_t GIntervals::range(int chromid) const
{
	int64_t range = 0;

	for (const_iterator iinterv = begin(); iinterv < end(); ++iinterv) {
		if (iinterv->chromid == chromid) 
			range += iinterv->range();
	}
	return range;
}

// intervs are expected to be already sorted
void GIntervals::unify_overlaps(bool unify_touching_intervals)
{
	if (empty())
		return;

	uint64_t cur_idx = 0;

	for (uint64_t i = 1; i < size(); i++) {
		if (operator[](cur_idx).chromid != operator[](i).chromid || operator[](cur_idx).end < operator[](i).start || (!unify_touching_intervals && operator[](cur_idx).end == operator[](i).start))
			operator[](++cur_idx) = operator[](i);
		// unite overlapping intervals
		else if (operator[](cur_idx).end < operator[](i).end)
			operator[](cur_idx).end = operator[](i).end;
	}
	erase(begin() + cur_idx + 1, end());
}

void GIntervals::unify(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs)
{
	const_iterator iintervs[] = { intervs1.begin(), intervs2.begin() };
	const_iterator intervends[] = { intervs1.end(), intervs2.end() };
	int last_chromid[] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	res_intervs.reserve(intervs1.size() + intervs2.size());

	// merge the two intervs into one vector sorted by chrom and start coord (complexity: O(n))
	while (iintervs[0] != intervends[0] && iintervs[1] != intervends[1]) {
		if (iintervs[0]->chromid == iintervs[1]->chromid)
			idx = iintervs[0]->start < iintervs[1]->start ? 0 : 1;
		else if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
			idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
			last_chromid[0] = iintervs[0]->chromid;
			last_chromid[1] = iintervs[1]->chromid;
		}

		res_intervs.push_back(*iintervs[idx]);
		++iintervs[idx];
	}

	for (int i = 0; i < 2; i++) {
		for (const_iterator iinterv = iintervs[i]; iinterv != intervends[i]; ++iinterv)
			res_intervs.push_back(*iinterv);
	}

	// having the whole set of merged (i.e. sorted) intervals let's unify the overlapping intervals (complexity: O(n))
	res_intervs.unify_overlaps();
}

void GIntervals::intersect(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs)
{
	// use local state to track virtual start positions
	// instead of copying entire input vectors
	const_iterator iintervs[2] = { intervs1.begin(), intervs2.begin() };
	const_iterator intervends[2] = { intervs1.end(), intervs2.end() };
	// Virtual start positions - allows "modifying" interval starts without copying
	int64_t virt_start[2] = { 0, 0 };
	bool use_virt[2] = { false, false };
	int last_chromid[2] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	while (iintervs[0] != intervends[0] && iintervs[1] != intervends[1]) {
		// Get effective start positions (use virtual if set, otherwise actual)
		int64_t eff_start0 = use_virt[0] ? virt_start[0] : iintervs[0]->start;
		int64_t eff_start1 = use_virt[1] ? virt_start[1] : iintervs[1]->start;

		if (iintervs[0]->chromid == iintervs[1]->chromid) {
			if (eff_start0 < eff_start1 && iintervs[0]->end <= eff_start1) {
				++iintervs[0];
				use_virt[0] = false;
			}
			else if (eff_start1 < eff_start0 && iintervs[1]->end <= eff_start0) {
				++iintervs[1];
				use_virt[1] = false;
			}
			else { // intervals intersect
				int64_t start = max(eff_start0, eff_start1);
				int64_t end = min(iintervs[0]->end, iintervs[1]->end);

				res_intervs.push_back(GInterval(iintervs[0]->chromid, start, end, 0));
				for (int i = 0; i < 2; i++) {
					if (iintervs[i]->end == end) {
						++iintervs[i];
						use_virt[i] = false;
					} else {
						virt_start[i] = end;
						use_virt[i] = true;
					}
				}
			}
		} else {
			if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
				idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
				last_chromid[0] = iintervs[0]->chromid;
				last_chromid[1] = iintervs[1]->chromid;
			} else {
				++iintervs[idx];
				use_virt[idx] = false;
			}
		}
	}
}

void GIntervals::diff(const GIntervals &intervs1, const GIntervals &intervs2, GIntervals &res_intervs)
{
	// Optimized version: use local state to track virtual start positions
	// instead of copying entire input vectors
	const_iterator iintervs[2] = { intervs1.begin(), intervs2.begin() };
	const_iterator intervends[2] = { intervs1.end(), intervs2.end() };
	// Virtual start positions - allows "modifying" interval starts without copying
	int64_t virt_start[2] = { 0, 0 };
	bool use_virt[2] = { false, false };
	int last_chromid[2] = { -1, -1 };
	int idx = 0;

	res_intervs.clear();
	while (iintervs[0] != intervends[0] && iintervs[1] != intervends[1]) {
		// Get effective start positions (use virtual if set, otherwise actual)
		int64_t eff_start0 = use_virt[0] ? virt_start[0] : iintervs[0]->start;
		int64_t eff_start1 = use_virt[1] ? virt_start[1] : iintervs[1]->start;

		if (iintervs[0]->chromid == iintervs[1]->chromid) {
			if (eff_start0 < eff_start1 && iintervs[0]->end <= eff_start1) {
				res_intervs.push_back(GInterval(iintervs[0]->chromid, eff_start0, iintervs[0]->end, 0));
				++iintervs[0];
				use_virt[0] = false;
			}
			else if (eff_start1 < eff_start0 && iintervs[1]->end <= eff_start0) {
				++iintervs[1];
				use_virt[1] = false;
			}
			else { // intervals intersect
				int64_t intersect_start = max(eff_start0, eff_start1);
				int64_t intersect_end = min(iintervs[0]->end, iintervs[1]->end);

				if (eff_start0 < intersect_start)
					res_intervs.push_back(GInterval(iintervs[0]->chromid, eff_start0, intersect_start, 0));

				for (int i = 0; i < 2; i++) {
					if (iintervs[i]->end == intersect_end) {
						++iintervs[i];
						use_virt[i] = false;
					} else {
						virt_start[i] = intersect_end;
						use_virt[i] = true;
					}
				}
			}
		} else {
			if (last_chromid[0] != iintervs[0]->chromid || last_chromid[1] != iintervs[1]->chromid) {
				idx = compare_by_start_coord(*iintervs[0], *iintervs[1]) ? 0 : 1;
				last_chromid[0] = iintervs[0]->chromid;
				last_chromid[1] = iintervs[1]->chromid;
			} else {
				if (!idx) {
					int64_t eff_start = use_virt[0] ? virt_start[0] : iintervs[0]->start;
					res_intervs.push_back(GInterval(iintervs[0]->chromid, eff_start, iintervs[0]->end, 0));
				}
				++iintervs[idx];
				use_virt[idx] = false;
			}
		}
	}

	// Append remaining intervals from intervs1 (with virtual start if applicable)
	if (iintervs[0] != intervends[0]) {
		if (use_virt[0]) {
			res_intervs.push_back(GInterval(iintervs[0]->chromid, virt_start[0], iintervs[0]->end, 0));
			++iintervs[0];
		}
		for (const_iterator iinterv = iintervs[0]; iinterv != intervends[0]; ++iinterv)
			res_intervs.push_back(*iinterv);
	}
}

const GInterval *GIntervals::containing_interval(const GInterval &interv)
{
	// run binary search
	GIntervals::const_iterator istart_interval = begin();
	GIntervals::const_iterator iend_interval = end();

	while (iend_interval - istart_interval > 1) {
		GIntervals::const_iterator imid_interval = istart_interval + (iend_interval - istart_interval) / 2;

		if (imid_interval->do_overlap(interv))
			return imid_interval->do_contain(interv) ? &*imid_interval : NULL;

		// is mid_interval < interval?
		if (GIntervals::compare_by_start_coord(*imid_interval, interv))
			istart_interval = imid_interval;
		else
			iend_interval = imid_interval;
	}

	if (iend_interval - istart_interval == 1 && istart_interval->do_overlap(interv))
		return istart_interval->do_contain(interv) ? &*istart_interval : NULL;

	return NULL;
}

void GIntervals::read(const GenomeChromKey &chromkey, istream &tab, int nostrand)
{
	string chrom;
	int64_t start, end;
	int strand;

	strand = 1;
	tab >> chrom;
	while(tab) {
		tab >> start >> end;
		if(!nostrand) {
			tab >> strand;
		}
		GInterval interval(chromkey.chrom2id(chrom.c_str()), start, end, (char)strand);
		interval.verify(chromkey);
		push_back(interval);
		tab >> chrom;
	}
}

void GIntervals::read_bed(const GenomeChromKey &chromkey, istream &bed)
{
	string chrom;
	int64_t start, end;
	char strandcode;
	float score;
	string name;
	int strand = 0;

	bed >> chrom;
	while(bed) {
		bed >> start >> end >> name >> score >> strandcode;
		try {
			strand = GInterval::char2strand(strandcode);
		} catch (TGLException &e) {
			TGLError<GInterval>(GInterval::BAD_STRAND, "Reading interval (%s, %ld, %ld, %c): %s", chrom.c_str(), start, end, strandcode, e.msg());
		}
		GInterval interval(chromkey.chrom2id(chrom.c_str()), start, end, (char)strand);
		interval.verify(chromkey);
		push_back(interval);
		while(bed.get() != '\n') {}
		bed >> chrom;
	}
}

void GIntervals::begin_chrom_iter(int chromid)
{
	build_chrom_map();
	m_cur_chromid = chromid;
	m_iter_chrom_index = 0;
	if ((uint64_t)chromid < m_chrom2itr.size())
		m_iinterval = m_chrom2itr[chromid];
	else
		m_iinterval = end();
}

GIntervals::const_iterator GIntervals::get_chrom_begin() const
{
	build_chrom_map();
	return m_chrom2itr[m_iinterval->chromid];
}

GIntervals::const_iterator GIntervals::get_chrom_end() const
{
	build_chrom_map();
	return (uint64_t)m_iinterval->chromid + 1 < m_chrom2itr.size() ? m_chrom2itr[m_iinterval->chromid + 1] : end();
}

void GIntervals::write(const GenomeChromKey &chromkey, ostream &tab)
{
	for(const_iterator i = begin(); i != end(); i++) {
		tab << chromkey.id2chrom(i->chromid)
		<< "\t" << i->start
		<< "\t" << i->end
		<< "\t" << (int)i->strand << "\n";
	}
}
