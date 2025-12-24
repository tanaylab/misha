#include <cstdint>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>

#include "GenomeTrackSparse.h"
#include "TrackIndex.h"

const int GenomeTrackSparse::RECORD_SIZE = 2 * sizeof(int64_t) + sizeof(float);
constexpr size_t GenomeTrackSparse::kSparseRecBytes;

GenomeTrackSparse::GenomeTrackSparse() :
	GenomeTrack1D(SPARSE),
	m_loaded(false),
	m_num_records(0),
	m_last_min_pos(numeric_limits<double>::quiet_NaN())
{}

void GenomeTrackSparse::read_header_at_current_pos_(BufferedFile &bf)
{
	int32_t format_signature = 0;
	if (bf.read(&format_signature, sizeof(format_signature)) != sizeof(format_signature))
		TGLError<GenomeTrackSparse>("Corrupt sparse header in %s", bf.file_name().c_str());
	if (format_signature >= 0)
		TGLError<GenomeTrackSparse>("Invalid sparse header signature in %s", bf.file_name().c_str());
}

void GenomeTrackSparse::init_read(const char *filename, int chromid)
{
	m_loaded = false; // Critical for read_file_into_mem()
	uint64_t header_start = 0;
	uint64_t total_bytes = 0;

	// Check for indexed format FIRST
	const std::string track_dir = GenomeTrack::get_track_dir(filename);
	const std::string idx_path = track_dir + "/track.idx";

	struct stat idx_st;
	if (stat(idx_path.c_str(), &idx_st) == 0) {
		// --- INDEXED PATH ---
		const std::string dat_path  = track_dir + "/track.dat";

		// "Smart Handle" logic - reopen if path or mode changed
		if (!m_dat_open || m_dat_path != dat_path || m_dat_mode != "rb") {
			m_bfile.close();
			if (m_bfile.open(dat_path.c_str(), "rb"))
				TGLError<GenomeTrackSparse>("Cannot open %s: %s", dat_path.c_str(), strerror(errno));
			m_dat_open = true;
			m_dat_path = dat_path;
			m_dat_mode = "rb";
		}

		auto idx   = get_track_index(track_dir);
		if (!idx)
			TGLError<GenomeTrackSparse>("Failed to load track index for %s", track_dir.c_str());

		auto entry = idx->get_entry(chromid);
		if (!entry || entry->length == 0) {
			// Chromosome not in index or empty contig - treat as empty
			m_num_records = 0;
			m_chromid = chromid;
			return;
		}

		if (m_bfile.seek(entry->offset, SEEK_SET))
			TGLError<GenomeTrackSparse>("Failed to seek to offset %llu in %s",
				(unsigned long long)entry->offset, dat_path.c_str());

		header_start = entry->offset;
		read_header_at_current_pos_(m_bfile);
		total_bytes = entry->length;
	} else {
		// --- PER-CHROMOSOME PATH ---
		m_bfile.close(); // Close previous file (if any)
		m_dat_open = false;
		read_type(filename); // This opens m_bfile and reads header

		header_start = 0;
		total_bytes = m_bfile.file_size();
	}

	// --- COMMON LOGIC ---
	const uint64_t header_size  = m_bfile.tell() - header_start;
	const double n = (total_bytes - header_size) / (double)kSparseRecBytes;

	if (n != (int64_t)n)
		TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s (n=%f, (int64_t)n=%lld)", filename, n, (long long)(int64_t)n);

	m_num_records = (int64_t)n;
	m_chromid = chromid;
}

void GenomeTrackSparse::init_write(const char *filename, int chromid)
{
	m_bfile.close();
	m_loaded = false;
	write_type(filename);
	m_chromid = chromid;
}

void GenomeTrackSparse::read_file_into_mem()
{
	if (m_loaded)
		return;

	m_intervals.resize(m_num_records);
	m_vals.resize(m_num_records);

	if (m_num_records == 0) {
		m_icur_interval = m_intervals.begin();
		m_loaded = true;
		return;
	}

	// Bulk read 
	const size_t total_bytes = m_num_records * kSparseRecBytes;
	std::vector<char> buffer(total_bytes);

	uint64_t bytes_read = m_bfile.read(buffer.data(), total_bytes);
	if (bytes_read != total_bytes) {
		if (m_bfile.error())
			TGLError<GenomeTrackSparse>("Failed to read a sparse track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s (expected %zu bytes, got %llu)",
			m_bfile.file_name().c_str(), total_bytes, (unsigned long long)bytes_read);
	}

	// Parse buffer into intervals and values
	const char *ptr = buffer.data();
	for (int64_t i = 0; i < m_num_records; ++i) {
		GInterval &interval = m_intervals[i];

		// Read start (int64_t)
		memcpy(&interval.start, ptr, sizeof(int64_t));
		ptr += sizeof(int64_t);

		// Read end (int64_t)
		memcpy(&interval.end, ptr, sizeof(int64_t));
		ptr += sizeof(int64_t);

		// Read val (float)
		memcpy(&m_vals[i], ptr, sizeof(float));
		ptr += sizeof(float);

		if (isinf(m_vals[i])) {
			m_vals[i] = numeric_limits<float>::quiet_NaN();
		}

		interval.chromid = m_chromid;

		if (interval.start < 0 || interval.start >= interval.end || (i && interval.start < m_intervals[i - 1].end)) {
			TGLError<GenomeTrackSparse>("Invalid format of a sparse track file %s", m_bfile.file_name().c_str());
		}
	}

	m_icur_interval = m_intervals.begin();
	m_loaded = true;
}

void GenomeTrackSparse::read_interval(const GInterval &interval)
{
	m_last_avg = m_last_nearest = m_last_min = m_last_max = m_last_stddev = m_last_sum = numeric_limits<float>::quiet_NaN();
	if (m_functions[MAX_POS])
		m_last_max_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[MIN_POS])
		m_last_min_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[EXISTS])
		m_last_exists = 0;
	if (m_functions[SIZE])
		m_last_size = 0;
	if (m_functions[SAMPLE])
		m_last_sample = numeric_limits<float>::quiet_NaN();
	if (m_functions[SAMPLE_POS])
		m_last_sample_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[FIRST])
		m_last_first = numeric_limits<float>::quiet_NaN();
	if (m_functions[FIRST_POS])
		m_last_first_pos = numeric_limits<double>::quiet_NaN();
	if (m_functions[LAST])
		m_last_last = numeric_limits<float>::quiet_NaN();
	if (m_functions[LAST_POS])
		m_last_last_pos = numeric_limits<double>::quiet_NaN();

	if (m_use_quantile)
		m_sp.reset();

	read_file_into_mem();

	if (m_intervals.empty())
		return;

	if (m_intervals.front().start >= interval.end) {
		m_last_nearest = m_vals.front();
		return;
	}

	if (m_intervals.back().end <= interval.start) {
		m_last_nearest = m_vals.back();
		return;
	}

	if (check_first_overlap(m_icur_interval, interval)) {
		calc_vals(interval);
	} else if (m_icur_interval + 1 < m_intervals.end() && check_first_overlap(m_icur_interval + 1, interval)) {
		++m_icur_interval;
		calc_vals(interval);
	} else {
		// run the binary search
		GIntervals::const_iterator istart_interval = m_intervals.begin();
		GIntervals::const_iterator iend_interval = m_intervals.end();

		while (iend_interval - istart_interval > 1) {
			GIntervals::const_iterator imid_interval = istart_interval + (iend_interval - istart_interval) / 2;

			if (check_first_overlap(imid_interval, interval)) {
				m_icur_interval = imid_interval;
				calc_vals(interval);
				break;
			}

			// is mid_interval < interval?
			if (GIntervals::compare_by_start_coord(*imid_interval, interval))
				istart_interval = imid_interval;
			else
				iend_interval = imid_interval;
		}

		if (iend_interval - istart_interval == 1 && check_first_overlap(istart_interval, interval)) {
			m_icur_interval = istart_interval;
			calc_vals(interval);
		}

		if (iend_interval - istart_interval == 1)
			m_last_nearest = iend_interval == m_intervals.end() || interval.dist2interv(*istart_interval) <= interval.dist2interv(*iend_interval) ?
					m_vals[istart_interval - m_intervals.begin()] : m_vals[iend_interval - m_intervals.begin()];
	}
}

double GenomeTrackSparse::last_max_pos() const
{
	return m_last_max_pos;
}

double GenomeTrackSparse::last_min_pos() const
{
	return m_last_min_pos;
}

void GenomeTrackSparse::write_next_interval(const GInterval &interval, float val)
{
	uint64_t size = 0;
	size += m_bfile.write(&interval.start, sizeof(interval.start));
	size += m_bfile.write(&interval.end, sizeof(interval.end));
	size += m_bfile.write(&val, sizeof(val));

	if ((int)size != RECORD_SIZE) {
		if (m_bfile.error())
			TGLError<GenomeTrackSparse>("Failed to write a sparse track file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		TGLError<GenomeTrackSparse>("Failed to write a sparse track file %s", m_bfile.file_name().c_str());
	}
}

const GIntervals &GenomeTrackSparse::get_intervals()
{
	read_file_into_mem();
	return m_intervals;
}

const vector<float> &GenomeTrackSparse::get_vals()
{
	read_file_into_mem();
	return m_vals;
}
