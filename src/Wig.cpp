#include <errno.h>

#include "port.h"
#include "Wig.h"

string Wig::FIELDS_STRS[NUM_FIELDS] = { "chrom=", "start=", "step=", "span=" };

// Bytes read from the file at a time
static const size_t WIG_BLOCK_SIZE = 1 << 20;

// The characters that separate fields: isspace() of the C locale ('\n' ends a line, and '\r' is
// removed before a line is split)
static const struct SpaceChars {
	bool is[256];

	SpaceChars() {
		memset(is, 0, sizeof(is));
		is[(unsigned char)' '] = is[(unsigned char)'\t'] = is[(unsigned char)'\v'] = is[(unsigned char)'\f'] = true;
	}
} s_space;

void Wig::seek(int64_t pos)
{
	if (m_buf.empty())
		m_buf.resize(WIG_BLOCK_SIZE + 1);
	m_bfile.seek(pos, SEEK_SET);
	m_buf_fpos = pos;
	m_buf_pos = m_buf_end = 0;
	m_file_end = false;
}

// The next line, without its '\n', and whether it had one (the last line of a file may not).
// line[len] can be overwritten. Returns false at the end of the file.
bool Wig::next_line(char *&line, size_t &len, bool &terminated)
{
	while (1) {
		char *start = m_buf.data() + m_buf_pos;
		char *nl = (char *)memchr(start, '\n', m_buf_end - m_buf_pos);

		if (nl) {
			line = start;
			len = nl - start;
			terminated = true;
			m_buf_pos += len + 1;
			return true;
		}

		if (m_file_end) {
			if (m_buf_pos == m_buf_end)
				return false;
			line = start;
			len = m_buf_end - m_buf_pos;
			terminated = false;
			m_buf_pos = m_buf_end;
			return true;
		}

		// Keep the partial line and read more after it; one byte stays free for a NUL
		size_t left = m_buf_end - m_buf_pos;
		memmove(m_buf.data(), start, left);
		m_buf_fpos += m_buf_pos;
		m_buf_pos = 0;
		m_buf_end = left;
		if (m_buf.size() < left + WIG_BLOCK_SIZE + 1)
			m_buf.resize(left + WIG_BLOCK_SIZE + 1);

		uint64_t n = m_bfile.read(m_buf.data() + m_buf_end, WIG_BLOCK_SIZE);
		if (m_bfile.error())
			TGLError<Wig>(FILE_ERROR, "Failed to read WIG file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
		if (!n)
			m_file_end = true;
		m_buf_end += n;
	}
}

// Splits the next line that is not empty into m_fields, which is left empty at the end of the
// file. Returns the number of lines consumed, counted as split_line_by_space_chars counts them:
// every line, and one more when the end of the file is reached.
int Wig::split_fields()
{
	int num_lines = 0;
	char *line;
	size_t len;
	bool terminated;

	m_fields.clear();
	while (next_line(line, len, terminated)) {
		++num_lines;

		// '\r' is dropped wherever it is, not only before '\n'
		if (memchr(line, '\r', len)) {
			char *out = line;
			for (size_t i = 0; i < len; ++i) {
				if (line[i] != '\r')
					*out++ = line[i];
			}
			len = out - line;
		}

		char *p = line;
		char *end = line + len;
		bool trailing_space = len && s_space.is[(unsigned char)line[len - 1]];
		while (1) {
			while (p < end && s_space.is[(unsigned char)*p])
				++p;
			if (p == end)
				break;
			char *field = p;
			while (p < end && !s_space.is[(unsigned char)*p])
				++p;
			m_fields.push_back(Field{ field, (size_t)(p - field) });
			*p = '\0';  // a separator, the '\n' or the byte after the last line
			if (p < end)
				++p;
		}
		// As split_line_by_space_chars: space after the last field adds one empty field
		if (!m_fields.empty() && trailing_space) {
			*end = '\0';
			m_fields.push_back(Field{ end, 0 });
		}

		if (!m_fields.empty())
			return num_lines;
		if (!terminated)
			return num_lines;
	}
	return num_lines + 1;
}

void Wig::init(const GenomeChromKey &chromkey, const string &filename, bool ignore_unknown_chroms)
{
	m_chromkey = (GenomeChromKey *)&chromkey;
	m_ignore_unknown_chroms = ignore_unknown_chroms;
	m_chrom_fpos.clear();
	m_chrom_lineno.clear();
	m_chrom_fpos.resize(chromkey.get_num_chroms(), -1);
	m_chrom_lineno.resize(chromkey.get_num_chroms(), -1);
	m_unknown_chroms.clear();
	m_num_clipped = m_num_dropped = 0;
	m_first_clipped.clear();

	// read the file to determine the positions of first chrom appearences
	int last_chromid = -1;

	if (m_bfile.open(filename.c_str(), "r"))
		TGLError<Wig>(FILE_ERROR, "Failed to open WIG file %s: %s", m_bfile.file_name().c_str(), strerror(errno));
	m_has_last_chrom = false;
	seek(0);

	Rec rec;
	int64_t lineno = 1;
	int64_t last_lineno;
	long last_fpos;
	int last_rec_type = -1;
	int chrom_rec_type = -1;

	while (1) {
		last_fpos = tell();
		last_lineno = lineno;

		if (!read_record(rec, lineno))
			break;

		if (last_rec_type != -1 && rec.type == BEDGRAPH_REC && last_rec_type != BEDGRAPH_REC)
			TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: looks like a bizar mixture of WIG / BedGraph formats", m_bfile.file_name().c_str());

		if (rec.type == FIXED_STEP_REC || rec.type == VAR_STEP_REC || rec.type == BEDGRAPH_REC) {
			int chromid = -1;

			if (rec.type == FIXED_STEP_REC)
				chromid = rec.fixed_step.chromid;
			else if (rec.type == VAR_STEP_REC)
				chromid = rec.var_step.chromid;
			else
				chromid = rec.bedgraph.chromid;

			chrom_rec_type = rec.type;

			if (chromid == -1) { // unknown chrom
				// a block of a known chromosome that resumes after it is out of order, as after any other chromosome
				last_chromid = -1;
				continue;
			}

			if (chromid != last_chromid) {
				last_chromid = chromid;
				if (m_chrom_fpos[chromid] != -1)
					TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: file is not sorted by chromosomes", m_bfile.file_name().c_str());
				m_chrom_fpos[chromid] = last_fpos;
				m_chrom_lineno[chromid] = last_lineno;
			}
		} else if (rec.type == VAL_REC) {
			if (chrom_rec_type != FIXED_STEP_REC)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by fixedStep header",
						m_bfile.file_name().c_str(), lineno - 1);
		} else if (rec.type == COORD_VAL_REC) {
			if (chrom_rec_type != VAR_STEP_REC)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by variableStep header",
						m_bfile.file_name().c_str(), lineno - 1);
		}
	}
}

// a + b for a, b >= 0, capped at INT64_MAX: a record that far is past the chromosome end and dropped
static inline int64_t add_capped(int64_t a, int64_t b)
{
	return b > INT64_MAX - a ? INT64_MAX : a + b;
}

bool Wig::clip(int chromid, int64_t start, int64_t end, int64_t lineno)
{
	int64_t chromsize = m_chromkey->get_chrom_size(chromid);

	if (m_first_clipped.empty()) {
		char buf[1000];
		snprintf(buf, sizeof(buf), "line %lld: (%s, %lld, %lld), chromosome size %lld", (long long)lineno,
		         m_chromkey->id2chrom(chromid).c_str(), (long long)start, (long long)end, (long long)chromsize);
		m_first_clipped = buf;
	}
	if (start >= chromsize) {
		++m_num_dropped;
		return true;
	}
	++m_num_clipped;
	return false;
}

uint64_t Wig::get_num_matched_chroms() const
{
	uint64_t num_matched = 0;

	for (vector<long>::const_iterator ifpos = m_chrom_fpos.begin(); ifpos != m_chrom_fpos.end(); ++ifpos) {
		if (*ifpos != -1)
			++num_matched;
	}
	return num_matched;
}

bool Wig::get_data(int chromid, GIntervals &intervals)
{
	intervals.clear();

	if (m_chrom_fpos[chromid] == -1)
		return false;

	int64_t lineno = m_chrom_lineno[chromid];
	int64_t span = -1;
	int64_t start = -1;
	int64_t step = 0;
	Rec rec;

	seek(m_chrom_fpos[chromid]);

	while (1) {
		if (!read_record(rec, lineno))
			break;

		if (rec.type == FIXED_STEP_REC) {
			if (rec.fixed_step.chromid != chromid)
				break;

			start = rec.fixed_step.start;
			step = rec.fixed_step.step;
			span = rec.fixed_step.span;

			if (span > step)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: span exceeds step", m_bfile.file_name().c_str(), lineno - 1);
		} else if (rec.type == VAR_STEP_REC) {
			if (rec.var_step.chromid != chromid)
				break;

			start = -1;
			step = -1;
			span = rec.var_step.span;
		} else if (rec.type == VAL_REC) {
			if (start < 0)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by fixedStep header",
						m_bfile.file_name().c_str(), lineno - 1);

			int64_t chromsize = m_chromkey->get_chrom_size(chromid);
			int64_t end = add_capped(start, span);
			if (end <= chromsize || !clip(chromid, start, end, lineno - 1))
				intervals.push_back(GInterval(chromid, start, min(end, chromsize), 0, GInterval::cast2udata(rec.val)));
			start = add_capped(start, step);
		} else if (rec.type == COORD_VAL_REC) {
			if (span < 0)
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld: value is not preceeded by variableStep header",
						m_bfile.file_name().c_str(), lineno - 1);

			int64_t chromsize = m_chromkey->get_chrom_size(chromid);
			int64_t coord = rec.coord_val.coord;
			int64_t end = add_capped(coord, span);
			if (end <= chromsize || !clip(chromid, coord, end, lineno - 1))
				intervals.push_back(GInterval(chromid, coord, min(end, chromsize), 0, GInterval::cast2udata(rec.coord_val.val)));
		} else if (rec.type == BEDGRAPH_REC) {
			if (rec.bedgraph.chromid != chromid)
				break;

			if (rec.bedgraph.start >= rec.bedgraph.end)
				TGLError<Wig>(BAD_COORD, "WIG file %s, line %ld: start coordinate exceeds end coordinate", m_bfile.file_name().c_str(), lineno - 1);

			int64_t chromsize = m_chromkey->get_chrom_size(chromid);
			if (rec.bedgraph.end <= chromsize || !clip(chromid, rec.bedgraph.start, rec.bedgraph.end, lineno - 1))
				intervals.push_back(GInterval(chromid, rec.bedgraph.start, min(rec.bedgraph.end, chromsize), 0, GInterval::cast2udata(rec.bedgraph.val)));
		}
	}

	for (size_t i = 1; i < intervals.size(); ++i) {
		if (intervals[i - 1].end > intervals[i].start)
			TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s: file contains overlapping intervals", m_bfile.file_name().c_str());
	}

	return true;
}

bool Wig::read_record(Rec &rec, int64_t &lineno)
{
	vector<Field> &fields = m_fields;

	while (1) {
		lineno += split_fields();

		if (fields.empty())
			break;

		// UCSC header lines
		if (fields[0] == "track" || fields[0] == "browser" || fields[0].str[0] == '#')
			continue;

		if (fields[0] == "variableStep") {
			if (fields.size() < 2 || fields.size() > 3 ||
					!fields[1].starts_with(FIELDS_STRS[CHROM_FIELD]) ||
					(fields.size() == 3 && !fields[2].starts_with(FIELDS_STRS[SPAN_FIELD])))
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);

			rec.var_step.chromid = str2chromid(fields[1].str + FIELDS_STRS[CHROM_FIELD].length(), lineno - 1);
			rec.var_step.span = fields.size() == 3 ? str2span(fields[2].str + FIELDS_STRS[SPAN_FIELD].length(), lineno - 1) : 1;
			rec.type = VAR_STEP_REC;
			return true;
		}

		if (fields[0] == "fixedStep") {
			if (fields.size() < 4 || fields.size() > 5 ||
					!fields[1].starts_with(FIELDS_STRS[CHROM_FIELD]) ||
					!fields[2].starts_with(FIELDS_STRS[START_FIELD]) ||
					!fields[3].starts_with(FIELDS_STRS[STEP_FIELD]) ||
					(fields.size() == 5 && !fields[4].starts_with(FIELDS_STRS[SPAN_FIELD])))
				TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);

			rec.fixed_step.chromid = str2chromid(fields[1].str + FIELDS_STRS[CHROM_FIELD].length(), lineno - 1);
			rec.fixed_step.start = str2start(fields[2].str + FIELDS_STRS[START_FIELD].length(), lineno - 1);
			rec.fixed_step.step = str2step(fields[3].str + FIELDS_STRS[STEP_FIELD].length(), lineno - 1);
			rec.fixed_step.span = fields.size() == 5 ? str2span(fields[4].str + FIELDS_STRS[SPAN_FIELD].length(), lineno - 1) : 1;
			rec.type = FIXED_STEP_REC;
			return true;
		}

		if (fields.size() == 1) {
			rec.val = str2val(fields[0].str, lineno - 1);
			rec.type = VAL_REC;
			return true;
		}

		if (fields.size() == 2) {
			rec.coord_val.coord = str2coord(fields[0].str, 1, lineno - 1);
			rec.coord_val.val = str2val(fields[1].str, lineno - 1);
			rec.type = COORD_VAL_REC;
			return true;
		}

		if (fields.size() == 4) {
			// Treat as bedGraph regardless of chrom-name prefix. The
			// earlier "fields[0] starts with chr" gate spuriously
			// rejected bedGraph from genomes whose contig names don't
			// start with chr (Ensembl-style, non-mammalian, etc.).
			//
			// We still need to distinguish a real bedGraph row from a
			// malformed WIG line that happens to split into 4 tokens
			// (e.g. a value glued onto a fixedStep header by missing
			// newlines: "0.5fixedStep chrom=X start=1 step=1"). Such
			// lines carry "=" inside their fields; a real bedGraph row
			// never does. Use that as the cheap disambiguator: if any
			// field contains '=', fall through to the generic
			// "invalid format" error path below.
			bool looks_like_bedgraph = true;
			for (int i = 0; i < 4 && looks_like_bedgraph; ++i) {
				if (memchr(fields[i].str, '=', fields[i].len))
					looks_like_bedgraph = false;
			}
			if (looks_like_bedgraph) {
				rec.bedgraph.chromid = str2chromid(fields[0].str, lineno - 1);
				rec.bedgraph.start = str2coord(fields[1].str, 0, lineno - 1);
				rec.bedgraph.end = str2coord(fields[2].str, 0, lineno - 1);
				rec.bedgraph.val = str2val(fields[3].str, lineno - 1);
				rec.type = BEDGRAPH_REC;
				return true;
			}
		}

		TGLError<Wig>(INVALID_FORMAT, "Invalid format of WIG file %s, line %ld", m_bfile.file_name().c_str(), lineno - 1);
	}

	return false;
}
