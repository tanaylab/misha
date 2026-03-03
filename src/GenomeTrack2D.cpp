#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "TGLException.h"
#include "GenomeTrack2D.h"
#include "TrackIndex2D.h"

void GenomeTrack2D::init_read(const char *filename, int chromid1, int chromid2)
{
	m_loaded = false;
	m_pair_has_data = false;
	m_base_offset_2d = 0;

	// Check for indexed format FIRST
	const std::string track_dir = GenomeTrack::get_track_dir(filename);
	const std::string idx_path = track_dir + "/track.idx";

	struct stat idx_st;
	if (stat(idx_path.c_str(), &idx_st) == 0) {
		// Try to load as a 2D index
		std::shared_ptr<TrackIndex2D> idx;
		try {
			idx = TrackIndex2D::get_track_index_2d(track_dir);
		} catch (...) {
			idx = nullptr;
		}

		if (idx && idx->is_loaded()) {
			// --- INDEXED PATH ---
			const std::string dat_path = track_dir + "/track.dat";

			// Smart handle: reopen file only if path changed
			if (!m_dat_open || m_dat_path != dat_path) {
				m_bfile.close();
				if (m_bfile.open(dat_path.c_str(), "rb"))
					TGLError<GenomeTrack2D>("Cannot open %s: %s", dat_path.c_str(), strerror(errno));
				m_dat_open = true;
				m_dat_path = dat_path;
			}

			auto entry = idx->get_entry(chromid1, chromid2);
			if (!entry || entry->length == 0) {
				// Pair not in index or empty -- no data for this pair
				m_chromid1 = chromid1;
				m_chromid2 = chromid2;
				m_pair_has_data = false;
				// Do NOT close m_bfile -- we keep it open for the next pair
				return;
			}

			// Seek to the start of this pair's data
			if (m_bfile.seek(entry->offset, SEEK_SET))
				TGLError<GenomeTrack2D>("Failed to seek to offset %llu in %s",
					(unsigned long long)entry->offset, dat_path.c_str());

			// Read and validate the format signature
			int format_signature;
			if (m_bfile.read(&format_signature, sizeof(format_signature)) != sizeof(format_signature))
				TGLError<GenomeTrack2D>("Failed to read format signature in %s at offset %llu",
					dat_path.c_str(), (unsigned long long)entry->offset);

			// Validate type matches
			if (format_signature != FORMAT_SIGNATURES[m_type])
				TGLError<GenomeTrack2D>("Track type mismatch in indexed 2D track at %s (expected %d, got %d)",
					dat_path.c_str(), FORMAT_SIGNATURES[m_type], format_signature);

			m_chromid1 = chromid1;
			m_chromid2 = chromid2;
			m_pair_has_data = true;
			m_base_offset_2d = entry->offset;
			// File position is now right after the format signature,
			// ready for the subclass (GenomeTrackRects) to call unserialize()
			return;
		}
		// If idx failed to load as 2D index, fall through to per-pair path
		// (the track.idx might be a 1D index -- different magic header)
	}

	// --- PER-PAIR PATH (existing logic) ---
	m_bfile.close();
	m_dat_open = false;

	if (!access(filename, R_OK) || errno != ENOENT) {
		read_type(filename);
		m_pair_has_data = m_bfile.opened();
	} else {
		m_pair_has_data = false;
	}

	m_chromid1 = chromid1;
	m_chromid2 = chromid2;
}

void GenomeTrack2D::init_write(const char *filename, int chromid1, int chromid2)
{
	m_bfile.close();
	m_loaded = false;
	m_dat_open = false;
	m_pair_has_data = false;
	write_type(filename);
	m_chromid1 = chromid1;
	m_chromid2 = chromid2;
}
