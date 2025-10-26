/*
 * GenomeTrackInfo.cpp
 *
 *  Created on: Feb 20, 2012
 *      Author: hoichman
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "rdbutils.h"
#include "GenomeTrack.h"
#include "GenomeTrackFixedBin.h"
#include "TrackIndex.h"

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrackinfo(SEXP _track, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		// check the arguments
		if (!Rf_isString(_track) || Rf_length(_track) != 1)
			verror("Track argument is not a string");

		const char *full_track_str = CHAR(STRING_ELT(_track, 0));

		enum { TYPE, DIM, SIZE, FORMAT, NUM_COLS };

		IntervUtils iu(_envir);
		SEXP answer, names, rtype, rdim, rsize, rformat;

		string trackpath(track2path(_envir, full_track_str));
		GenomeTrack::Type type = GenomeTrack::get_type(trackpath.c_str(), iu.get_chromkey());

		// Check if track uses indexed format
		string idx_path = trackpath + "/track.idx";
		struct stat idx_stat;
		bool is_indexed = (stat(idx_path.c_str(), &idx_stat) == 0);

		if (type == GenomeTrack::FIXED_BIN) {
			enum { BINSIZE = NUM_COLS, NUM_FIXED_BIN_COLS };

            answer = rprotect_ptr(RSaneAllocVector(VECSXP, NUM_FIXED_BIN_COLS));
            names = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_FIXED_BIN_COLS));

			GenomeTrackFixedBin gtrack;
			GIntervals all_genome_intervs;
			iu.get_all_genome_intervs(all_genome_intervs);
			char filename[FILENAME_MAX];
            SEXP rbinsize;

			string resolved = GenomeTrack::find_existing_1d_filename(iu.get_chromkey(), trackpath, all_genome_intervs.front().chromid);
			snprintf(filename, sizeof(filename), "%s/%s", trackpath.c_str(), resolved.c_str());

			gtrack.init_read(filename, all_genome_intervs.front().chromid);
            rbinsize = rprotect_ptr(Rf_ScalarInteger(gtrack.get_bin_size()));
			SET_VECTOR_ELT(answer, BINSIZE, rbinsize);
			SET_STRING_ELT(names, BINSIZE, Rf_mkChar("bin.size"));
		} else {
            answer = rprotect_ptr(RSaneAllocVector(VECSXP, NUM_COLS));
            names = rprotect_ptr(RSaneAllocVector(STRSXP, NUM_COLS));
		}

        rtype = rprotect_ptr(Rf_mkString(GenomeTrack::TYPE_NAMES[type]));
		SET_VECTOR_ELT(answer, TYPE, rtype);
		SET_STRING_ELT(names, TYPE, Rf_mkChar("type"));

        rdim = rprotect_ptr(Rf_ScalarInteger(GenomeTrack::is_1d(type) ? 1 : 2));
		SET_VECTOR_ELT(answer, DIM, rdim);
		SET_STRING_ELT(names, DIM, Rf_mkChar("dimensions"));

		int64_t totsize = 0;
		struct stat buf;

		if (is_indexed) {
			// For indexed tracks, sum track.dat + track.idx
			string dat_path = trackpath + "/track.dat";
			if (stat(dat_path.c_str(), &buf) == 0) {
				totsize += buf.st_size;
			}
			totsize += idx_stat.st_size;
		} else {
			// For per-chromosome tracks, sum all per-chromosome files
			vector<string> filenames;
			get_chrom_files(trackpath.c_str(), filenames);
			for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename) {
				string fullpath(trackpath + "/" + *ifilename);
				if (stat(fullpath.c_str(), &buf))
					verror("Cannot stat %s: %s", fullpath.c_str(), strerror(errno));
				totsize += buf.st_size;
			}
		}

        rsize = rprotect_ptr(Rf_ScalarReal(totsize));
		SET_VECTOR_ELT(answer, SIZE, rsize);
		SET_STRING_ELT(names, SIZE, Rf_mkChar("size.in.bytes"));

		// Add format field
		rformat = rprotect_ptr(Rf_mkString(is_indexed ? "indexed" : "per-chromosome"));
		SET_VECTOR_ELT(answer, FORMAT, rformat);
		SET_STRING_ELT(names, FORMAT, Rf_mkChar("format"));

		Rf_setAttrib(answer, R_NamesSymbol, names);

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
