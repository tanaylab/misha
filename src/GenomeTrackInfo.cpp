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

using namespace std;
using namespace rdb;

extern "C" {

SEXP gtrackinfo(SEXP _track, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		// check the arguments
		if (!isString(_track) || length(_track) != 1)
			verror("Track argument is not a string");

		const char *full_track_str = CHAR(STRING_ELT(_track, 0));

		enum { TYPE, DIM, SIZE, NUM_COLS };

		IntervUtils iu(_envir);
		SEXP answer;
		SEXP names;

		string trackpath(track2path(_envir, full_track_str));
		GenomeTrack::Type type = GenomeTrack::get_type(trackpath.c_str(), iu.get_chromkey());

		if (type == GenomeTrack::FIXED_BIN) {
			enum { BINSIZE = NUM_COLS, NUM_FIXED_BIN_COLS };

			rprotect(answer = allocVector(VECSXP, NUM_FIXED_BIN_COLS));
			rprotect(names = allocVector(STRSXP, NUM_FIXED_BIN_COLS));

			GenomeTrackFixedBin gtrack;
			GIntervals all_genome_intervs;
			iu.get_all_genome_intervs(all_genome_intervs);
			char filename[FILENAME_MAX];

			sprintf(filename, "%s/%s", trackpath.c_str(), GenomeTrack::get_1d_filename(iu.get_chromkey(), all_genome_intervs.front().chromid).c_str());

			gtrack.init_read(filename, all_genome_intervs.front().chromid);
			SET_VECTOR_ELT(answer, BINSIZE, allocVector(INTSXP, 1));
			INTEGER(VECTOR_ELT(answer, BINSIZE))[0] = gtrack.get_bin_size();
			SET_STRING_ELT(names, BINSIZE, mkChar("bin.size"));
		} else {
			rprotect(answer = allocVector(VECSXP, NUM_COLS));
			rprotect(names = allocVector(STRSXP, NUM_COLS));
		}

		SET_VECTOR_ELT(answer, TYPE, allocVector(STRSXP, 1));
		SET_STRING_ELT(VECTOR_ELT(answer, TYPE), 0, mkChar(GenomeTrack::TYPE_NAMES[type]));
		SET_STRING_ELT(names, TYPE, mkChar("type"));

		SET_VECTOR_ELT(answer, DIM, allocVector(INTSXP, 1));
		INTEGER(VECTOR_ELT(answer, DIM))[0] = GenomeTrack::is_1d(type) ? 1 : 2;
		SET_STRING_ELT(names, DIM, mkChar("dimensions"));

		int64_t totsize = 0;
		vector<string> filenames;
		struct stat buf;

		get_chrom_files(trackpath.c_str(), filenames);
		for (vector<string>::const_iterator ifilename = filenames.begin(); ifilename != filenames.end(); ++ifilename) {
			string fullpath(trackpath + "/" + *ifilename);
			if (stat(fullpath.c_str(), &buf))
				verror("Cannot stat %s: %s", fullpath.c_str(), strerror(errno));
			totsize += buf.st_size;
		}
		SET_VECTOR_ELT(answer, SIZE, allocVector(REALSXP, 1));
		REAL(VECTOR_ELT(answer, SIZE))[0] = totsize;
		SET_STRING_ELT(names, SIZE, mkChar("size.in.bytes"));

		setAttrib(answer, R_NamesSymbol, names);

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	}
	return R_NilValue;
}

}
