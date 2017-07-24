#ifndef _XOPEN_SOURCE
	#define _XOPEN_SOURCE 500
#endif
#include <ftw.h>

#include "rdbutils.h"

using namespace std;
using namespace rdb;

static vector<string> s_tracks;
static vector<string> s_intervs;

static int pick_tracks_n_intervals(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf) {
	const char *fname = fpath + ftwbuf->base;
	unsigned fname_len = strlen(fname);

	if (tflag == FTW_D) {
		if (fname_len > rdb::TRACK_FILE_EXT.size() && !strcmp(fname + fname_len - rdb::TRACK_FILE_EXT.size(), rdb::TRACK_FILE_EXT.c_str())) {
			s_tracks.push_back(fpath);
			return FTW_SKIP_SUBTREE;
		}

		if (fname_len > rdb::INTERV_FILE_EXT.size() && !strcmp(fname + fname_len - rdb::INTERV_FILE_EXT.size(), rdb::INTERV_FILE_EXT.c_str())) {
			s_intervs.push_back(fpath);
			return FTW_SKIP_SUBTREE;
		}

		// skip directories that have "." in their names (if not skipped these names will make a mess while trying to resolve a track/interv name)
		for (const char *p = fpath + ftwbuf->base; *p; p++) {
			if (*p == '.')
				return FTW_SKIP_SUBTREE;
		}
	}

	if (tflag == FTW_F && fname_len > rdb::INTERV_FILE_EXT.size() && !strcmp(fname + fname_len - rdb::INTERV_FILE_EXT.size(), rdb::INTERV_FILE_EXT.c_str()))
		s_intervs.push_back(fpath);

	return FTW_CONTINUE;
}

extern "C" {

SEXP gfind_tracks_n_intervals(SEXP _dir, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		// check the arguments
		if (!isString(_dir) || length(_dir) != 1)
			verror("Dir argument argument is not a string");

		const char *dir = CHAR(STRING_ELT(_dir, 0));
		IntervUtils iu(_envir);

        s_tracks.clear();
        s_intervs.clear();
		nftw(dir, pick_tracks_n_intervals, 50, FTW_ACTIONRETVAL);

		SEXP answer;
		SEXP tracks;
		SEXP intervs;

		rprotect(tracks = allocVector(STRSXP, s_tracks.size()));
		for (vector<string>::iterator itrack = s_tracks.begin(); itrack < s_tracks.end(); ++itrack) {
			itrack->resize(itrack->size() - rdb::TRACK_FILE_EXT.size());
			for (size_t pos = itrack->find_first_of('/'); pos != string::npos; pos = itrack->find_first_of('/', pos + 1))
				itrack->at(pos) = '.';
			SET_STRING_ELT(tracks, itrack - s_tracks.begin(), mkChar(itrack->c_str() + strlen(dir) + 1));
		}

		rprotect(intervs = allocVector(STRSXP, s_intervs.size()));
		for (vector<string>::iterator iinterv = s_intervs.begin(); iinterv < s_intervs.end(); ++iinterv) {
			iinterv->resize(iinterv->size() - rdb::INTERV_FILE_EXT.size());
			for (size_t pos = iinterv->find_first_of('/'); pos != string::npos; pos = iinterv->find_first_of('/', pos + 1))
				iinterv->at(pos) = '.';
			SET_STRING_ELT(intervs, iinterv - s_intervs.begin(), mkChar(iinterv->c_str() + strlen(dir) + 1));
		}

		rprotect(answer = allocVector(VECSXP, 2));
		SET_VECTOR_ELT(answer, 0, tracks);
		SET_VECTOR_ELT(answer, 1, intervs);

		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
	}
	return R_NilValue;
}

}
