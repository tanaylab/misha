/*
 * fsync a staged directory tree, or a single directory entry, from R.
 *
 * The staged writers (R/track-create-atomic.R, GenomeTrackModify's Stager)
 * commit with rename(), which is atomic with respect to other *processes*
 * but says nothing about stable storage: after a power loss or a kernel
 * panic the committed name can point at partially written data. Closing
 * that hole needs an fsync of every staged file, of the staging directory,
 * and - after the rename - of the parent directory that now holds the new
 * name. R has no fsync, hence this .Call.
 *
 * Cost, measured 2026-08-24 on a 1.28 GB dense track (77 s to create):
 * nil on the lab NFS (the client already flushes at close()), nil on
 * tmpfs, ~0.35 s per GB on local ext4. Under 1% either way, so this is
 * unconditional rather than an option.
 *
 * Exposed as .gdb.fsync(path, recursive).
 */

#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include "rdbutils.h"

using namespace std;
using namespace rdb;

// fsync one already-existing path. A directory fsync is best-effort: some
// filesystems answer EINVAL, which means "not supported", not "lost data".
static void fsync_one(const string &path, bool is_dir)
{
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        if (is_dir && errno == EACCES)   // unreadable dir: nothing we can do
            return;
        verror("Failed to open %s for fsync: %s", path.c_str(), strerror(errno));
    }

    int rc = fsync(fd);
    int err = errno;
    close(fd);

    if (rc != 0 && !(is_dir && err == EINVAL))
        verror("Failed to fsync %s: %s", path.c_str(), strerror(err));
}

// Depth-first: files before the directory that names them, children before
// the parent, so a surviving directory entry always points at synced data.
static void fsync_tree(const string &dir)
{
    vector<string> files, subdirs;

    DIR *d = opendir(dir.c_str());
    if (!d)
        verror("Failed to open directory %s: %s", dir.c_str(), strerror(errno));

    // Collect first, act after closedir: verror() unwinds and would leak the
    // handle otherwise.
    for (struct dirent *e = readdir(d); e; e = readdir(d)) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, ".."))
            continue;

        string path = dir + "/" + e->d_name;
        struct stat st;
        if (lstat(path.c_str(), &st))
            continue;                    // vanished under us; nothing to sync

        if (S_ISDIR(st.st_mode))
            subdirs.push_back(path);
        else if (S_ISREG(st.st_mode))
            files.push_back(path);
    }
    closedir(d);

    for (vector<string>::const_iterator i = files.begin(); i != files.end(); ++i)
        fsync_one(*i, false);
    for (vector<string>::const_iterator i = subdirs.begin(); i != subdirs.end(); ++i)
        fsync_tree(*i);

    fsync_one(dir, true);
}

extern "C" {

SEXP gdb_fsync(SEXP _path, SEXP _recursive, SEXP _envir)
{
    try {
        if (!Rf_isString(_path) || Rf_length(_path) < 1)
            verror("'path' must be a character vector");
        if (!Rf_isLogical(_recursive) || Rf_length(_recursive) != 1)
            verror("'recursive' must be a single logical");

        const bool recursive = LOGICAL(_recursive)[0] == TRUE;
        const int n = Rf_length(_path);

        for (int i = 0; i < n; ++i) {
            const char *path = CHAR(STRING_ELT(_path, i));
            if (path == NULL || path[0] == '\0')
                continue;

            struct stat st;
            if (stat(path, &st))
                verror("Failed to stat %s: %s", path, strerror(errno));

            if (recursive && S_ISDIR(st.st_mode))
                fsync_tree(path);
            else
                fsync_one(path, S_ISDIR(st.st_mode));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    }
    return R_NilValue;
}

}
