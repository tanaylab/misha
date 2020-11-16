#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#ifdef __linux__
#include <sys/sendfile.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>

#include "FileUtils.h"
#include "TGLException.h"

struct FD {
    int fd{-1};
    ~FD() {
        if (fd != -1)
            close(fd);
    }
};

#ifndef __linux__
struct BUFFER {
    size_t size;
    char *buffer;
    BUFFER(size_t s) : size((s > 16384)? s : 16384), buffer(new char[size]) { }
    ~BUFFER() {
        delete [] buffer;
    }
};
#endif

void FileUtils::copy_file(const char *src, const char *tgt)
{
    FD srcfd;
    FD tgtfd;
    struct stat srcstat;

    if ((srcfd.fd = open(src, O_RDONLY, 0)) == -1)
        TGLError(errno, "Error opening file %s for reading: %s", src, strerror(errno));
    if (fstat(srcfd.fd, &srcstat) == -1)
        TGLError(errno, "Error trying to stat file %s: %s", src, strerror(errno));
    if ((tgtfd.fd = creat(tgt, srcstat.st_mode)) == -1)
        TGLError(errno, "Error opeining file %s for writing: %s", tgt, strerror(errno));
#ifdef __linux__
    if (sendfile(tgtfd.fd, srcfd.fd, NULL, srcstat.st_size) == -1)
        TGLError(errno, "Error copying file %s to %s: %s\n", src, tgt, strerror(errno));
#else
    BUFFER buffer(srcstat.st_blksize);
    if (buffer.buffer == NULL)
        TGLError(errno, "Error allocating buffer: %s\n", strerror(errno));

    ssize_t nread;
    while ((nread = read(srcfd.fd, buffer.buffer, buffer.size)) > 0) {
        const char *p = buffer.buffer;
        size_t remaining = nread;
        while (remaining > 0) {
            ssize_t nwritten = write(tgtfd.fd, p, remaining);
            if (nwritten <= 0)
                TGLError(errno, "Error writing file %s while copying from %s: %s\n", tgt, src, strerror(errno));
            p += nwritten;
            remaining -= nwritten;
        }
    }
    if (nread < 0)
        TGLError(errno, "Error reading file %s while copying to %s: %s\n", src, tgt, strerror(errno));
#endif
}

void FileUtils::move_file(const char *src, const char *tgt)
{
    if (rename(src, tgt) == -1) {
        if (errno == EXDEV) {
            FileUtils::copy_file(src, tgt);
            if (unlink(src) == -1) {
                auto olderrno = errno;
                unlink(tgt);
                TGLError(olderrno, "Error removing file %s: %s", src, strerror(olderrno));
            }
        } else
            TGLError(errno, "Error moving file %s to %s: %s\n", src, tgt);
    }
}
