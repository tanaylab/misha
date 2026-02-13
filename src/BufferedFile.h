#ifndef BUFFEREDFILE_H_
#define BUFFEREDFILE_H_

#include <cstdint>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

class BufferedFile {
public:
	BufferedFile() { init(131072); }  // 128KB default buffer — reduced from 2MB to cut per-kid memory (390 tracks x 128KB = 49MB vs 780MB)
	BufferedFile(unsigned bufsize) { init(bufsize); }
	~BufferedFile();

	// returns 0 on success, -1 on failure.
    // if lock is true, fcntl lock is acquired according to the mode (see: fcntl, F_SETLKW);
    // file is unlocked on close() or destructor.
	int open(const char *path, const char *mode, bool lock = false);

	// see fclose for return value
	int close();

	// see fgetc for return value
	int getc();

	// see fread for return value
	uint64_t read(void *ptr, uint64_t size);

	// see fwrite for return value
	uint64_t write(const void *ptr, uint64_t size);

	// see ftell for return value
	int64_t tell() const { return m_virt_pos; }

	// see fseek for return value
	int seek(int64_t offset, int whence);

	int error() const { return !m_fp || ferror(m_fp); }

	int eof() const { return m_eof; }

	bool opened() const { return m_fp != NULL; }

	const std::string &file_name() const { return m_filename; }

	int64_t file_size() const { return m_file_size; }

    // truncates the file beyond the current position
    int truncate();

    int stat(struct stat *fs) { return fstat(fileno(m_fp), fs); }

	// in case of error throws an exception
	static int64_t file_size(const char *path);

protected:
	FILE       *m_fp{NULL};
	int         m_eof;
	std::string m_filename;
	char       *m_buf{NULL};
	unsigned    m_bufsize;
	int64_t     m_file_size;
	int64_t     m_virt_pos;
	int64_t     m_phys_pos;
	int64_t     m_sbuf_pos;
	int64_t     m_ebuf_pos;

	void init(unsigned bufsize);

private:
	// restrict copying of class (otherwise need reference counters for m_buf and m_fp)
	BufferedFile(const BufferedFile &) = delete;
	BufferedFile &operator=(const BufferedFile &obj) = delete;
};


class BufferedFiles : public std::vector<BufferedFile *> {
public:
	BufferedFiles() : std::vector<BufferedFile *>() {}
	BufferedFiles(size_type n) : std::vector<BufferedFile *>(n, NULL) {}

	~BufferedFiles() {
		for (iterator ifile = begin(); ifile != end(); ++ifile)
			delete *ifile;
	}

private:
	BufferedFiles(const BufferedFiles &) = delete;
	BufferedFiles &operator=(const BufferedFiles &) = delete;
};

//--------------------------- IMPLEMENTATION ------------------------------------

inline BufferedFile::~BufferedFile()
{
	close();
	delete []m_buf;
}

inline void BufferedFile::init(unsigned bufsize) {
	m_fp = NULL;
	m_eof = true;
	m_buf = NULL;
	m_bufsize = bufsize;
	m_file_size = 0;
	m_virt_pos = -1;
	m_phys_pos = 0;
	m_sbuf_pos = 0;
	m_ebuf_pos = 0;
	m_buf = new char[m_bufsize];
}

inline int BufferedFile::getc()
{
	// is the new read already cached?
	if (m_virt_pos >= m_sbuf_pos && m_virt_pos + 1 <= m_ebuf_pos)
		return m_buf[m_virt_pos++ - m_sbuf_pos];

	char c;
	if (!read(&c, 1))
		return -1;
	return c;
}

inline uint64_t BufferedFile::read(void *ptr, uint64_t size)
{
    // is the new read already cached?
    if (m_virt_pos >= m_sbuf_pos && m_virt_pos + (long)size <= m_ebuf_pos) {
        memcpy(ptr, m_buf + m_virt_pos - m_sbuf_pos, size);
        m_virt_pos += size;
        return size;
    }

    // we must perform a new read
    if (m_phys_pos != m_virt_pos)
        fseeko(m_fp, m_virt_pos, SEEK_SET);

    // For very large requests (>= buffer size), bypass cache and read directly
    if (size >= m_bufsize) {
#ifdef POSIX_FADV_SEQUENTIAL
        int fd = fileno(m_fp);
        posix_fadvise(fd, m_virt_pos, (off_t)size, POSIX_FADV_SEQUENTIAL);
#endif
        uint64_t bytes_read = fread(ptr, 1, size, m_fp);
        m_virt_pos += bytes_read;
        m_phys_pos = m_virt_pos;
        // Invalidate buffer window after large direct read
        m_sbuf_pos = m_ebuf_pos = 0;
        if (!bytes_read && feof(m_fp))
            m_eof = true;
        return bytes_read;
    }

    // If request is smaller than buffer size, perform cached read with modest prefetch
    if (size < m_bufsize) {
        uint64_t to_read = size;
        const uint64_t modest_prefetch_small = 128 * 1024;  // small reads benefit
        if (size < 100 * 1024 && to_read < modest_prefetch_small)
            to_read = modest_prefetch_small;

        if (to_read > m_bufsize)
            to_read = m_bufsize;

        // Advise kernel about access pattern
#ifdef POSIX_FADV_SEQUENTIAL
        int fd = fileno(m_fp);
        if (to_read >= 256 * 1024)
            posix_fadvise(fd, m_virt_pos, (off_t)to_read, POSIX_FADV_SEQUENTIAL);
        else
            posix_fadvise(fd, m_virt_pos, (off_t)to_read, POSIX_FADV_RANDOM);
#endif

        uint64_t bytes_read = fread(m_buf, 1, to_read, m_fp);

        m_phys_pos = m_virt_pos + bytes_read;
        m_sbuf_pos = m_virt_pos;
        m_ebuf_pos = m_phys_pos;

        uint64_t ret = bytes_read;
        if (ret > size)
            ret = size;

        m_virt_pos += ret;
        memcpy(ptr, m_buf, ret);
        if (!bytes_read && feof(m_fp))
            m_eof = true;
        return ret;
    }

    // Fallback (should not be reached often)
    uint64_t bytes_read = fread(ptr, 1, size, m_fp);
    m_virt_pos += bytes_read;
    m_phys_pos = m_virt_pos;
    if (!bytes_read && feof(m_fp))
        m_eof = true;
    return bytes_read;
}

inline uint64_t BufferedFile::write(const void *ptr, uint64_t size)
{
	if (m_phys_pos != m_virt_pos) {
		fseeko(m_fp, m_virt_pos, SEEK_SET);
		m_phys_pos = m_virt_pos;
	}

	uint64_t retv = fwrite(ptr, 1, size, m_fp);

	if (retv) {
		// if the cached buffer overlaps the written region => trash the buffer
		if (std::max(m_virt_pos, m_sbuf_pos) < std::min(m_virt_pos + (long)retv, m_ebuf_pos))
			m_sbuf_pos = m_ebuf_pos = 0;

		m_virt_pos += retv;
		m_phys_pos = m_virt_pos;

		if (m_file_size < m_virt_pos)
			m_file_size = m_virt_pos;
	}

	return retv;
}

inline int BufferedFile::seek(int64_t offset, int whence)
{
	switch (whence) {
	case SEEK_END:
		if (offset < 0 || m_file_size - 1 < offset) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos = m_file_size - offset - 1;
		break;
	case SEEK_CUR:
		if (m_virt_pos + offset < 0 || m_virt_pos + offset > m_file_size - 1) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos += offset;
		break;
	case SEEK_SET:
		if (offset < 0 || offset > m_file_size) {
			errno = EINVAL;
			return -1;
		}
		m_virt_pos = offset;
		break;
	default:
		errno = EINVAL;
		return -1;
	}

	m_eof = m_virt_pos == m_file_size;
	return 0;
}

#endif
