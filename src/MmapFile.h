#ifndef MMAP_FILE_H
#define MMAP_FILE_H

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdint>
#include <cstddef>
#include <string>

class MmapFile {
public:
    MmapFile() : m_data(nullptr), m_size(0), m_fd(-1) {}
    ~MmapFile() { close(); }

    MmapFile(const MmapFile&) = delete;
    MmapFile& operator=(const MmapFile&) = delete;

    MmapFile(MmapFile&& other) noexcept
        : m_data(other.m_data), m_size(other.m_size), m_fd(other.m_fd) {
        other.m_data = nullptr;
        other.m_size = 0;
        other.m_fd = -1;
    }

    MmapFile& operator=(MmapFile&& other) noexcept {
        if (this != &other) {
            close();
            m_data = other.m_data;
            m_size = other.m_size;
            m_fd = other.m_fd;
            other.m_data = nullptr;
            other.m_size = 0;
            other.m_fd = -1;
        }
        return *this;
    }

    bool open(const std::string &path, bool sequential = false) {
        m_fd = ::open(path.c_str(), O_RDONLY);
        if (m_fd < 0) return false;
        struct stat st;
        if (fstat(m_fd, &st) < 0) { ::close(m_fd); m_fd = -1; return false; }
        m_size = st.st_size;
        if (m_size == 0) { ::close(m_fd); m_fd = -1; return true; }
        int flags = MAP_PRIVATE;
        // No MAP_POPULATE: with many tracks per call (e.g. ~50 motif tracks
        // × ~38MB per chromosome), eager page-in adds seconds of page-table
        // walking per chromosome transition even with a warm cache.
        m_data = static_cast<uint8_t*>(mmap(nullptr, m_size, PROT_READ, flags, m_fd, 0));
        if (m_data == MAP_FAILED) { m_data = nullptr; ::close(m_fd); m_fd = -1; return false; }
        ::close(m_fd); m_fd = -1;  // fd not needed after mmap (naryn pattern)
        // MADV_SEQUENTIAL on NFS triggers a large synchronous read-ahead
        // window per page fault — devastating for tile-clustered access
        // patterns (e.g. gextract over hundreds of motif tracks with a
        // tiled iterator), where we touch only a few KB out of each
        // ~50 MB chrom file. Measured cold-NFS: SEQUENTIAL=159s,
        // NORMAL=118s, WILLNEED=117s, RANDOM=4.7s for the same workload.
        // RANDOM is also no slower than SEQUENTIAL on warm cache.
        madvise(m_data, m_size, sequential ? MADV_SEQUENTIAL : MADV_RANDOM);
        return true;
    }

    void close() {
        if (m_data) { munmap(m_data, m_size); m_data = nullptr; }
        if (m_fd >= 0) { ::close(m_fd); m_fd = -1; }
        m_size = 0;
    }

    bool is_open() const { return m_data != nullptr; }
    const uint8_t* data() const { return m_data; }
    size_t size() const { return m_size; }

    template<typename T>
    const T& read_at(size_t offset) const {
        return *reinterpret_cast<const T*>(m_data + offset);
    }

    const uint8_t* range(size_t offset, size_t len) const {
        return m_data + offset;
    }

private:
    uint8_t* m_data;
    size_t m_size;
    int m_fd;
};

#endif // MMAP_FILE_H
