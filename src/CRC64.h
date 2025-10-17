/*
 * CRC64.h
 *
 * CRC64-ECMA implementation for genome index checksums
 * Polynomial: 0xC96C5795D7870F42
 */

#ifndef CRC64_H_
#define CRC64_H_

#include <cstdint>
#include <cstddef>

namespace misha {

// CRC64-ECMA polynomial
static const uint64_t CRC64_POLY = 0xC96C5795D7870F42ULL;

// Lookup table for CRC64 (generated at first use)
class CRC64 {
public:
    CRC64() {
        init_table();
    }

    // Compute CRC64 checksum for a buffer
    uint64_t compute(const unsigned char *data, size_t len) {
        uint64_t crc = 0xFFFFFFFFFFFFFFFFULL;

        for (size_t i = 0; i < len; ++i) {
            uint8_t byte = data[i];
            crc = (crc >> 8) ^ m_table[(crc ^ byte) & 0xFF];
        }

        return ~crc;
    }

    // Compute CRC64 incrementally (for large datasets)
    uint64_t compute_incremental(uint64_t crc, const unsigned char *data, size_t len) {
        for (size_t i = 0; i < len; ++i) {
            uint8_t byte = data[i];
            crc = (crc >> 8) ^ m_table[(crc ^ byte) & 0xFF];
        }
        return crc;
    }

    // Initialize incremental computation
    uint64_t init_incremental() {
        return 0xFFFFFFFFFFFFFFFFULL;
    }

    // Finalize incremental computation
    uint64_t finalize_incremental(uint64_t crc) {
        return ~crc;
    }

private:
    uint64_t m_table[256];

    void init_table() {
        for (int i = 0; i < 256; ++i) {
            uint64_t crc = i;
            for (int j = 0; j < 8; ++j) {
                if (crc & 1) {
                    crc = (crc >> 1) ^ CRC64_POLY;
                } else {
                    crc = crc >> 1;
                }
            }
            m_table[i] = crc;
        }
    }
};

// Global instance
static CRC64 g_crc64;

// Convenience function
inline uint64_t compute_crc64(const unsigned char *data, size_t len) {
    return g_crc64.compute(data, len);
}

} // namespace misha

#endif /* CRC64_H_ */
