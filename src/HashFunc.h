#ifndef HASHFUNC_H_
#define HASHFUNC_H_

#include <cstdint>
#include <functional>

#ifndef BSWAP_8
#define	BSWAP_8(x)	((x) & 0xff)
#endif

#ifndef BSWAP_16
#define	BSWAP_16(x)	((BSWAP_8(x) << 8) | BSWAP_8((x) >> 8))
#endif

#ifndef BSWAP_32
#define	BSWAP_32(x)	((BSWAP_16(x) << 16) | BSWAP_16((x) >> 16))
#endif

#ifndef BSWAP_64
#define	BSWAP_64(x)	((BSWAP_32(x) << 32) | BSWAP_32((x) >> 32))
#endif

namespace std
{
	template<> struct hash<std::pair<uint64_t, uint64_t>>
	{
		uint64_t operator()(const std::pair<uint64_t, uint64_t> &v) const {
			uint64_t h = v.first * 0x9e3779b97f4a7c15ULL;
			h ^= v.second * 0x517cc1b727220a95ULL;
			return h ^ (h >> 32);
		}
	};


	template<> struct hash<std::pair<uint32_t, uint32_t>>
	{
		uint64_t operator()(const std::pair<uint32_t, uint32_t> &v) const {
			return v.first ^ BSWAP_32(v.second);
		}
	};
}

#endif /* HASHFUNC_H_ */
