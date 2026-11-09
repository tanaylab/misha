#ifndef GENOMEUTILS_H_
#define GENOMEUTILS_H_

#include <string>

extern char s_complementary_basepair[256];
static struct Complementary_basepair_initializer {
	Complementary_basepair_initializer();
	static bool s_initialized;
} __complementary_basepair_initializer;

// unsigned char: a plain char is signed here, and a byte >= 128 would index the table below 0.
inline char basepair2complementary(char c) { return s_complementary_basepair[(unsigned char)c]; }

std::string seq2reverse_complementary(const std::string &seq);
std::string seq2reverse(const std::string &seq);
std::string seq2complementary(const std::string &seq);

#endif
