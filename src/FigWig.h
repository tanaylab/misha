/*
 * FigWig.h
 *
 * A reader of figwig files, version 1 (https://github.com/jmschrei/figwig, docs/format.rst):
 * the tables, the metadata and the values of position-order chunks (dense, sparse, matrix and
 * runs, in every element type) and of track-order blocks, codec 0 (stored) and 1 (zlib). The
 * summaries, zoom levels and the regional summary index are not needed for the values and are
 * not read.
 */

#ifndef FIGWIG_H_
#define FIGWIG_H_

#include <cstdint>
#include <string>
#include <vector>

class FigWig {
public:
	struct Attr {
		std::string key;
		std::string text;   // a string as it is, anything else as its JSON text
		bool        null;   // the value is JSON null
	};

	struct Track {
		std::string       name;
		std::string       dtype;
		int               group;
		int               slot;    // the track's place among its group's tracks
		std::vector<Attr> attrs;   // the keys of the record other than name, dtype and group
	};

	struct Chrom {
		std::string name;
		uint64_t    size;
		uint64_t    first_tile;
	};

	FigWig() {}
	~FigWig();
	FigWig(const FigWig &) = delete;
	FigWig &operator=(const FigWig &) = delete;

	void open(const char *path);

	const std::vector<Chrom> &chroms() const { return m_chroms; }
	const std::vector<Track> &tracks() const { return m_tracks; }
	bool                      track_order() const { return m_track_order; }
	uint32_t                  tile() const { return m_tile; }
	bool                      has_genome() const { return m_has_genome; }
	const std::string        &genome() const { return m_genome; }
	int                       group_size(int group) const { return m_group_tracks[group].size(); }

	// Position order: the values of tile 'tile' (numbered over the genome) of group 'group', for
	// the tracks of the group in 'slots', into out[i * tile() + offset] for slots[i]. Returns false,
	// leaving out alone, when the tile holds no chunk: every value of the group is 0 there.
	bool read_tile(int group, uint64_t tile, const std::vector<int> &slots, double *out);

	// Track order: the nonzero bases of a track on a chromosome, a block at a time, in ascending
	// position. Returns the number of blocks; read_block() decodes one of them.
	uint64_t num_blocks(int track, int chromid);
	void     read_block(int track, int chromid, uint64_t iblock, std::vector<uint32_t> &positions, std::vector<double> &values);

private:
	struct LiveIndex {   // position order: the live tiles of a group
		bool                  loaded{false};
		std::vector<uint32_t> tiles;
		std::vector<uint8_t>  kinds;
		std::vector<uint32_t> sizes;
		std::vector<uint32_t> raws;
		std::vector<uint32_t> crcs;
		std::vector<uint64_t> offsets;
	};

	struct Block {
		uint32_t chrom;
		uint32_t first;
		uint32_t last;
		uint32_t n;
		uint64_t offset;
		uint32_t size;
		uint32_t raw;
		uint32_t kind;
		uint32_t crc;
	};

	struct BlockIndex {  // track order: the blocks of a track, and where each chromosome's start
		bool                  loaded{false};
		std::vector<Block>    blocks;
		std::vector<uint64_t> chrom_start;  // n_chroms + 1 entries
	};

	std::string m_path;
	int         m_fd{-1};
	uint64_t    m_length{0};
	uint64_t    m_chroms_offset{0};
	uint64_t    m_index_offset{0};
	uint32_t    m_tile{0};
	uint16_t    m_codec{0};
	bool        m_track_order{false};
	uint32_t    m_num_groups{0};
	uint64_t    m_num_tiles{0};
	bool        m_has_genome{false};
	std::string m_genome;

	std::vector<Chrom>            m_chroms;
	std::vector<Track>            m_tracks;
	std::vector<std::vector<int>> m_group_tracks;  // track indices of each group, in slot order
	std::vector<uint64_t>         m_directory;
	std::vector<LiveIndex>        m_live;
	std::vector<BlockIndex>       m_blocks;

	std::vector<unsigned char> m_filebuf;    // a chunk or block as it lies in the file
	std::vector<unsigned char> m_decoded;    // its decoded bytes

	void read_at(uint64_t offset, void *buf, uint64_t size, const char *what);
	void parse_metadata(const std::string &json);
	void load_live(int group);
	void load_blocks(int track);
	// the decoded bytes of a chunk or block into m_decoded
	void decode(uint64_t offset, uint32_t size, uint32_t raw, uint32_t crc, bool stored, const char *what);
	[[noreturn]] void corrupt(const char *fmt, ...) const;
};

#endif /* FIGWIG_H_ */
