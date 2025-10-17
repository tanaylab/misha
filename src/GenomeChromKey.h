#ifndef GENOMECHROMKEY_H_
#define GENOMECHROMKEY_H_

#include <cstdint>
#include <stdint.h>
#include <unordered_map>
#include <vector>
#include "HashFunc.h"
#include "TGLException.h"

using namespace std;

// -------------------- GenomeChromKey  -----------------------
// GenomeChromKey manages the mapping between chromosome id and chromosome name.
// !!!!!!!!! IN CASE OF ERROR THIS CLASS THROWS TGLException  !!!!!!!!!!!!!!!!

class GenomeChromKey {
public:
	enum Errors { CHROM_EXISTS, CHROM_NOEXISTS, ID_NOEXISTS, FILE_READ_FAILED, BAD_FILE_FORMAT, NUM_ERRORS };

	GenomeChromKey() : m_id(0) {}

	int           chrom2id(const string &chrom) const;
	int           chrom2id(const char *chrom) const;
	const string &id2chrom(int id) const;
	uint64_t      get_chrom_size(int id) const;
	uint64_t        get_num_chroms() const { return m_id2chrom.size(); }

	void read_chroms_sizes_file(const char *fname);

	// returns id of the new chromosome
	int add_chrom(const string &chrom, uint64_t size);
	void add_chrom_alias(const string &alias, int id);
    void get_aliases(int id, vector<string> &aliases) const;

private:
	struct Chrom {
		string   name;
		uint64_t size;

		Chrom(const string &_name, uint64_t _size) : name(_name), size(_size) {}
		bool operator<(const Chrom &chrom) const { return name < chrom.name; }
	};

	typedef std::unordered_map<string, int> Name2id;
	typedef std::vector<Chrom> Id2chrom;

	Name2id    m_name2id;
	Name2id    m_alias2id;
	Id2chrom   m_id2chrom;
	int        m_id;
};


// --------------------- implementation -----------------------

inline int GenomeChromKey::add_chrom(const string &name, uint64_t size)
{
	if (m_name2id.find(name) != m_name2id.end())
		TGLError<GenomeChromKey>(CHROM_EXISTS, "Chromosome %s already exists", name.c_str());
	m_name2id[name] = m_id;
	m_id2chrom.push_back(Chrom(name, size));
	return m_id++;
}

inline int GenomeChromKey::chrom2id(const string &name) const
{
	Name2id::const_iterator iname2id = m_name2id.find(name);
	if (iname2id != m_name2id.end())
		return iname2id->second;
	Name2id::const_iterator ialias2id = m_alias2id.find(name);
	if (ialias2id != m_alias2id.end())
		return ialias2id->second;
	TGLError<GenomeChromKey>(CHROM_NOEXISTS, "Chromosome \"%s\" does not exist", name.c_str());
	return 0;
}

inline int GenomeChromKey::chrom2id(const char *name) const
{
	Name2id::const_iterator iname2id = m_name2id.find(name);
	if (iname2id != m_name2id.end())
		return iname2id->second;
	Name2id::const_iterator ialias2id = m_alias2id.find(name);
	if (ialias2id != m_alias2id.end())
		return ialias2id->second;
	TGLError<GenomeChromKey>(CHROM_NOEXISTS, "Chromosome \"%s\" does not exist", name);
	return 0;
}

inline const string &GenomeChromKey::id2chrom(int id) const
{
	if (id >= (int)m_id2chrom.size()){
		TGLError<GenomeChromKey>(ID_NOEXISTS, "Id %d cannot be mapped to any chromosome", id);
	}
	return m_id2chrom[id].name;
}

inline uint64_t GenomeChromKey::get_chrom_size(int id) const
{
	if (id >= (int)m_id2chrom.size())
		TGLError<GenomeChromKey>(ID_NOEXISTS, "Id %d cannot be mapped to any chromosome", id);
	return m_id2chrom[id].size;
}

inline void GenomeChromKey::add_chrom_alias(const string &alias, int id)
{
	if (alias.empty())
		return;
	if (m_name2id.find(alias) != m_name2id.end())
		return;
	Name2id::const_iterator ialias = m_alias2id.find(alias);
	if (ialias != m_alias2id.end()) {
		if (ialias->second == id)
			return;
		return;
	}
	m_alias2id[alias] = id;
}

inline void GenomeChromKey::get_aliases(int id, vector<string> &aliases) const
{
	aliases.clear();
	for (Name2id::const_iterator it = m_alias2id.begin(); it != m_alias2id.end(); ++it) {
		if (it->second == id)
			aliases.push_back(it->first);
	}
}

#endif /* GENOMECHROMKEY_H_ */
