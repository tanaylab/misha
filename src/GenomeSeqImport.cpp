#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include <set>
#include <string>
#include <vector>

#include "rdbutils.h"
#include "BufferedFile.h"

using namespace std;
using namespace rdb;

// The contig name the indexed import gives a FASTA header (GenomeSeqMultiImport.cpp)
string sanitize_fasta_header(const string &header);

namespace {

class SeqWriter {
public:
	void open(const string &path) {
		if (m_file.open(path.c_str(), "w"))
			verror("Failed to open file %s: %s", path.c_str(), strerror(errno));
	}

	void put(char c) {
		m_buf.push_back(c);
		if (m_buf.size() >= 1048576)
			flush_buf();
	}

	void close() {
		flush_buf();
		// a full disk shows here, not at write() (see BufferedFile::flush)
		if (m_file.flush() != 0 || m_file.close() != 0)
			verror("Error while writing %s: %s", m_file.file_name().c_str(), strerror(errno));
	}

private:
	void flush_buf() {
		if (m_buf.empty())
			return;
		m_file.write(&m_buf.front(), m_buf.size());
		if (m_file.error())
			verror("Error while writing %s: %s", m_file.file_name().c_str(), strerror(errno));
		m_buf.clear();
	}

	BufferedFile m_file;
	vector<char> m_buf;
};

bool file_exists(const string &path) {
	struct stat st;
	return stat(path.c_str(), &st) == 0;
}

void rename_file(const string &from, const string &to) {
	if (rename(from.c_str(), to.c_str()))
		verror("Failed to rename %s to %s: %s", from.c_str(), to.c_str(), strerror(errno));
}

}

extern "C" {

// Imports a FASTA file into the per-chromosome .seq files of _seqdir. A file of one record (or of
// sequence lines without a header) gives one chromosome, named _name (the file name, as before). A
// file of several records gives one chromosome per record, named after its header as the indexed
// format names it. Returns the names of the chromosomes, in file order.
SEXP gseqimport(SEXP _fasta, SEXP _seqdir, SEXP _name, SEXP _envir)
{
	try {
		RdbInitializer rdb_init;

		if (!Rf_isString(_fasta) || Rf_length(_fasta) != 1)
			verror("Fasta argument is not a string");

		if (!Rf_isString(_seqdir) || Rf_length(_seqdir) != 1)
			verror("Seq directory argument is not a string");

		if (!Rf_isString(_name) || Rf_length(_name) != 1)
			verror("Name argument is not a string");

		const char *fasta_fname = CHAR(STRING_ELT(_fasta, 0));
		string seqdir(CHAR(STRING_ELT(_seqdir, 0)));
		string name(CHAR(STRING_ELT(_name, 0)));
		BufferedFile fasta_file;

		if (fasta_file.open(fasta_fname, "r"))
			verror("Failed to open file %s: %s", fasta_fname, strerror(errno));

		// The first record goes to a temporary file: its name is known only once the file
		// turns out to hold one record (the file name) or more (its header).
		string tmp_path = seqdir + "/." + name + ".seq.tmp";
		SeqWriter writer;
		writer.open(tmp_path);

		vector<string> names;       // the header of each record
		set<string> taken;
		bool leading_seq = false;   // sequence lines before the first header

		auto start_record = [&](const string &header) {
			string contig = sanitize_fasta_header(header);

			if (names.empty() && !leading_seq) {
				names.push_back(contig);   // the first record, in the temporary file
				return;
			}

			if (names.empty() || leading_seq)
				verror("FASTA file %s has sequence before its first header and more than one record", fasta_fname);

			writer.close();
			if (names.size() == 1) {
				// a second record: the first one is a chromosome of its own too, named after its header
				string path = seqdir + "/" + names[0] + ".seq";
				if (file_exists(path))
					verror("Duplicate contig name '%s' (FASTA file %s)", names[0].c_str(), fasta_fname);
				rename_file(tmp_path, path);
				taken.insert(names[0]);
			}

			string path = seqdir + "/" + contig + ".seq";
			if (taken.count(contig) || file_exists(path))
				verror("Duplicate contig name '%s' (FASTA file %s)", contig.c_str(), fasta_fname);
			writer.open(path);
			names.push_back(contig);
			taken.insert(contig);
		};

		string line;
		int c;

		do {
			line.clear();
			while ((c = fasta_file.getc()) != EOF && c != '\n')
				line.push_back((char)c);

			if (line.empty() || line[0] == ';')
				continue;

			if (line[0] == '>') {
				start_record(line);
				continue;
			}

			for (char ch : line) {
				if (isalpha(ch) || ch == '-') {
					if (names.empty())
						leading_seq = true;
					writer.put(ch);
				} else if (!isspace(ch))
					verror("Invalid format of FASTA file %s", fasta_fname);
			}
		} while (c != EOF);

		if (fasta_file.error())
			verror("Error while reading %s: %s", fasta_fname, strerror(errno));

		writer.close();

		vector<string> chroms;
		if (names.size() <= 1) {
			rename_file(tmp_path, seqdir + "/" + name + ".seq");
			chroms.push_back(name);
		} else
			chroms = names;

		SEXP answer = rprotect_ptr(RSaneAllocVector(STRSXP, chroms.size()));
		for (size_t i = 0; i < chroms.size(); ++i)
			SET_STRING_ELT(answer, i, RSaneMkChar(chroms[i].c_str()));
		return answer;
	} catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }
	return R_NilValue;
}

}
