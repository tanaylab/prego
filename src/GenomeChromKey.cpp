#include "port.h"
BASE_CC_FILE

#include <errno.h>
#include <stdlib.h>
#include <fstream>

#include "strutil.h"

#include "GenomeChromKey.h"

void GenomeChromKey::read_chroms_sizes_file(const char *fname)
{
	vector<string> fields;
	ifstream fs(fname);
	vector<Chrom> chroms;

	m_name2id.clear();
	m_id2chrom.clear();
	m_id = 0;

	while (1) {
		split_line(fs, fields);

		if (fs.eof())
			break;

		if (fs.fail())
			TGLError<GenomeChromKey>(FILE_READ_FAILED, "Failed to read chrom sizes file %s: %s", fname, strerror(errno));

		if (fields.size() != 2)
			TGLError<GenomeChromKey>(FILE_READ_FAILED, "Invalid format of chrom sizes file %s", fname);

		char *endptr;
		int64_t chrom_len = strtoll(fields[1].c_str(), &endptr, 10);
		if (*endptr)
			TGLError<GenomeChromKey>(BAD_FILE_FORMAT, "Reading chrom sizes file %s: size of chromosome %s is not an integer", fname, fields[0].c_str());

		if (chrom_len < 0)
			TGLError<GenomeChromKey>(BAD_FILE_FORMAT, "Reading chrom sizes file %s: size of chromosome %s is out of range", fname, fields[0].c_str());

		chroms.push_back(Chrom((string("chr") + fields[0]).c_str(), chrom_len));
	}
	fs.close();

	sort(chroms.begin(), chroms.end());
	for (vector<Chrom>::const_iterator i = chroms.begin(); i != chroms.end(); ++i)
		add_chrom(i->name, i->size);
}
