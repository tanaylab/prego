#ifndef FASTQ_H_
#define FASTQ_H_

#include "BufferedFile.h"
#include "SeqQual.h"
#include "TGLException.h"

//----------------- IN CASE OF ERROR THIS CLASS THROWS TGLException ------------------

class Fastq {
public:
	enum Errors { FILE_ERROR, INVALID_FORMAT, BAD_QUAL };

	Fastq(const string &filename) { init(filename); }
	Fastq() {}

	void init(const string &filename);

	// returns true if record was read, false if eof reached
	bool read_record();

	const string  &last_seq() const { return m_seq; }
	const SeqQual &last_quality() const { return m_quality; }
	const string  &last_index_seq() const { return m_index_seq; }
	const string  &last_id() const { return m_id; }
	unsigned char  get_phred_base() const { return 33; }

private:
	BufferedFile m_bfile;
	string       m_buf;
	string       m_seq;
	SeqQual      m_quality;
	string       m_index_seq;
	string       m_id;
	int          m_seq_len;
	int          m_lineno;
};

#endif /* FASTQ_H_ */
