#include "Fastq.h"

void Fastq::init(const string &filename)
{
	if (m_bfile.open(filename.c_str(), "r"))
		TGLError<Fastq>(FILE_ERROR, "Failed to open FASTQ file %s: %s", m_bfile.file_name().c_str(), strerror(errno));

	m_seq_len = 0;
	m_lineno = 1;
}

bool Fastq::read_record()
{
	int record_lineno = 0;
	string *str = &m_buf;

	m_buf.clear();
	m_seq.clear();

	while (1) {
		int c = (char)m_bfile.getc();

		if (m_bfile.error())
			TGLError<Fastq>(FILE_ERROR, "Error while reading FASTQ file %s: %s", m_bfile.file_name().c_str(), strerror(errno));

		if (c == '\r')
			c = '\n';

		if ((m_bfile.eof() || c == '\n') && record_lineno == 3 && str->length()) {
			if (m_buf.length() != m_seq.length())
				TGLError<Fastq>(BAD_QUAL, "File %s, line %d: quality string length differs from sequence length", m_bfile.file_name().c_str(), m_lineno);

			m_quality.init(m_buf, get_phred_base());

			if (c == '\n')
				m_lineno++;
			return true;
		}

		if (m_bfile.eof()) {
			if (record_lineno)
				TGLError<Fastq>(INVALID_FORMAT, "Invalid format of FASTQ file %s", m_bfile.file_name().c_str());
			break;
		}

		if (c == '\n') {
			m_lineno++;
			if (!str->empty()) {
				record_lineno++;

				if (record_lineno == 1) {
					size_t pos = str->rfind(':');

					if (pos == string::npos || pos == str->length() - 1)
						TGLError<Fastq>(INVALID_FORMAT, "Invalid format of FASTQ file %s", m_bfile.file_name().c_str());

					m_id.assign(*str, 0, pos);
					m_index_seq.assign(*str, pos, string::npos);
					str = &m_seq;
				} else
					str = &m_buf;
				str->clear();
			}
		} else
			str->push_back(str == &m_seq ? (char)toupper(c) : c);
	}
	return false;
}
