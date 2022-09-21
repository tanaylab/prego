#include "port.h"
BASE_CC_FILE
#include "FastaParse.h"
#include "strutil.h"

FastaParse::~FastaParse()
{
	if(!m_ext_storage) {
		delete m_seq;
	}
}
void FastaParse::first()
{
	string slack;
	char code;
	m_file >> code;
	while(m_file && code == '#') {
		getline(m_file, slack, '\n');//eat comment
		m_file >> code;
	}
	if(code != '>') {
		Rcpp::Rcerr << "no header line in fasta parsig, code was " << code << endl;
		m_is_valid = false;
	} else {
		next();
		m_is_valid = true;
	}
}
bool FastaParse::next()
{
	if(!m_file) { 
		m_is_valid = false;
		return(false); 
	}
	split_line(m_file,m_aux_fields, ' ');
	if(m_aux_fields.size() == 0) {
		m_is_valid = false;
		return(false);
	}
	m_desc = m_aux_fields[0];
	*m_seq = "";

	Rcpp::Rcerr << "got desc " << m_desc << endl;

	//read linse until ">"

	char code;
	string line;
	m_file >> code;
	while(m_file && code != '>') {
		*m_seq += code;
		getline(m_file, line, '\n');
		*m_seq += line;
		m_file >> code;
	}
	if(m_capitalize) {
		for(string::iterator i = m_seq->begin(); i != m_seq->end(); i++) {
			if(*i > 'Z') *i -= 'a'-'A';
		}
	}
	return(true);
}

