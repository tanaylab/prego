#include "port.h"
BASE_CC_FILE
#include "GenomeSequence.h"
#include "FastaParse.h"

void GenomeSequence::read_key_fasta(istream &key_tab)
{
	string chrom, fn;
	key_tab >> chrom;
	while(key_tab) {
		key_tab >> fn;
		ASSERT(m_ChromId.find(chrom) == m_ChromId.end(), "Double chromosome " << chrom << " on key");
		int cid  = m_Seq.size();
		m_ChromId[chrom] = cid;
		m_Seq.resize(cid + 1);
		m_ChromName.push_back(chrom);
		m_ChromFn.resize(cid + 1);
		m_ChromFn[cid] = fn;
	
//currently we read all chromos, this should be changed to something more memory aware
		read_chrom_fasta(cid);
	
		key_tab >> chrom;
	}
}

void GenomeSequence::read_chrom_fasta(int cid)
{
	FastaParse seq(m_ChromFn[cid], &(m_Seq[cid]), m_capitalize);
}
