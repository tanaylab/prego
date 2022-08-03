#ifndef biodata_GenomesAlignment_h
#define biodata_GenomesAlignment_h 1

//#define REV_DBG 1

#include <vector>
#include <string>
#include <map>

#include "PhyloTree.h"
#include "MafParse.h"


class GenomesAlignment {

protected:

	const PhyloTree &m_phylo;

	vector<string> m_chroms;

	int m_base_aln_coor;

public:

	int get_base_aln_coor(){
		return(m_base_aln_coor);
	}

	const vector<string> &get_chroms() {
		return(m_chroms);
	}


//coordinates are relative to first species genome (sp=0)
	virtual string::const_iterator get_seq(int sp, const string &chrom, int fr) = 0;

	GenomesAlignment(const PhyloTree &phy) :
		m_phylo(phy)
	{
		m_base_aln_coor=0;
	}

	/*int  get_base_aln_coor(){
		return m_base_aln_coor;
	}*/

	virtual int get_max_pos(const string &crhom) = 0;
	virtual int get_refgenome_pos(int pos) const = 0;
	virtual int get_alnpos_by_refgenome(int refpos) = 0;

	virtual void set_focus_interval(const string &chr, int fr, int to) = 0;
};

class GAlignsMafsCache : public GenomesAlignment {

protected:

	map<string, string> m_chrom_maf_fn;

	string m_cur_chrom;

	vector<string> m_seq;

	int m_max_chr_coord;

	int m_max_missing_align;

	int m_spat_bin;

//This is indexing alignment coordinates to ref genome coordinates, in jumps of m_spat_bin
	vector<int> m_ref_genome_bins;
	vector<int> m_alnpos_to_ref_genome_bins;

	struct AlnFrag {
		string chrom;
		int first_coord;
		int strand;
		AlnFrag(const string &c, int fc, int s) :
			chrom(c),
			first_coord(fc),
			strand(s)
			{}
	};

//This is containing the first and last coordinate of each fragment, in alignment coorindates
	vector<int> m_frag_start;
	vector<int> m_frag_end;

//This is indexing alignment coordinate to fragment ids, in jumps of m_spat_bin
	vector<int> m_alnpos_frag;

//This is defining the mapping between fragment and the various genomes coordinates
	vector<vector<AlnFrag> > m_frag_to_genome;

	int m_base_ref_coord;
	int m_max_ref_coord;

	int m_base_global_aln_coord;
	int m_curr_global_aln_coord;
	int m_max_global_aln_coord;
	int m_use_ref_genome_coords;

public:

	int step_align_coord(int fr_aln_coord, int fr_chr_coord, int to_chr_coord) const;
	GAlignsMafsCache(const PhyloTree &phy, ifstream &key_tab, int force_full, int spat_bin = 1000);

//in alignment coordinates!
	virtual string::const_iterator get_seq(int sp, const string &chrom, int fr);

	virtual int get_refgenome_pos(int pos) const;
	virtual int get_alnpos_by_refgenome(int refpos);





	int get_max_pos(const string &chrom) {
		if(chrom != m_cur_chrom) {
			//reinit chrom
			read_chrom(chrom);
			m_cur_chrom = chrom;
		}

		return(m_seq[0].size());
	}
	virtual void set_focus_interval(const string &chr, int fr, int to);
	void reset();

	void print_all_seq(){
		cerr << "DEBUG:: full seq  " << m_phylo.get_max_sp_node_id() << " species"  << endl;
		for (int i=0; i<m_phylo.get_max_sp_node_id(); i++){
			cerr << i <<":\t" ;
			cerr << m_seq[i] << endl;
		}
		cerr << "end sequences" << endl;
	}

	void set_ref_use_genome_coords(int use_ref_genome_coords){
		m_use_ref_genome_coords = use_ref_genome_coords;
	}

private:

	void read_chrom(const string &chrom);
	void read_chrom_ref_genome_coord(const string &chrom);
	void read_chrom_global_aln_coord(const string &chrom);
	void create_maps();
	void append_seq_by_maf_entry(MafParse& maf , int max_spi , int active_sp , int offset , int suffix_offset);
	void append_seq_of_N (string& str, const string& to_append);
	int get_num_active_sp (MafParse& maf , int max_spi);
	int get_ref_sp_id (MafParse& maf, int max_spi);
};

#endif // GenomesAlignment
