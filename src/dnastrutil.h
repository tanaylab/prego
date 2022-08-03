#ifndef seqpack_strutil_h
#define seqpack_strutil_h 1

#include <vector>
#include <string>
#include <list>

bool dnapat_match(const string &mot, string::const_iterator targ);
int count_dnapat_match(string::const_iterator i0, string::const_iterator to_i, const string &mot);
int count_dnapat_match(const string &targ, const string &mot);

bool dna_is_poly_x(const string &targ);
bool dna_is_poly_xy(const string &targ);

void find_dup(string::const_iterator start, int gap, 
				int &dist, int &where, int &hit_l, int &hit_pos, int dup_horiz);
void find_rc_dup(string::const_iterator start, int gap, 
	int &dist, int &where, 
	int &hit_l, int &hit_pos, int dup_horiz);
void find_max_stem(string::const_iterator start,  string::const_iterator end,
				int &stem, int &where1, int &where2);

void find_match(string::const_iterator start, string::const_iterator targ,
			int gap, int &dist, int &where, int dup_horiz);
void find_match_uc(string::const_iterator start, string::const_iterator targ,
			int gap, int &dist, int &where, int dup_horiz);

void center_stem(string::const_iterator start, int &dist);
void inward_stem(string::const_iterator left, string::const_iterator right, int &dist);
float get_gapless_identity(const string &seq1, const string &seq2);

int repeat_dist(string::const_iterator min, string::const_iterator max, int offset, int l);
void rpt_stat(string::const_iterator min, string::const_iterator max, 
			int offset, int l, int &rpt_dist, int horiz);
void dna_match_rev(string::const_iterator left, string::const_iterator right, int &dist);
void dna_match(string::const_iterator left, string::const_iterator right, int &dist);
void best_gap_positions(const string &s1, const string &s2,
	int pos, int horiz, int gap_length,  
	list<int> &best_pos, int &best_mis);

void find_perfect_submatch(string::const_iterator pat, 
		string::const_iterator base, 
		int l, int search_l, int &hit_l);

int get_dinuc_code(string::const_iterator i);
int get_trinuc_code(string::const_iterator i);
int get_dinuc_rc(int din);
void get_best_tandems(string::const_iterator base, int l, const list<int> &cands, int &best_miss, list<int> &tands);

void synonm_code_randomize(const string &seq, int from, int to, 
			int strand, int horiz, float gccont, string &targ);

int get_nuc_rc(int nuc);

#endif // seqpack_strutil_h
