#include "port.h"
BASE_CC_FILE
#include "strutil.h"
#include "dnastrutil.h"
#include "Random.h"

const int dinuc_revcomp[] = {
15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0
};
const char revcomp[] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 't', 0, 'g', 0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 'a', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

const int nuc_revcomp[] = {3,2,1,0};

bool dnapat_match(const string &mot, string::const_iterator targ)
{
	for(string::const_iterator mot_i = mot.begin();
	    mot_i != mot.end(); 
	    mot_i++) {
		if(*mot_i != '*' && *mot_i != *targ) {
			return(false);
		}
		targ++;
	}
	return(true);
}
int count_dnapat_match(string::const_iterator i0,
			string::const_iterator i_to, const string &mot)
{
	string openmot;
	for(string::const_iterator i = mot.begin(); i != mot.end(); i++) {
		if(*i < 'A') {
			for(int c = 0; c < (*i - '0'); c++) {
				openmot += "*";
			}
		} else {
			openmot += *i;
		}
	}
	int count = 0;
	for(; i0 != i_to; i0++) {
		string::const_iterator i = i0;
		string::const_iterator j = openmot.begin(); 
		for(; j != openmot.end(); j++) {
			if(*j != '*' && *j != *i) {
				break;
			}
			i++;
			if(i == i_to) {
				break;
			}
		}
		if(j == openmot.end()) {
			count++;
		}
	}
	return(count);
}
int count_dnapat_match(const string &targ, const string &mot)
{
	string openmot;
	for(string::const_iterator i = mot.begin(); i != mot.end(); i++) {
		if(*i < 'A') {
			for(int c = 0; c < (*i - '0'); c++) {
				openmot += "*";
			}
		} else {
			openmot += *i;
		}
	}
	int count = 0;
	for(string::const_iterator i0 = targ.begin(); i0 != targ.end(); i0++) {
		string::const_iterator i = i0;
		string::const_iterator j = openmot.begin(); 
		for(; j != openmot.end(); j++) {
			if(*j != '*' && *j != *i) {
				break;
			}
			i++;
			if(i == targ.end()) {
				break;
			}
		}
		if(j == openmot.end()) {
			count++;
		}
	}
	return(count);
}

bool dna_is_poly_x(const string &targ)
{
	if(targ.length() == 0) {
		return(true);
	}
	for(string::const_iterator i = targ.begin() + 1; i != targ.end(); i++) {
		if(*i != *(i-1)) {
			return(false);
		}
	}
	return(true);
}

bool dna_is_poly_xy(const string &targ)
{
	if(targ.length() < 2) {
		return(true);
	}
	for(string::const_iterator i = targ.begin() + 2; i != targ.end(); i++) {
		if(*i != *(i-2)) {
			return(false);
		}
	}
	return(true);
}

//compute the best nongap duplicate of start..start+gap in start-dup_horiz..start+gap+dup_horiz
//and put the hamming dist in dist and the relative pos in where.
//where will be -1..-horiz if the best hit is to the left
//0..horiz-1 if the best hit is to the right
//horiz..horiz+gap if the best hit is derived by wrap around at start+where

void find_dup(string::const_iterator start, int gap, 
	int &dist, int &where, 
	int &hit_l, int &hit_pos, int dup_horiz)
{

/*	cerr << "look for gap " << gap << " ";
	for(int i = 0; i < gap; i++) {
		cerr << *(start+i);
	}
	cerr << endl;*/
	static const int UCDIFF = 'a'-'A';

	dist = gap + 1;
	hit_l = 0;
	hit_pos = -1;

	string::const_iterator wrap1 = start - gap + 1;
	string::const_iterator wrap2 = start + 1;
	for(int k = 0; k < gap - 1; k++) {
		string::const_iterator i = wrap1;
		string::const_iterator j = wrap2;
		int cdist = 0;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < wrap1 + gap) {
			if(*i != *j && *i != *j+UCDIFF && *i != *j-UCDIFF) {
				cdist++;
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
			hit_pos = dup_horiz + k;
		}
		if(cdist < dist) {
			dist = cdist;
			where = dup_horiz + k;
		}
		wrap1++;
		wrap2++;
	}
	string::const_iterator source = start - gap;
	string::const_iterator horiz = start - dup_horiz - gap;
	string::const_iterator end = start + gap;

	while(source > horiz) {
		string::const_iterator i = start;
		string::const_iterator j = source;
		int cdist = 0;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < end) {

			if(*i != *j && *i != *j+UCDIFF && *i != *j-UCDIFF) {
				cdist++;
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
			hit_pos = source - start + gap - 1;
		}
		if(cdist < dist) {
			dist = cdist;
			where = source - start + gap - 1;
		}
		source--;
	}

	horiz = start + dup_horiz + gap;
	source = start + gap;
	
	while(source < horiz) {
		string::const_iterator i = start;
		string::const_iterator j = source;
		int cdist = 0;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < end) {
			if(*i != *j) {
				cdist++;
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
			hit_pos = source - end;
		}
		if(cdist < dist) {
			dist = cdist;
			where = source - end;
		}
		source++;
	}
	
/*	cerr << "best was dist " << dist << " ";
	for(int i = 0; i < gap; i++) {
		if(where > 0) {
			cerr << *(end+where+i);
		} else {
			cerr << *(start+where-gap+i);
		}
	}
	cerr << endl;*/
}

void find_rc_dup(string::const_iterator start, int gap, 
	int &dist, int &where, 
	int &hit_l, int &hit_pos, int dup_horiz)
{

/*	cerr << "look for gap " << gap << " ";
	for(int i = 0; i < gap; i++) {
		cerr << *(start+i);
	}
	cerr << endl;*/
	static const int UCDIFF = 'a'-'A';

	dist = gap + 1;
	hit_l = 0;
	hit_pos = -1;
	string rcseq;

	for(string::const_iterator i = start + gap - 1; i >= start; i--) {
		rcseq.push_back(revcomp[(unsigned char)(*i)]);
	}
	
	string::const_iterator source = start - gap;
	string::const_iterator horiz = start - dup_horiz - gap;
	string::const_iterator end = rcseq.end();

	while(source > horiz) {
		string::const_iterator i = rcseq.begin();
		string::const_iterator j = source;
		int cdist = 0;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < end) {

			if(*i != *j && *i != *j+UCDIFF && *i != *j-UCDIFF) {
				cdist++;
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
			hit_pos = source - start + gap - 1;
		}
		if(cdist < dist) {
			dist = cdist;
			where = source - start + gap - 1;
		}
		source--;
	}

	horiz = start + dup_horiz + gap;
	source = start + gap;
	
	while(source < horiz) {
		string::const_iterator i = rcseq.begin();
		string::const_iterator j = source;
		int cdist = 0;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < end) {
			if(*i != *j) {
				cdist++;
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
			hit_pos = source - start - gap;
		}
		if(cdist < dist) {
			dist = cdist;
			where = source - start - gap;
		}
		source++;
	}
}

void find_perfect_submatch(string::const_iterator pat, 
		string::const_iterator base, 
		int l, int search_l, int &hit_l)
{
	static const int UCDIFF = 'a'-'A';

	hit_l = 0;
	string::const_iterator source = base;
	string::const_iterator horiz = base + search_l;
	string::const_iterator end = pat + l;

	while(source < horiz) {
		string::const_iterator i = pat;
		string::const_iterator j = source;
		int cur_hit_l = 0;
		int best_hit_l = 0;
		while(i < end) {

			if(*i != *j && *i != *j+UCDIFF && *i != *j-UCDIFF) {
				if(cur_hit_l > best_hit_l) {
					best_hit_l = cur_hit_l;
				}
				cur_hit_l = 0;
			} else {
				cur_hit_l++;
			}
			i++;
			j++;
		}
		if(cur_hit_l > best_hit_l) {
			best_hit_l = cur_hit_l;
		}
		if(best_hit_l > hit_l) {
			hit_l = best_hit_l;
		}
		source++;
	}
}

void find_match(string::const_iterator start, string::const_iterator targ,
			int gap, int &dist, int &where, int dup_horiz)
{

	string::const_iterator source = start - gap;
	string::const_iterator horiz = start - dup_horiz - gap;
	string::const_iterator end = targ + gap;

	dist = gap + 1;

	while(source > horiz) {
		string::const_iterator i = targ;
		string::const_iterator j = source;
		int cdist = 0;
		while(i < end) {
			if(*i != *j) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - start + gap;
		}
		source--;
	}

	horiz = start + dup_horiz;
	source = start;
	
	while(source < horiz) {
		string::const_iterator i = targ;
		string::const_iterator j = source;
		int cdist = 0;
		while(i < end) {
			if(*i != *j) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - start;
		}
		source++;
	}
}

void find_match_uc(string::const_iterator start, string::const_iterator targ,
			int gap, int &dist, int &where, int dup_horiz)
{

	static const int UCDIFF = 'a'-'A';

	string::const_iterator source = start - gap;
	string::const_iterator horiz = start - dup_horiz - gap;
	string::const_iterator end = targ + gap;

	dist = gap + 1;

	while(source > horiz) {
		string::const_iterator i = targ;
		string::const_iterator j = source;
		int cdist = 0;
		while(i < end) {
			if((*i - *j) % UCDIFF != 0) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - start + gap;
		}
		source--;
	}

	horiz = start + dup_horiz;
	source = start;
	
	while(source < horiz) {
		string::const_iterator i = targ;
		string::const_iterator j = source;
		int cdist = 0;
		while(i < end) {
			if((*i - *j) % UCDIFF != 0) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j++;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - start;
		}
		source++;
	}
}

void find_max_stem(string::const_iterator start,  string::const_iterator end,
				int &stem, int &where1, int &where2)
{

	string::const_iterator cstart = start;

	stem = 0;

	while(cstart < end - 6) {
		string::const_iterator cend = end - 1;
		while(cend > cstart + 6) {
			string::const_iterator i = cstart;
			string::const_iterator j = cend;
			while(i < j-1) {
				if(*i != revcomp[(unsigned char)(*j)]) {
					break;
				}
				i++;
				j--;
			}
			if(i - cstart > stem) {
				stem = i - cstart;
				where1 = cstart- start;
				where2 = cend - start;
			}
			cend--;
		}
		cstart++;
	}
}


float get_gapless_identity(const string &seq1, const string &seq2)
{
	float same = 0;
	int n = 0;
	string::const_iterator j = seq2.begin();
	for(string::const_iterator i = seq1.begin(); i != seq1.end(); i++) {
		if(*j != '-' && *i != '-') {
			if(*j == *i) {
				same++;
			}
			n++;
		}
		j++;
	}
	return(same/n);
}

int repeat_dist(string::const_iterator min, string::const_iterator max, int offset, int l)
{
	string::const_iterator i = min + offset;
	string::const_iterator j = min + offset;

	int dist = 0;
	for(int k = 0; k < l; k++) {
		if((*j) > 'Z') {
			return(0);
		}
		j++;
	}
	while(i >= min && j < max
	&& (*i) < 'Z' && (*j) < 'Z') {
		dist++;
		i--;
		j++;
	}
	return(dist);
}
void center_stem(string::const_iterator center, int &dist)
{

	dist = 0;

	string::const_iterator i = center;
	string::const_iterator j = center-1;
	
	while(*i == revcomp[(unsigned char)(*j)]) {
		dist++;
		i++;
		j--;
	}
}
void inward_stem(string::const_iterator left, string::const_iterator right, int &dist)
{

	dist = 0;

	while(*left == revcomp[(unsigned char)(*right)]) {
		dist++;
		left++;
		right--;
	}
}
void dna_match(string::const_iterator left, string::const_iterator right, int &dist)
{
	static const int UCDIFF = 'a'-'A';

	dist = 0;

	while((*left - *right) % UCDIFF == 0) {
		dist++;
		left++;
		right++;
	}
}
void dna_match_rev(string::const_iterator left, string::const_iterator right, int &dist)
{

	static const int UCDIFF = 'a'-'A';
	dist = 0;

	while((*left - *right) % UCDIFF == 0) {
		dist++;
		left--;
		right--;
	}
}
void rpt_stat(string::const_iterator min, string::const_iterator max, 
				int offset, int l, int &rpt_dist, int horiz)
{
	string::const_iterator i = min + offset;
	string::const_iterator j = min + offset;

	rpt_dist = -1;
	for(int k = 0; k < l; k++) {
		if((*j) > 'Z') {
			return;
		}
		j++;
	}
	rpt_dist = 0;
	while(rpt_dist < horiz && i >= min && j < max) {
		if((*i) > 'Z' || (*j) > 'Z') {
			break;
		}
		rpt_dist++;
		i--;
		j++;
	}
}

//the gap is assimed to be in seq1 !
void best_gap_positions(const string &s1, const string &s2,
	int cent_pos, int horiz, int gap_length,  
	list<int> &best_pos, int &best_mis) 
{
	best_mis = 2*horiz + 1;

	const int UC_DIFF = 'a' - 'A';

/*	int min_pos = 0;
	if(pos > horiz) {
		min_pos = (horiz - pos);
	}
	int max_pos = 2*horiz;
	if(max_pos - min_pos > s1.length()) {
		max_pos  = (s1.length() + min_pos);
	}
	string::const_iterator seq1 = s1.begin() + pos - horiz + min_pos;
	string::const_iterator seq2 = s2.begin() + pos - horiz + min_pos;*/

	int max_pos = horiz;
	string::const_iterator seq1 = s1.begin() + cent_pos;
	while(seq1 < s1.end() && max_pos < 2*horiz) {
		seq1++;
		if(*seq1 != '-') {
			max_pos++;
		}
	}
	int min_pos = horiz;
	seq1 = s1.begin() + cent_pos;
	string::const_iterator seq2 = s2.begin() + cent_pos;
	while(seq1 > s1.begin() && min_pos > 0) {
		seq1--;
		seq2--;
		if(*seq1 != '-') {
			min_pos--;
		}
	}

//	cerr << "will scan, minpos " << min_pos << " max p " << max_pos << " l " << gap_length << " maf_l = " << s1.length() << " pos " << cent_pos << endl;

	for(int pos = min_pos; pos < max_pos; pos++) {
//		cerr << "pos = " << pos << endl;
		string::const_iterator j = seq2;
		string::const_iterator i = seq1; 
		while(*i == '-') {
			i++;
		}
		while(*j == '-') {
			j++;
		}
		int mis = 0;
		int cnt = min_pos;
		while(cnt < pos) {
			if((*j - *i) % UC_DIFF != 0) {
				mis++;
				if(mis > best_mis) {
					break;
				}
			} 
			j++;
			i++;
			cnt++;
			while(*i == '-') {
				i++;
			}
			while(*j == '-') {
				j++;
			}
		}
//		cerr << "mis before " << pos << " was " << mis << " bst " << best_mis << endl;
		if(mis > best_mis) {
			continue;
		}
		int add = 0;	//jump gap positions
		while(add < gap_length) {
			j++;
			while(*j == '-') {
				j++;
			}
			add++;
		}
		while(cnt < max_pos) {
			if((*j - *i) % UC_DIFF != 0) {
				mis++;
				if(mis > best_mis) {
					break;
				}
			}
			j++;
			i++;
			cnt++;
			while(*i == '-') {
				i++;
			}
			while(*j == '-') {
				j++;
			}
		}
//		cerr << "mis aft " << pos << " was " << mis << " bst " << best_mis << endl;
		if(mis < best_mis) {
			best_pos.clear();	
			best_mis = mis;
			best_pos.push_back(pos);
		} else if(mis == best_mis) {
			best_pos.push_back(pos);
		}
	}
}

int get_dinuc_code(string::const_iterator i)
{
	int cd = 0;
	switch(*i) {
		case 'A': cd = 0; break;
		case 'a': cd = 0; break;
		case 'C': cd = 4; break;
		case 'c': cd = 4; break;
		case 'G': cd = 8; break;
		case 'g': cd = 8; break;
		case 'T': cd = 12; break;
		case 't': cd = 12; break;
		default: return(-1);
	}
	switch(*(i+1)) {
		case 'A': break;
		case 'a': break;
		case 'C': cd += 1; break;
		case 'c': cd += 1; break;
		case 'G': cd += 2; break;
		case 'g': cd += 2; break;
		case 'T': cd += 3; break;
		case 't': cd += 3; break;
		default: return(-1);
	}
	return(cd);
}
int get_trinuc_code(string::const_iterator i)
{
	int cd = 0;
	switch(*i) {
		case 'A': cd = 0; break;
		case 'a': cd = 0; break;
		case 'C': cd = 16; break;
		case 'c': cd = 16; break;
		case 'G': cd = 32; break;
		case 'g': cd = 32; break;
		case 'T': cd = 48; break;
		case 't': cd = 48; break;
		default: return(-1);
	}
	switch(*(i+1)) {
		case 'A': cd += 0; break;
		case 'a': cd += 0; break;
		case 'C': cd += 4; break;
		case 'c': cd += 4; break;
		case 'G': cd += 8; break;
		case 'g': cd += 8; break;
		case 'T': cd += 12; break;
		case 't': cd += 12; break;
		default: return(-1);
	}
	switch(*(i+2)) {
		case 'A': break;
		case 'a': break;
		case 'C': cd += 1; break;
		case 'c': cd += 1; break;
		case 'G': cd += 2; break;
		case 'g': cd += 2; break;
		case 'T': cd += 3; break;
		case 't': cd += 3; break;
		default: return(-1);
	}
	return(cd);
}

int get_dinuc_rc(int din)
{
	return(dinuc_revcomp[din]);
}

void get_best_tandems(string::const_iterator base, int l, const list<int> &cands, int &best_miss, list<int> &tands)
{
	int best_match = 0;
	for(list<int>::const_iterator offs = cands.begin();
	    offs != cands.end();
	    offs++) {
		int match5 = 0;
		int match3 = 0;
		for(string::const_iterator i = base + *offs - l; i < base + *offs; i++) {
			if(*i == *(i+l)) {
				match5++;
			}
		}
		for(string::const_iterator i = base + *offs; i < base + *offs + l; i++) {
			if(*i == *(i+l)) {
				match3++;
			}
		}
		int match = (match5 < match3 ? match3 : match5);
		if(match > best_match) {
			best_match = match;
			tands.clear();
			tands.push_back(*offs);
		} else if(match == best_match) {
			tands.push_back(*offs);
		}
	}
	best_miss = l - best_match;
}

void synonm_code_randomize(const string &seq, int from, int to, 
			int strand, int horiz, float gc_cont, string &targ)
{

	//first copy to template:

	cerr << "will randomize " << from << " " << to << endl;
	targ.resize(2*horiz + to - from);
	copy(seq.begin() + from - horiz, seq.begin() + to + horiz, targ.begin());

	float gc_thresh3 = (gc_cont/2)/(1-gc_cont/2);
	for(string::iterator i = targ.begin() + horiz; 
	    i < targ.end() - horiz; i+=3) {
		//get codon, randomize alternatives
		string::iterator j = i+1;
		string::iterator k = i+2;
	//	cerr << "Codon = " << *i << *j << *k << endl;
		if(*i > 'Z') {
			*i -= 'z'-'Z';
		}
		if(*j > 'Z') {
			*j -= 'z'-'Z';
		}
		if(*k > 'Z') {
			*k -= 'z'-'Z';
		}
		float r = Random::fraction();
		int mode = -1;
		switch (*i) {
		case 'A':
			switch (*j) {
			case 'A':
//Lys 	AAA, AAG
//Asn 	AAU, AAC 	
				mode = 2;
			case 'C':
//Thr 	ACU, ACC, ACA, ACG
				mode = 4;
			case 'G':
//Arg 	AGA, AGG
//Ser 	AGU, AGC
				mode = 2;
				break;
			case 'T':
//Met 	AUG/Start
//Ile 	AUU, AUC, AUA 	
				if(*k != 'G') {
					if(r < gc_thresh3) {
						*k = 'C';
					} else if(Random::fraction() > 0.5) {
						*k = 'T';
					} else {
						*k = 'A';
					}
				}
				mode = -1;
				break;
			default: break;
			}
			break;
		case 'C':
			switch (*j) {
			case 'A':
//Gln 	CAA, CAG 	
//His 	CAU, CAC 	
				mode = 2;
				break;
			case 'C':
//Pro 	CCU, CCC, CCA, CCG
				mode = 4;
				break;
			case 'G':
//Arg 	CGU, CGC, CGA, CGG
				mode = 4;
				break;
			case 'T':
//Leu 	CUU, CUC, CUA, CUG
				mode = 4;
				break;
			default: break;
			}
			break;
		case 'G':
			switch (*j) {
			case 'A':
//Asp 	GAU, GAC 	
//Glu 	GAA, GAG 	
				mode = 2;
				break;
			case 'C':
//Ala 	GCU, GCC, GCA, GCG 	
				mode = 4;
				break;
			case 'G':
//Gly 	GGU, GGC, GGA, GGG 	
				mode = 4;
				break;
			case 'T':
//Val 	GUU, GUC, GUA, GUG
				mode = 4;
				break;
			default: break;
			}
			break;
		case 'T':
			switch (*j) {
			case 'A':
//Tyr 	UAU, UAC
//STOP 	UAG, UAA
				mode = 2;
				break;
			case 'C':
//Ser 	UCU, UCC, UCA, UCG
				mode = 4;
				break;
			case 'G':
//Cys 	UGU, UGC 	
//Trp 	UGG
//UGA	stop
				if(*k == 'T' || *k == 'C') {
					if(r > 0.5) {
						*k = 'T';
					} else {
						*k = 'C';
					}
				}
				break;
			case 'T':
//Phe 	UUU, UUC
//Leu	UUA, UUG
				mode = 2;
				break;
			default: break;
			}
			break;
		}
		if(mode == 4) {
			if(r < 0.5) {
				if(r > gc_cont/2) {
					*k = 'A';
				} else {
					*k = 'C';
				}
			} else {
				if(r > 0.5 + gc_cont/2) {
					*k = 'T';
				} else {
					*k = 'G';
				}
			}
		} else if(mode == 2) {
			if(*k == 'A' || *k == 'G') {
				if(r>gc_cont) {
					*k = 'A';
				} else {
					*k = 'G';
				}
			} else if(*k == 'T' || *k == 'C') {
				if(r>gc_cont) {
					*k = 'T';
				} else {
					*k = 'C';
				}
			}
		}
//		cerr << "Rand = " << *i << *j << *k << endl;
	}
}

/*void find_rc_dup(string::const_iterator start, int gap, 
			int &dist, int &where, int dup_horiz)
{

	string::const_iterator source = start - gap;
	string::const_iterator horiz = start - dup_horiz - gap;
	string::const_iterator end = start + gap;

	dist = gap + 1;
	
	while(source > horiz) {
		string::const_iterator i = start;
		string::const_iterator j = source+gap;
		int cdist = 0;
		while(i < end) {
			if(*i != revcomp[*j]) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j--;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - start + gap;
		}
		source--;
	}

	horiz = start + dup_horiz + gap;
	source = start + gap;
	
	while(source < horiz) {
		string::const_iterator i = start;
		string::const_iterator j = source + gap;
		int cdist = 0;
		while(i < end) {
			if(*i != revcomp[*j]) {
				cdist++;
			}
			if(cdist > dist) {
				break;
			}
			i++;
			j--;
		}
		if(cdist <= dist) {
			dist = cdist;
			where = source - end;
		}
		source++;
	}
}*/

int get_nuc_rc(int nuc)
{
	return nuc_revcomp[nuc];
}