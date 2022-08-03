#ifndef SEQQUAL_H_
#define SEQQUAL_H_

#include <math.h>
#include <string>
#include <vector>

using namespace std;

class SeqQual {
public:
	// phred_base can be either 33 for Phred+33 scores (Illumina 1.8+ format) or 64 for Phread+64 (Illumina 1.3+)
	SeqQual() {}
	SeqQual(const string &quality, int phred_base) { init(quality, phred_base); }

	void init(const string &quality, int phred_base);

	const vector<unsigned char> &quality() const { return m_quality; }
	float log_precision() const { return m_log_precision; }
	float precision() const { return m_precision; }
	float base_prob(int pos) const { return s_qual2base_prob[m_quality[pos]]; }
	float mis_prob(int pos) const { return s_qual2mis_prob[m_quality[pos]]; }
	float log_base_prob(int pos) const { return s_log_qual2base_prob[m_quality[pos]]; }
	float log_mis_prob(int pos) const { return s_log_qual2mis_prob[m_quality[pos]]; }

private:
	vector<unsigned char> m_quality; // m_quality[i] is the Phred score of a base at position i
	float                 m_precision;
	float                 m_log_precision;

	static float s_qual2base_prob[256];  // translates the quality at position i to probability of base being correct
	static float s_qual2mis_prob[256];   // translates the quality at position i to probability of mismatch
	static float s_log_qual2base_prob[256];
	static float s_log_qual2mis_prob[256];

	struct Qual2prob_init {
		Qual2prob_init();
	};

	static Qual2prob_init s_qual2prob_init;
};


//---------------------------- IMPLEMENTATION ----------------------------------

inline void SeqQual::init(const string &quality, int phred_base)
{
	m_log_precision = 0;
	m_quality.resize(quality.length());
	for (unsigned i = 0; i < quality.length(); ++i) {
		m_quality[i] = quality[i] - phred_base;
		m_log_precision += log_base_prob(i);
	}
	m_precision = exp(m_log_precision);
}

#endif /* SEQQUAL_H_ */
