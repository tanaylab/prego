#include <math.h>

#include "SeqQual.h"

float SeqQual::s_qual2base_prob[256];
float SeqQual::s_qual2mis_prob[256];
float SeqQual::s_log_qual2base_prob[256];
float SeqQual::s_log_qual2mis_prob[256];

SeqQual::Qual2prob_init SeqQual::s_qual2prob_init;

SeqQual::Qual2prob_init::Qual2prob_init()
{
	for (unsigned char qual = 0; qual < 0xff; ++qual) {
		// precision cannot drop below 0.25
		s_qual2base_prob[qual] = max((float)0.25, 1 - powf(10, -0.1 * qual));
		s_qual2mis_prob[qual] = (1 - s_qual2base_prob[qual]) / 3.;
		s_log_qual2base_prob[qual] = log(s_qual2base_prob[qual]);
		s_log_qual2mis_prob[qual] = log(s_qual2mis_prob[qual]);
	}
}

