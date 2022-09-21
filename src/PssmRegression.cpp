#include "port.h"
BASE_CC_FILE
#include "LeastSquare.h"
#include "PssmRegression.h"

PssmRegression::PssmRegression(const vector<string> &loci, const string &init_mot, int bidir,
                               const ds_bitvec &relev, vector<float> &chip_stat, float eps,
                               int no_neg_mode, float min_rms_for_star)
    : m_bidirect(bidir), m_sequences(loci), m_relevant_vars(relev), m_pos_coefs(7), m_epsilon(eps),
      m_no_neg_mode(no_neg_mode), m_min_rms_for_star(min_rms_for_star) {

    m_interv_data.resize(m_relevant_vars.on_bits() + 1);
    m_coefs.resize(m_relevant_vars.on_bits() + 1, vector<float>(6));

    int vcount = 1;
    for (int vid = 0; vid < m_relevant_vars.size(); vid++) {
        if (m_relevant_vars[vid]) {
            m_interv_data[vcount] = chip_stat[vid];
            vcount++;
        }
    }
    m_chars.resize(init_mot.size(), vector<float>('T' + 1));

    int pos = 0;
    float unif = 0.02;
    float pref = 0.94;
    m_is_wildcard.resize(init_mot.size(), false);
    for (string::const_iterator i = init_mot.begin(); i != init_mot.end(); i++) {
        if (*i == '*') {
            fill(m_chars[pos].begin(), m_chars[pos].end(), 0.25);
            m_is_wildcard[pos] = true;
        } else {
            fill(m_chars[pos].begin(), m_chars[pos].end(), unif);
            m_chars[pos][*i] = pref;
        }
        pos++;
    }
}

float PssmRegression::go() {
    // iter on pos

    float prev_rms = _REAL(MAX) / 2;
    m_cur_rms = _REAL(MAX) / 3;
    int max_pos = m_chars.size();
    cerr << "starting regression iterations, epsilon = " << m_epsilon << endl;
    while ((prev_rms - m_cur_rms) > m_epsilon) {
        prev_rms = m_cur_rms;
        for (int pos = 0; pos < max_pos; pos++) {
            compute_wgts(pos);
            reopti_pos(pos);
            cerr << "pos " << pos << " RMS " << m_cur_rms << endl;
        }
        cerr << "delta in current iter = " << m_cur_rms - prev_rms;
    }
    int vcount = 1;
    int max_vid = m_relevant_vars.size();
    ofstream plot("regdone_plot.txt");
    for (int vid = 0; vid < max_vid; vid++) {
        if (!m_relevant_vars[vid]) {
            continue;
        }
        vector<float> &coefs = m_coefs[vcount];
        float x = coefs[1] * m_chars[max_pos - 1]['A'] + coefs[2] * m_chars[max_pos - 1]['G'] +
                  coefs[3] * m_chars[max_pos - 1]['C'] + coefs[4] * m_chars[max_pos - 1]['T'];
        plot << x << "\t" << -m_interv_data[vcount] << endl;
        vcount++;
    }
    return (m_cur_rms);
}

void PssmRegression::compute_wgts(int pos) {
    int vcount = 1;
    int max_vid = m_relevant_vars.size();
    for (int vid = 0; vid < max_vid; vid++) {
        if (!m_relevant_vars[vid]) {
            continue;
        }
        vector<float> &coefs = m_coefs[vcount];
        update_coefs(coefs, vid, pos);
        vcount++;
    }
}

void PssmRegression::update_coefs(vector<float> &coefs, int vid, int pos) {
    coefs[1] = 0;
    coefs[2] = 0;
    coefs[3] = 0;
    coefs[4] = 0;
    coefs[5] = 1; // constant

    string::const_iterator i = m_sequences[vid].begin();
    string::const_iterator max_i = m_sequences[vid].end() - m_chars.size();
    for (; i < max_i; i++) {
        string::const_iterator j = i;
        float prod = 1;
        int cur_pos = 0;
        char cur_nuc = '*';
        for (vector<vector<float>>::const_iterator p = m_chars.begin(); p < m_chars.end(); p++) {
            if (!(*j)) {
                prod = 0;
                break;
            }
            if (cur_pos != pos) {
                if (*j == 'N' || *j == '*') {
                    prod *= 0.25;
                } else {
                    prod *= (*p)[*j];
                }
            } else {
                cur_nuc = *j;
            }
            cur_pos++;
            j++;
        }
        if (pos == -1) {
            coefs[0] += prod;
        } else {
            if (prod != 0) {
                switch (cur_nuc) {
                case 'A':
                    coefs[1] += prod;
                    break;
                case 'G':
                    coefs[2] += prod;
                    break;
                case 'C':
                    coefs[3] += prod;
                    break;
                case 'T':
                    coefs[4] += prod;
                    break;
                default:
                    break;
                }
            }
        }
        if (m_bidirect) {
            float rprod = 1;
            j = i;
            int cur_pos = m_chars.size() - 1;
            for (vector<vector<float>>::reverse_iterator p = m_chars.rbegin(); p != m_chars.rend();
                 p++) {
                if (!(*j)) {
                    rprod = 0;
                    break;
                }
                if (cur_pos == pos) {
                    cur_nuc = *j;
                } else {
                    switch (*j) {
                    case 'A':
                        rprod *= (*p)['T'];
                        break;
                    case 'T':
                        rprod *= (*p)['A'];
                        break;
                    case 'C':
                        rprod *= (*p)['G'];
                        break;
                    case 'G':
                        rprod *= (*p)['C'];
                        break;
                    case '*':
                        rprod *= 0.25;
                        break;
                    case 'N':
                        rprod *= 0.25;
                        break;
                    default:
                        break;
                    }
                }
                cur_pos--;
                j++;
            }
            if (pos == -1) {
                coefs[0] += prod;
            } else {
                if (rprod != 0) {
                    switch (cur_nuc) {
                    case 'A':
                        coefs[4] += rprod;
                        break;
                    case 'G':
                        coefs[3] += rprod;
                        break;
                    case 'C':
                        coefs[2] += rprod;
                        break;
                    case 'T':
                        coefs[1] += rprod;
                        break;
                    default:
                        break;
                    }
                }
            }
        }
    }
    //	cerr << "coefs " << vid << " - " << coefs[1] << " " << coefs[2] << " " << coefs[3] << " " <<
    //coefs[4] << endl;
}

void PssmRegression::reopti_pos(int pos) {
    float chisq = 0;

    //	for(int i = 1; i < m_coefs.size(); i++) {
    //		cerr << "TEST\t" << pos << "\t" << m_coefs[i][1] << "\t" << m_coefs[i][2] << "\t" <<
    //m_coefs[i][3] << "\t" << m_coefs[i][4] << "\t" << m_interv_data[i] << endl;
    //
    //	}
    bool is_star = m_is_wildcard[pos];

    cerr << "prev coefs " << m_chars[pos]['A'] << " " << m_chars[pos]['G'] << " "
         << m_chars[pos]['C'] << " " << m_chars[pos]['T'] << endl;

    bool cont = true;
    int max_c = 5;
    vector<char> chr_list(6);
    chr_list[1] = 'A';
    chr_list[2] = 'G';
    chr_list[3] = 'C';
    chr_list[4] = 'T';
    chr_list[5] = 'N';
    while (cont) {
        cerr << "least square with max c " << max_c << endl;
        least_square(m_coefs, m_interv_data, m_pos_coefs, max_c, chisq);
        float rms = sqrt(chisq / m_coefs.size());
        if (is_star && (m_cur_rms - rms) < m_min_rms_for_star) {
            cerr << "Will not generalize star since rms delta is too small\n";
            return;
        }
        cerr << "rms was " << sqrt(chisq / m_coefs.size()) << endl;
        cont = false;
        for (int i = 1; i <= max_c; i++) {
            if (m_no_neg_mode && chr_list[i] != 'N' && m_pos_coefs[i] < 0) {
                cerr << "remove char " << chr_list[i] << "\t coef " << m_pos_coefs[i] << endl;
                char c = chr_list[i];
                chr_list[i] = chr_list[max_c];
                chr_list[max_c] = c;
                for (vector<vector<float>>::iterator j = m_coefs.begin(); j != m_coefs.end(); j++) {
                    (*j)[i] = (*j)[max_c];
                }
                max_c--;
                cont = true;
            }
        }
    }

    m_is_wildcard[pos] = false;
    m_cur_rms = sqrt(chisq / m_coefs.size());
    // float tot = m_pos_coefs[1] + m_pos_coefs[2] +
    // 				m_pos_coefs[3] + m_pos_coefs[4];

    fill(m_chars[pos].begin(), m_chars[pos].end(), 0);

    m_base_coef = 0;
    for (int i = 1; i <= max_c; i++) {
        if (chr_list[i] != 'N') {
            cerr << "set " << chr_list[i] << " to  " << m_pos_coefs[i] << endl;
            m_chars[pos][chr_list[i]] = m_pos_coefs[i];
        } else {
            m_base_coef = m_pos_coefs[i];
        }
    }
    //	m_chars[pos]['A'] = m_pos_coefs[1];
    //	m_chars[pos]['G'] = m_pos_coefs[2];
    //	m_chars[pos]['C'] = m_pos_coefs[3];
    //	m_chars[pos]['T'] = m_pos_coefs[4];

    cerr << "after opti pos " << chisq << endl;
    cerr << "new coefs " << m_chars[pos]['A'] << " " << m_chars[pos]['G'] << " "
         << m_chars[pos]['C'] << " " << m_chars[pos]['T'] << " const " << m_base_coef << endl;
}

float PssmRegression::test(const ds_bitvec &test_set, const vector<float> &chip_data,
                           ostream &out) {
    // test coefs
    int max_vid = m_relevant_vars.size();

    vector<float> coefs(6);
    float rms = 0;
    for (int vid = 0; vid < max_vid; vid++) {
        if (test_set[vid]) {
            coefs[0] = 0;
            update_coefs(coefs, vid, -1);
            out << coefs[0] << "\t" << chip_data[vid] << endl;
            float delta = (coefs[0] - chip_data[vid]);
            rms += delta * delta;
        }
    }
    rms /= test_set.on_bits();
    rms = sqrt(rms);
    return (rms);
}

void PssmRegression::update_yvals(const ds_bitvec &relevant_set, vector<float> &chip_data) {
    // test coefs
    int max_vid = m_relevant_vars.size();

    vector<float> coefs(6);
    for (int vid = 0; vid < max_vid; vid++) {
        if (relevant_set[vid]) {
            coefs[0] = 0;
            update_coefs(coefs, vid, -1);
            //			cerr << "Corr\t" << vid << "\t" << chip_data[vid] << "\t" << coefs[0] <<
            //endl;
            chip_data[vid] += coefs[0] + m_base_coef;
        }
    }
}

float PssmRegression::output_pssm(ostream &out, int psid) {
    int max_pos = m_chars.size();
    float alltot = 1;
    for (int pos = 0; pos < max_pos; pos++) {
        float tot = m_chars[pos]['A'] + m_chars[pos]['C'] + m_chars[pos]['G'] + m_chars[pos]['T'];
        alltot *= tot;

        out << psid << "\t" << pos << "\t" << m_chars[pos]['A'] / tot << "\t"
            << m_chars[pos]['C'] / tot << "\t" << m_chars[pos]['G'] / tot << "\t"
            << m_chars[pos]['T'] / tot << endl;
    }
    return (alltot);
}
