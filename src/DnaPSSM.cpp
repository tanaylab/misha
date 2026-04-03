#include <cstdint>
#include "port.h"
#include "DnaPSSM.h"
#include "Random.h"
#include <R_ext/Arith.h>

#include <algorithm>

// O(1) lookup tables for base encoding — replaces switch statements in hot loops.
// BASE_ENCODE:      A/a->0, C/c->1, G/g->2, T/t->3, else->-1
// COMPLEMENT_ENCODE: A/a->3(T), C/c->2(G), G/g->1(C), T/t->0(A), else->-1
// NEUTRAL_CHAR:     1 for N/n/*, 0 otherwise

#define X -1
const int8_t DnaLookupTables::BASE_ENCODE[256] = {
//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 0-15
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 16-31
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 32-47
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 48-63
    X,  0,  X,  1,  X,  X,  X,  2,  X,  X,  X,  X,  X,  X,  X,  X,  // 64-79  (A=65->0, C=67->1, G=71->2)
    X,  X,  X,  X,  3,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 80-95  (T=84->3)
    X,  0,  X,  1,  X,  X,  X,  2,  X,  X,  X,  X,  X,  X,  X,  X,  // 96-111 (a=97->0, c=99->1, g=103->2)
    X,  X,  X,  X,  3,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 112-127 (t=116->3)
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 128-143
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 144-159
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 160-175
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 176-191
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 192-207
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 208-223
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 224-239
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // 240-255
};

const int8_t DnaLookupTables::COMPLEMENT_ENCODE[256] = {
//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  3,  X,  2,  X,  X,  X,  1,  X,  X,  X,  X,  X,  X,  X,  X,  // A->3(T), C->2(G), G->1(C)
    X,  X,  X,  X,  0,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // T->0(A)
    X,  3,  X,  2,  X,  X,  X,  1,  X,  X,  X,  X,  X,  X,  X,  X,  // a->3, c->2, g->1
    X,  X,  X,  X,  0,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  // t->0
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
    X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
};
#undef X

const int8_t DnaLookupTables::NEUTRAL_CHAR[256] = {
//  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  // '*'=42->1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  // 'N'=78->1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  // 'n'=110->1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
};


void DnaProbVec::normalize()
{
	float sum = m_p[0] + m_p[1] + m_p[2] + m_p[3];

	m_p[0] /= sum;
	m_p[1] /= sum;
	m_p[2] /= sum;
	m_p[3] /= sum;

	if(m_p[0] != 0) {
		m_logp[0] = log(m_p[0]);
	} else {
		m_logp[0] = R_NegInf;
	}
	if(m_p[1] != 0) {
		m_logp[1] = log(m_p[1]);
	} else {
		m_logp[1] = R_NegInf;
	}
	if(m_p[2] != 0) {
		m_logp[2] = log(m_p[2]);
	} else {
		m_logp[2] = R_NegInf;
	}
	if(m_p[3] != 0) {
		m_logp[3] = log(m_p[3]);
	} else {
		m_logp[3] = R_NegInf;
	}
}
void DnaProbVec::normalize_log()
{
	float sum = m_logp[0];
       	log_sum_log(sum, m_logp[1]);
       	log_sum_log(sum, m_logp[2]);
       	log_sum_log(sum, m_logp[3]);

	// cerr << "normalize, sum = " << sum << " 0 " << m_logp[0] << " 1 " << m_logp[1] << endl;

	m_logp[0] -= sum;
	m_logp[1] -= sum;
	m_logp[2] -= sum;
	m_logp[3] -= sum;

	m_p[0] = exp(m_logp[0]);
	m_p[1] = exp(m_logp[1]);
	m_p[2] = exp(m_logp[2]);
	m_p[3] = exp(m_logp[3]);
}

float DnaProbVec::get_entropy() const
{
	return(-m_p[0]*m_logp[0]
		- m_p[1]*m_logp[1]
		- m_p[2]*m_logp[2]
		- m_p[3]*m_logp[3]);
}

float DnaProbVec::get_max_log_prob() const
{
	float logp = m_logp[0];
	if(m_logp[1] > logp) {
		logp = m_logp[1];
	}
	if(m_logp[2] > logp) {
		logp = m_logp[2];
	}
	if(m_logp[3] > logp) {
		logp = m_logp[3];
	}

	return(logp);
}

ostream &operator<<(ostream &out, DnaProbVec &pvec)
{
	out << int(pvec.get_prob('A')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('C')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('G')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('T')*1000)/1000.0 << endl;
	return(out);
}

ostream &operator<<(ostream &out, const DnaProbVec &pvec)
{
	out << int(pvec.get_prob('A')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('C')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('G')*1000)/1000.0 << "\t"
	<< int(pvec.get_prob('T')*1000)/1000.0 << endl;
	return(out);
}


const DnaPSSM &DnaPSSM::operator=(const DnaPSSM &other)
{
	m_chars = other.m_chars;
	m_min_range = other.m_min_range;
	m_max_range = other.m_max_range;
	m_bidirect = other.m_bidirect;
	return(*this);
}
void DnaPSSM::resize(int sz)
{
	m_chars.resize(sz);
}

void DnaPSSM::init_from_seed(const string &seed, float prior)
{
	m_chars.resize(seed.length());
	vector<DnaProbVec>::iterator p = m_chars.begin();
	vector<float> back(4, prior);
	for(string::const_iterator i = seed.begin(); i != seed.end(); i++) {
		p->reset(back);
		switch(*i) {
			case 'a':
			case 'A': p->set_direct_prob(0, 1 - 3*prior); break;
			case 'c':
			case 'C': p->set_direct_prob(1, 1 - 3*prior); break;
			case 'g':
			case 'G': p->set_direct_prob(2, 1 - 3*prior); break;
			case 't':
			case 'T': p->set_direct_prob(3, 1 - 3*prior); break;
		}
		p->normalize();
		p++;
	}
}

float DnaPSSM::get_max_ll() const
{
	float logp = 0;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		logp += p->get_max_log_prob();
	}
	return(logp);

}

void DnaPSSM::calc_like(const string &target, float &logp) const
{
	string::const_iterator i = target.begin();
	logp = 0;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*i];
		if(code < 0) {
			logp = R_NegInf;
			return;
		}
		logp += p->get_log_prob_from_code(code);
		i++;
	}
}
void DnaPSSM::calc_like_rc(const string &target, float &logp) const
{
	string::const_iterator i = target.begin();
	logp = 0;
	for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
	    p != m_chars.rend();
	    p++) {
		int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*i];
		if(code < 0) {
			logp = R_NegInf;
			return;
		}
		logp += p->get_log_prob_from_code(code);
		i++;
	}
}
void DnaPSSM::calc_like(string::const_iterator &j, float &logp) const
{
	logp = 0;
	string::const_iterator i = j;
	for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
	    p < m_chars.end();
	    p++) {
		int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*i];
		if(code < 0) {
			logp = R_NegInf;
			return;
		}
		logp += p->get_log_prob_from_code(code);
		i++;
	}
}

static const float c_log_quarter = -1.38629;

void DnaPSSM::calc_like_rc(string::const_iterator &j, float &logp) const
{
	logp = 0;
	string::const_iterator i = j;
	for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
	    p != m_chars.rend();
	    p++) {
		int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*i];
		if(code < 0) {
			logp = R_NegInf;
			return;
		}
		logp += p->get_log_prob_from_code(code);
		i++;
	}
}

string::const_iterator DnaPSSM::max_like_match(const string &target,
				float &best_logp, int &best_dir, const bool &combine_strands) const
{
	if(target.length() < m_chars.size()) {
		best_logp = R_NegInf;
		return(target.begin());
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	string::const_iterator best_pos;
	best_logp = R_NegInf;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += p->get_avg_log_prob();
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}

		float total_logp = logp;

		if(m_bidirect) {
			float rlogp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					rlogp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					rlogp += p->get_avg_log_prob();
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						rlogp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}

			if (combine_strands) {
				// Combine probabilities by adding them in log space
				log_sum_log(total_logp, rlogp);
				best_dir = 0; // Indicate combined strands
			} else {
				// select best strand
				if(rlogp > logp) {
					total_logp = rlogp;
					best_dir = -1;
				} else {
					best_dir = 1;
				}
			}
		} else {
			best_dir = 1;
		}

		if(total_logp > best_logp) {
			best_logp = total_logp;
			best_pos = i;
		}
	}
	return(best_pos);
}


void DnaPSSM::update_like_vec(const string &target,
			vector<float> &likes, vector<float> &deltas,
			vector<int1> &dirs)
{
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator delta = deltas.begin() + m_min_range;
	vector<float>::iterator like = likes.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += c_log_quarter;
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
		*dir = 1;
		if(m_bidirect) {
			float rlogp = 0;
			j = i;
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					rlogp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					rlogp += p->get_avg_log_prob();
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						rlogp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
			if(rlogp > logp) {
				logp = rlogp;
				*dir = -1;
			}
		}
		if(logp == R_NegInf) {
			*delta = R_NegInf;
			*like = R_NegInf;
		} else {
			*delta = -(*like);
			*delta += logp;
			*like = logp;
		}
		like++;
		delta++;
		dir++;
	}
}

void DnaPSSM::integrate_like_seg(const char *min_i, const char *max_i, float &energy) const
{
	energy = R_NegInf;
	for(const char *i = min_i;
	    i <= max_i;
	    i++) {
		const char *j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += p->get_avg_log_prob();
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
   		log_sum_log(energy, logp);
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp += p->get_avg_log_prob();
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
       		log_sum_log(energy, logp);
		}
	}
}
void DnaPSSM::integrate_like(const string &target, float &energy, vector<float> *spat_dist) const
{
	if(target.length() < m_chars.size()) {
		energy = R_NegInf;
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	energy = R_NegInf;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += c_log_quarter;
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
   		log_sum_log(energy, logp);
		if(spat_dist) {
			log_sum_log((*spat_dist)[i - target.begin()], logp);
		}
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp += p->get_avg_log_prob();
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
   			log_sum_log(energy, logp);
			if(spat_dist) {
				log_sum_log((*spat_dist)[i - target.begin()], logp);
			}
		}
	}
}

void DnaPSSM::count(string::const_iterator seq, float weight, int dir)
{
	if(dir == 1) {
		for(vector<DnaProbVec>::iterator i = m_chars.begin();
		    i != m_chars.end();
		    i++) {
			i->incr_weight(*seq, weight);
			++seq;
		}
	} else {
		for(vector<DnaProbVec>::reverse_iterator i = m_chars.rbegin();
		    i != m_chars.rend();
		    i++) {
			int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*seq];
			if(code >= 0) {
				i->direct_incr_weight(code, weight);
			}
			++seq;
		}
	}
}

void DnaPSSM::count_weighted(const string &target, vector<float> &wgts,
				vector<int1> &dirs, float thresh_wgt)
{
	//iterate on the correct range
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator wgt = wgts.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		//if logp is very small - ignore it
		if(*wgt < thresh_wgt) {
			wgt++;
			dir++;
			continue;
		}
		string::const_iterator j = i;
		if(*dir == 1) {
			for(vector<DnaProbVec>::iterator p = m_chars.begin();
			    p < m_chars.end();
			    p++) {
				if((*j) && *j != 'N' && *j !='*' && *j != 'n') {
					p->incr_weight(*j, *wgt);
				}
				j++;
			}
		} else {
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					p->direct_incr_weight(code, *wgt);
				}
				j++;
			}
		}
		wgt++;
		dir++;
	}
}
void DnaPSSM::count_log_weighted(const string &target, vector<float> &wgts,
				vector<int1> &dirs, float thresh_wgt)
{
	//iterate on the correct range
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	vector<float>::iterator wgt = wgts.begin() + m_min_range;
	vector<int1>::iterator dir = dirs.begin() + m_min_range;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		//if logp is very small - ignore it
		if(*wgt < thresh_wgt) {
			wgt++;
			dir++;
			continue;
		}
		string::const_iterator j = i;
		if(*dir == 1) {
			for(vector<DnaProbVec>::iterator p = m_chars.begin();
			    p < m_chars.end();
			    p++) {
				if((*j) && *j != 'N' && *j !='*' && *j != 'n') {
					p->incr_log_weight(*j, *wgt);
				}
				j++;
			}
		} else {
			for(vector<DnaProbVec>::reverse_iterator p =
							m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					p->direct_incr_log_weight(code, *wgt);
				}
				j++;
			}
		}
		wgt++;
		dir++;
	}
}

void DnaPSSM::normalize()
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->normalize();
	}
}
void DnaPSSM::normalize_logs()
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->normalize_log();
	}
}

void DnaPSSM::add_dirichlet_prior(float prior)
{
	if (prior <= 0)
		return;

	vector<float> freqs(4);
	for (vector<DnaProbVec>::iterator i = m_chars.begin();
	     i != m_chars.end();
	     ++i) {
		freqs[0] = i->get_direct_prob(0) + prior;
		freqs[1] = i->get_direct_prob(1) + prior;
		freqs[2] = i->get_direct_prob(2) + prior;
		freqs[3] = i->get_direct_prob(3) + prior;
		i->reset(freqs);
		i->normalize();
	}
}

void DnaPSSM::reset_prior(const vector<float> &prior)
{
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		i->reset(prior);
	}
}

//currently assuming same length profiles
float DnaPSSM::dot_product(DnaPSSM &arg)
{
	vector<DnaProbVec>::iterator j = arg.m_chars.begin();
	float prod = 1;
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		prod *= (*i).dot(*j);
		j++;
	}
	return(prod);
}
float DnaPSSM::log_dot_product(DnaPSSM &arg)
{	
	vector<DnaProbVec>::iterator j = arg.m_chars.begin();
	float prod = 1;
	for(vector<DnaProbVec>::iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		prod *= (*i).dot(*j);
		j++;
	}
	return(log(prod));
}

void DnaPSSM::permut_randomize()
{
	int max_i = m_chars.size();
	for(int count = 0; count < max_i*2; count++) {
		int i = int(Random::fraction() * max_i);
		int j = int(Random::fraction() * max_i);
		DnaProbVec tmp = m_chars[i];
		m_chars[i] = m_chars[j];
		m_chars[j] = tmp;
	}
}

void DnaPSSM::write_tab(ostream &pssmd, int id) const
{
	int pos = 0;
	for(vector<DnaProbVec>::const_iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		pssmd << id << "\t" << pos << "\t" << *i;
		pos++;
	}
}

const float DnaPSSM::CONSENSUS_SINGLE_THRESH = 0.6;
const float DnaPSSM::CONSENSUS_DOUBLE_THRESH = 0.85;

string DnaPSSM::get_consensus() const
{
	string output;
	vector<int> ps(4);
	for(vector<DnaProbVec>::const_iterator i = m_chars.begin();
	    i != m_chars.end();
	    i++) {
		ps[0] = int(1000 * i->get_prob('A'))*4;
		ps[1] = int(1000 * i->get_prob('C'))*4 + 1;
		ps[2] = int(1000 * i->get_prob('G'))*4 + 2;
		ps[3] = int(1000 * i->get_prob('T'))*4 + 3;

		sort(ps.begin(), ps.end());

		if(ps[3] > CONSENSUS_SINGLE_THRESH*4000) {
			int code = ps[3] % 4;
			switch(code) {
				case 0: output += 'A';
					break;
				case 1: output += 'C';
					break;
				case 2: output += 'G';
					break;
				case 3: output += 'T';
					break;
				default:
					break;
			}
			continue;
		}
		if(ps[3] + ps[2] >= 4000*CONSENSUS_DOUBLE_THRESH) {
			int code = (ps[3]%4) * 4 + ps[2]%4;
			switch(code) {
				case 1: output += 'M';	//AC
					break;
				case 2: output += 'R';	//AG
					break;
				case 3: output += 'W';	//AT
					break;
				case 4: output += 'M';	//CA
					break;
				case 6: output += 'S';	//CG
					break;
				case 7: output += 'Y';	//CT
					break;
				case 8: output += 'R';	//GA
					break;
				case 9: output += 'S';	//GC
					break;
				case 11: output += 'K';	//GT
					break;
				case 12: output += 'W';	//TA
					break;
				case 13: output += 'Y';	//TC
					break;
				case 14: output += 'K';	//TG
					break;
				default:
					output += 'e';
					break;
			}
			continue;
		}
		output += "*";
	}
	return(output);
}

ostream &operator<<(ostream &out, const DnaPSSM &pssm)
{
	// cerr << "[" << pssm.get_min_range() << "," << pssm.get_max_range() << "] dir=" << pssm.is_bidirect() << endl;
	for(uint64_t i = 0; i < pssm.size(); i++) {
		out << pssm[i];
	}
	out << endl;
	for(uint64_t i = 0; i < pssm.size(); i++) {
		out << pssm[i].get_log_prob('A') << "\t" << pssm[i].get_log_prob('C') << "\t" << pssm[i].get_log_prob('G') << "\t" << pssm[i].get_log_prob('T') << endl;
	}
	return(out);
}

void DnaPSSM::integrate_energy(const string &target, float &energy, vector<float> &spat_func, int spat_bin_size) const
{
	if(target.length() < m_chars.size()) {
		energy = R_NegInf;
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	energy = R_NegInf;
	int pos = 0;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		int spat_bin = int(pos/spat_bin_size);
		// Clamp to valid range
		if(spat_bin >= (int)spat_func.size()) {
			spat_bin = spat_func.size() - 1;
		}
		pos++;
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += p->get_avg_log_prob();
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
		logp += log(spat_func[spat_bin]);
		log_sum_log(energy, logp);
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator
							p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp += c_log_quarter;
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
			logp += log(spat_func[spat_bin]);
			log_sum_log(energy, logp);
		}
	}
}

void DnaPSSM::integrate_energy_logspat(const string &target, float &energy, vector<float> &spat_log_func, int spat_bin_size) const
{
	if(target.length() < m_chars.size()) {
		energy = R_NegInf;
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	energy = R_NegInf;
	int pos = 0;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		int spat_bin = int(pos/spat_bin_size);
		// Clamp to valid range
		if(spat_bin >= (int)spat_log_func.size()) {
			spat_bin = spat_log_func.size() - 1;
		}
		pos++;

		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += p->get_avg_log_prob();
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
		logp += spat_log_func[spat_bin];
		log_sum_log(energy, logp);

		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp += c_log_quarter;
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
			logp += spat_log_func[spat_bin];
			log_sum_log(energy, logp);
		}
	}
}

void DnaPSSM::integrate_energy_max_logspat(const string &target, float &energy, vector<float> &spat_log_func, int spat_bin_size) const
{
	if(target.length() < m_chars.size()) {
		energy = R_NegInf;
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}

	energy = R_NegInf;

	int pos = 0;
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		int spat_bin = int(pos/spat_bin_size);
		// Clamp to valid range
		if(spat_bin >= (int)spat_log_func.size()) {
			spat_bin = spat_log_func.size() - 1;
		}
		pos++;

		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p < m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += p->get_avg_log_prob();
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			j++;
		}
		logp += spat_log_func[spat_bin];

		if(m_bidirect) {
			float logp_rev = 0;
			j = i;
			for(vector<DnaProbVec>::const_reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp_rev = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp_rev += c_log_quarter;
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp_rev += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
			logp_rev += spat_log_func[spat_bin];
			log_sum_log(logp, logp_rev);
		}

		if(logp > energy) {
			energy = logp;
		}
	}
}

void DnaPSSM::like_thresh_match(const string &target, float thresh,
		list<int> &poss, list<float> &vals, list<int> &dirs)
{
	if(target.length() < m_chars.size()) {
		return;
	}

	string::const_iterator max_i = target.begin() + m_max_range;
	if(max_i > target.end() - m_chars.size()) {
		max_i = target.end() - m_chars.size();
	}
	for(string::const_iterator i = target.begin() + m_min_range;
	    i <= max_i;
	    i++) {
		string::const_iterator j = i;
		float logp = 0;
		for(vector<DnaProbVec>::const_iterator p = m_chars.begin();
		    p != m_chars.end();
		    p++) {
			if(!(*j)) {
				logp = R_NegInf;
				break;
			}
			if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
				logp += c_log_quarter;
			} else {
				int code = DnaLookupTables::BASE_ENCODE[(unsigned char)*j];
				if(code >= 0) {
					logp += p->get_log_prob_from_code(code);
				}
			}
			if(logp < thresh) {
				break;
			}
			j++;
		}
		if(logp > thresh) {
			poss.push_back(i - target.begin());
			dirs.push_back(1);
			vals.push_back(logp);
		}
		if(m_bidirect) {
			logp = 0;
			j = i;
			for(vector<DnaProbVec>::reverse_iterator p = m_chars.rbegin();
			    p != m_chars.rend();
			    p++) {
				if(!(*j)) {
					logp = R_NegInf;
					break;
				}
				if(DnaLookupTables::NEUTRAL_CHAR[(unsigned char)*j]) {
					logp += c_log_quarter;
				} else {
					int code = DnaLookupTables::COMPLEMENT_ENCODE[(unsigned char)*j];
					if(code >= 0) {
						logp += p->get_log_prob_from_code(code);
					}
				}
				j++;
			}
			if(logp > thresh) {
				poss.push_back(i - target.begin());
				dirs.push_back(-1);
				vals.push_back(logp);
			}
		}
	}
}
