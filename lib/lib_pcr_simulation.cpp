#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

void amplify(unordered_map<string, long> *tube_pcr, double efficiency_amplification) {
    for (auto it=tube_pcr->begin(); it!=tube_pcr->end(); it++) {
        it->second += R::rbinom(it->second, efficiency_amplification);
    }
}

int change_basepair(string *s_seq, int locus) {
    double r = ((double) rand() / (RAND_MAX));
    char basepair = (*s_seq)[locus];
    if (r < 0.25) {
        if (basepair != 'T') {
            (*s_seq)[locus] = 'T';
            return 0;
        } else {
            r += 0.25;
        }
    }
    if (r < 0.5) {
        if (basepair != 'C') {
            (*s_seq)[locus] = 'C';
            return 0;
        } else {
            r += 0.25;
        }
    }
    if (r < 0.75) {
        if (basepair != 'G') {
            (*s_seq)[locus] = 'G';
            return 0;
        } else {
            r += 0.25;
        }
    }
    if (r <= 1) {
        if (basepair != 'A') {
            (*s_seq)[locus] = 'A';
            return 0;
        } else {
            (*s_seq)[locus] = 'T';
        }
    }
    return 0;
}

string mutate_seq(string s_seq) {
    // get the mutation location
    double r = ((double) rand() / (RAND_MAX));
    int locus = round(r * s_seq.size());
    // change the basepair
    change_basepair(&s_seq, locus);
    // return the sequence
    return s_seq;
}

void mutate(unordered_map<string, long> *tube_pcr, double ratio_mutation) {
    for (auto it=tube_pcr->begin(); it!=tube_pcr->end(); it++) {
        // calculate mutation sequence
        string s_seq = it->first;
        double n_seq_mutant = R::rbinom(s_seq.size() * it->second, ratio_mutation);
        it->second -= n_seq_mutant;
        for (auto i=0; i<n_seq_mutant; i++) {
            // do the mutation, if the sequence exist
            string s_seq_new = mutate_seq(s_seq);
            auto it_mut = tube_pcr->find(s_seq_new);
            if (it_mut == tube_pcr->end()) {
                //   insert the sequence if it is not exist
                (*tube_pcr).insert(make_pair(s_seq_new, 1));
            } else {
                //   add the frequency to the exist sequence
                it_mut->second++;
            }
        }
        
    }
}

// [[Rcpp::export]]
List pcr_amplify(List temp, int cycle, double efficiency_amplification, double ratio_mutation) {

    CharacterVector seq_temp = temp["seq"];
    IntegerVector freq_temp = temp["freq"];

    unordered_map<string, long> tube_pcr;

    for (auto i=0; i<seq_temp.size(); i++) {
        tube_pcr.insert(make_pair(seq_temp[i], freq_temp[i]));
    }

    for (auto i=0; i<cycle; i++) {
        Rcout<<"PCR cycle "<<i<<endl;
        amplify(&tube_pcr, efficiency_amplification);
        mutate(&tube_pcr, ratio_mutation);
    }

    CharacterVector seq_res(tube_pcr.size());
    NumericVector freq_res(tube_pcr.size());

    int i=0;
    for (auto it = tube_pcr.begin(); it != tube_pcr.end(); it++) {
        seq_res[i] = it->first;
        freq_res[i] = it->second;
        i++;
    }

    return List::create(_["seq"] = seq_res, _["freq"] = freq_res);
}

/*** R
*/
