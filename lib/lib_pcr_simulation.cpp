// Creator: Wenjie SUN
// Email: sunwjie at gmail.com
// Creation date: 2023/6/19
// Description:
//   This file contains the functions to perform PCR amplification and mutation.
//   The functions are used in the pcr_simulation function in the lib_pcr_simulation.R file.
// Modification history:
//   Wenjie SUN 2023/6/19 Add comments to the file.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

// Function Name: amplify
// Inputs:
//   tube_pcr (unordered_map*): A pointer to an unordered map that contains barcode sequences as the keys and their counts as the values.
//   efficiency_amplification (double): The efficiency of PCR amplification, which is the probability of a single molecule being successfully amplified per cycle.
// Output: None (void)
// Description:
//   This function performs PCR amplification on the input barcode sequences stored in tube_pcr.
//   The function simulates the amplification process by using a binomial distribution to randomly select the number of molecules that will be successfully amplified based on the efficiency_amplification parameter. The amplified count is then added to each corresponding barcode's existing clone count in tube_pcr
// Note:
//   This function modifies the input unordered map object directly since it takes a pointer to it as an argument rather than making a copy.
void amplify(unordered_map<string, long> *tube_pcr, double efficiency_amplification) {
    for (auto it=tube_pcr->begin(); it!=tube_pcr->end(); it++) {
        it->second += R::rbinom(it->second, efficiency_amplification);
    }
}


// The function to change a specific locus basepair in a sequence randomly
// Inputs:
//   s_seq (string*): A pointer to a string that contains the sequence.
//   locus (int): The locus of the basepair to be changed.
// Output: None (void)
// Description:
//   This function changes the basepair at the specified locus in the input sequence.
//   Each basepair has a 25% (equal) chance of being selected to replace the basepair at the specified locus.
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

// The function to mutate a locus in a sequence randomly
// Inputs:
//   s_seq (string): A string that contains the sequence.
// Output: A string that contains the mutated sequence.
// Description:
//   This function mutates a random locus in the input sequence.
//   The function first randomly selects a locus in the input sequence.
//   Then, the function calls the change_basepair function to change the basepair at the selected locus.
string mutate_seq(string s_seq) {
    // get the mutation location
    double r = ((double) rand() / (RAND_MAX));
    int locus = round(r * s_seq.size());
    // change the basepair
    change_basepair(&s_seq, locus);
    // return the sequence
    return s_seq;
}

// The function to mutate sequences during PCR process
// Inputs:
//   tube_pcr (unordered_map*): A pointer to an unordered map that contains barcode sequences as the keys and their counts as the values.
//   ratio_mutation (double): The chance of mutation happend per base per cycle during PCR process.
// Output: None (void)
// Description:
//   This function mutates the sequences in the input unordered map object.
//   The function first calculates the number of molecules that will be mutated based on the ratio_mutation parameter.
//   Then, the function calls the mutate_seq function to mutate the sequences.
//   The function then inserts the mutated sequences into the input unordered map object.
// Note:
//   This function modifies the input unordered map object directly since it takes a pointer to it as an argument rather than making a copy.
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
// The function to perform PCR amplification and mutation
// Inputs:
//   temp (List): A list that contains barcode sequences and their counts.
//   cycle (int): The number of PCR cycles.
//   efficiency_amplification (double): The efficiency of PCR amplification.
//   ratio_mutation (double): The chance of mutation happend per base per cycle during PCR process.
// Output: A list that contains barcode sequences and their counts after PCR amplification and mutation.
// Description:
//   This function performs PCR amplification and mutation on the input barcode sequences.
//   The function first converts the input list object into an unordered map object.
//   Then, the function calls the amplify and mutate functions to perform PCR amplification and mutation.
//   The function then converts the unordered map object back to a list object and returns it.
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
