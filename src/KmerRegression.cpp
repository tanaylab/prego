#include "ProgressReporter.h"
#include <Rcpp.h>
#include <string.h>
#include <vector>
using namespace std;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
Rcpp::NumericMatrix kmer_matrix_cpp(Rcpp::CharacterVector sequences, Rcpp::CharacterVector kmers,
                                    int from_range = 0, Rcpp::Nullable<int> to_range = R_NilValue) {
    int nseq = sequences.size();
    int nkmer = kmers.size();
    Rcpp::NumericMatrix freq(nseq, nkmer);
    std::string dna[4] = {"A", "C", "G", "T"};
    ProgressReporter progress;

    progress.init(nseq, 5);
    for (int i = 0; i < nseq; i++) {
        for (int j = 0; j < nkmer; j++) {
            std::string seq = Rcpp::as<std::string>(sequences[i]);
            std::string kmer = Rcpp::as<std::string>(kmers[j]);

            int to_range_val;
            if (to_range.isNull()) {
                to_range_val = seq.length();
            } else {
                to_range_val = Rcpp::as<int>(to_range);
            }

            for (size_t pos = from_range; pos <= to_range_val - kmer.length(); pos++) {
                if (kmer.find("N") != std::string::npos) {
                    for (const std::string &base : dna) {
                        std::string replaced_kmer = kmer;
                        std::replace(replaced_kmer.begin(), replaced_kmer.end(), 'N', base[0]);
                        if (seq.substr(pos, replaced_kmer.length()) == replaced_kmer) {
                            freq(i, j)++;
                        }
                    }
                } else if (seq.substr(pos, kmer.length()) == kmer) {
                    freq(i, j)++;
                }
            }
        }
        progress.report(1);
    }

    colnames(freq) = kmers;

    return freq;
}
