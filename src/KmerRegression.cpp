#include "ProgressReporter.h"
#include <Rcpp.h>
#include <RcppParallel.h>
#include <string.h>
#include <vector>
using namespace std;

// [[Rcpp::plugins("cpp17")]]
struct KmerCounter : public RcppParallel::Worker {
    const Rcpp::StringVector sequences;
    const Rcpp::StringVector kmers;
    const int from_range;
    const int to_range_val;
    std::string dna[4] = {"A", "C", "G", "T"};
    RcppParallel::RMatrix<double> result;

    KmerCounter(const Rcpp::StringVector sequences, const Rcpp::StringVector kmers, int from_range,
                int to_range_val, Rcpp::NumericMatrix result)
        : sequences(sequences), kmers(kmers), from_range(from_range), to_range_val(to_range_val),
          result(result) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (size_t j = 0; j < (size_t)kmers.length(); j++) {
                std::string seq = Rcpp::as<std::string>(sequences[i]);
                std::string kmer = Rcpp::as<std::string>(kmers[j]);

                for (size_t pos = from_range; pos <= to_range_val - kmer.length(); pos++) {
                    if (kmer.find("N") != std::string::npos) {
                        for (const std::string &base : dna) {
                            std::string replaced_kmer = kmer;
                            std::replace(replaced_kmer.begin(), replaced_kmer.end(), 'N', base[0]);
                            if (seq.substr(pos, replaced_kmer.length()) == replaced_kmer) {
                                result(i, j)++;
                            }
                        }
                    } else if (seq.substr(pos, kmer.length()) == kmer) {
                        result(i, j)++;
                    }
                }
            }
        }
    }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix kmer_matrix_cpp(Rcpp::CharacterVector sequences,
                                             Rcpp::CharacterVector kmers, int from_range = 0,
                                             Rcpp::Nullable<int> to_range = R_NilValue) {
    int nseq = sequences.size();
    int nkmer = kmers.size();
    Rcpp::NumericMatrix freq(nseq, nkmer);

    int to_range_val;
    if (to_range.isNull()) {
        std::string seq = Rcpp::as<std::string>(sequences[0]);
        to_range_val = seq.length();
    } else {
        to_range_val = Rcpp::as<int>(to_range);
    }   

    KmerCounter counter(sequences, kmers, from_range, to_range_val, freq);
    RcppParallel::parallelFor(0, nseq, counter);

    colnames(freq) = kmers;

    return freq;
}
