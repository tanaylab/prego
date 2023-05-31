#include <Rcpp.h>
#include <RcppParallel.h>
#include <string.h>
#include <vector>
using namespace std;

// [[Rcpp::plugins("cpp17")]]
struct KmerCounter : public RcppParallel::Worker {
    // source vector
    const Rcpp::StringVector sequences;

    // kmer length
    const std::size_t kmer_length;

    // output vector of maps
    std::vector<std::unordered_map<std::string, int>> &output_maps;

    // positions to consider within each sequence
    const int from_range;
    const int to_range;

    KmerCounter(const Rcpp::CharacterVector sequences, const std::size_t kmer_length,
               std::vector<std::unordered_map<std::string, int>> &output_maps,
               const int from_range, const int to_range)
        : sequences(sequences), kmer_length(kmer_length), output_maps(output_maps),
          from_range(from_range), to_range(to_range) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::unordered_map<std::string, int> string_map;
            std::string seq = Rcpp::as<std::string>(sequences[i]);
            std::string str = seq.substr(from_range, to_range - from_range + 1);

            for (std::size_t j = 0; j <= str.size() - kmer_length; j++) {
                std::string sub_str = str.substr(j, kmer_length);
                string_map[sub_str]++;
            }

            output_maps[i] = string_map;
        }
    }
};

// [[Rcpp::export]]
Rcpp::IntegerMatrix kmer_matrix_cpp(Rcpp::CharacterVector sequences,
                                    int kmer_length, int from_range = 0,
                                    Rcpp::Nullable<int> to_range = R_NilValue) {
    int nseq = sequences.size();
    // Initialize the output vector of maps with empty maps for each sequence
    std::vector<std::unordered_map<std::string, int>> output_maps(nseq);

    int to_range_val;
    if (to_range.isNull()) {
        std::string seq = Rcpp::as<std::string>(sequences[0]);
        to_range_val = seq.length();
    } else {
        to_range_val = Rcpp::as<int>(to_range);
    }

    KmerCounter counter(sequences, kmer_length, output_maps, from_range, to_range_val);

    RcppParallel::parallelFor(0, nseq, counter);

    // Get the unique kmers and their corresponding column indices
    std::unordered_map<std::string, int> kmer_indices;
    for (const auto &string_map : output_maps) {
        for (const auto &pair : string_map) {
            if (kmer_indices.find(pair.first) == kmer_indices.end()) {
                kmer_indices[pair.first] = kmer_indices.size();
            }
        }
    }

    // Now create the output matrix
    Rcpp::IntegerMatrix output(sequences.size(), kmer_indices.size());

    for (std::size_t i = 0; i < (size_t)sequences.size(); i++) {
        for (const auto &pair : output_maps[i]) {
            if (kmer_indices.find(pair.first) != kmer_indices.end()) {
                output(i, kmer_indices[pair.first]) = pair.second;
            }
        }
    }

    // Set the column names to the kmers
    Rcpp::CharacterVector kmers(kmer_indices.size());
    for (const auto &pair : kmer_indices) {
        kmers[pair.second] = pair.first;
    }
    Rcpp::colnames(output) = kmers;

    return output;
}
