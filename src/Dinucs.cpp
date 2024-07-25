#include "port.h"

#include <Rcpp.h>
#include <algorithm>
#include <string.h>
#include <unordered_map>
#include <vector>
using namespace std;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
Rcpp::DataFrame n_nuc_distribution(Rcpp::StringVector sequences, int n, int size = 1000) {
    if (size < 1) {
        Rcpp::stop("Size must be at least 1");
    }

    if (n < 1) {
        Rcpp::stop("n must be at least 1");
    }

    // Prepare a vector of unordered_maps to hold n-nucleotide frequencies at each position
    std::vector<std::unordered_map<std::string, int>> nucleotide_counts(size);

    // Convert StringVector to vector of strings
    std::vector<std::string> seqs = Rcpp::as<std::vector<std::string>>(sequences);

    // Iterate over each sequence
    for (const auto &sequence : seqs) {
        // Validate that the sequence length is at least 'size'
        if ((int)sequence.size() < size) {
            Rcpp::stop("All sequences must have a length of at least 'size'");
        }

        // In each sequence, iterate over each position and extract the n-nucleotide
        for (int pos = 0; pos < size - n + 1; ++pos) {
            std::string n_nucleotide = sequence.substr(pos, n);

            // Update the counts in the unordered_map for the corresponding n-nucleotide and
            // position
            nucleotide_counts[pos][n_nucleotide]++;
        }
    }

    // After iterating over all sequences and positions, calculate the percentages
    int total_sequences = sequences.size();

    // Generate all possible n-nucleotides
    std::vector<std::string> all_n_nucleotides;
    std::function<void(std::string, int)> generate_n_nucleotides = [&](std::string current,
                                                                       int depth) {
        if (depth == n) {
            all_n_nucleotides.push_back(current);
            return;
        }
        for (char base : {'A', 'C', 'G', 'T'}) {
            generate_n_nucleotides(current + base, depth + 1);
        }
    };
    generate_n_nucleotides("", 0);

    // Create an Rcpp data frame with the required structure and fill it with the calculated
    // percentages
    int num_n_nucleotides = std::pow(4, n);
    Rcpp::NumericMatrix percentage_matrix(size, num_n_nucleotides);

    for (int pos = 0; pos < size; ++pos) {
        for (int n_nuc_i = 0; n_nuc_i < num_n_nucleotides; ++n_nuc_i) {
            std::string n_nucleotide = all_n_nucleotides[n_nuc_i];
            percentage_matrix(pos, n_nuc_i) =
                (double)nucleotide_counts[pos][n_nucleotide] / total_sequences;
        }
    }

    // Create the DataFrame
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("pos") = Rcpp::seq(1, size));

    // Add columns for each n-nucleotide
    for (int i = 0; i < num_n_nucleotides; ++i) {
        df.push_back(percentage_matrix(Rcpp::_, i), all_n_nucleotides[i]);
    }

    return df;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calc_sequences_dinuc_cpp(std::vector<std::string> sequences) {
    // Define all possible dinucleotides
    std::vector<std::string> dinucs = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                                       "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"};

    // Initialize result matrix
    int num_seqs = sequences.size();
    int num_dinucs = dinucs.size();
    Rcpp::NumericMatrix result(num_seqs, num_dinucs);

    // Process each sequence
    for (int i = 0; i < num_seqs; ++i) {
        std::unordered_map<std::string, int> dinuc_counts;
        const std::string &seq = sequences[i];

        // Count dinucleotides in the sequence
        for (size_t j = 0; j < seq.length() - 1; ++j) {
            std::string dinuc = seq.substr(j, 2);
            dinuc_counts[dinuc]++;
        }

        // Fill the result matrix
        for (int j = 0; j < num_dinucs; ++j) {
            result(i, j) = dinuc_counts[dinucs[j]];
        }
    }

    // Set column names
    Rcpp::CharacterVector colnames(num_dinucs);
    for (int i = 0; i < num_dinucs; ++i) {
        colnames[i] = dinucs[i];
    }
    result.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), colnames);

    return result;
}