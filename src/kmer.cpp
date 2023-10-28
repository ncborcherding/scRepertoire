// nucleotide kmer counting

#include <Rcpp.h>
#include <vector>
#include <string>

inline const unsigned short int toIndex(const char x) {
    switch(x) {
        case 'A':case 'a':
            return 0;
        case 'C':case 'c':
            return 1;
        case 'G':case 'g':
            return 2;
        default:
            return 3;
    }
}

inline const std::string toKmer(int index, int k) {
    std::string kmer = "";
    for (int i = 0; i < k; i++) {
        int baseIndex = index % 4; 
        char base;

        if (baseIndex == 0) {
            base = 'A';
        } else if (baseIndex == 1) {
            base = 'C';
        } else if (baseIndex == 2) {
            base = 'G';
        } else {
            base = 'T';
        }

        kmer = base + kmer;
        index /= 4;
    }
    return kmer;
}

// [[Rcpp::export]]
Rcpp::CharacterVector rcppGenerateUniqueNtMotifs(int k) {
    unsigned int numKmers = 1 << (k + k);
    Rcpp::CharacterVector motifs (numKmers);
    for (unsigned int i = 0; i < numKmers; i++) {
        motifs[i] = toKmer(i, k);
    }
    return motifs;
}

inline void kmerCount(std::vector<long double>& bins, const unsigned int& mask, const std::string& seq, const int k) {
    unsigned int kmer = 0;
    for (int i = 0; i < k - 1; i++) {
        kmer = (kmer << 2) | toIndex(seq[i]);
    }
    for (int i = k - 1; i < seq.size(); i++) {
        kmer = ((kmer << 2) & mask) | toIndex(seq[i]);
        bins[kmer]++;
    }
}

const long double sum(std::vector<long double>& v) {
    long double n = 0;
    for (long double num : v) {
        n += num;
    }
    return n;
}

Rcpp::NumericVector convertToNumericVector(std::vector<long double>& v, int n) {
    Rcpp::NumericVector converted (n);
    for (int i = 0; i < n; i++) {
        converted[i] = v[i];
    }
    return converted;
}

// [[Rcpp::export]]
const Rcpp::NumericVector rcppGetNtKmerPercent(const std::vector<std::string>& seqs, const int k) {
    const unsigned int mask = (1 << (k + k)) - 1;
    int numKmers = mask + 1;
    std::vector<long double> bins (numKmers);

    for (std::string seq : seqs) {
        kmerCount(bins, mask, seq, k);
    }

    long double sumOfNonZeroTerms = sum(bins);
    if (sumOfNonZeroTerms < 0.0001) { // == 0
        return Rcpp::NumericVector (numKmers, R_NaReal);
    }
    
    const double scaleFactor = 1 / sumOfNonZeroTerms;
    for (int i = 0; i < mask + 1; i++) {
        bins[i] *= scaleFactor;
    }

    return convertToNumericVector(bins, numKmers);
}
