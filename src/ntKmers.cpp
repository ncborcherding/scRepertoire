// 2-bit-based nucleotide kmer counting - unoptimized
// could use a kmercounter class with an uint_fast64_t[128] for toNtIndex instead of the switch statement
// by Qile Yang

#include <Rcpp.h>
#include <vector>
#include <string>
#include "scRepHelper.h"

inline unsigned short int toNtIndex(const char nt) {
    switch(nt) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default: return 3;
    }
}

constexpr char Nts[4] = {'A', 'C', 'G', 'T'};

inline char lastNt(unsigned int index) {
    return Nts[index & 3];
}

inline std::string toNtKmer(unsigned long int index, int k) {
    std::string kmer = "";
    for (int i = 0; i < k; i++) {
        kmer = lastNt(index) + kmer;
        index >>= 2;
    }
    return kmer;
}

inline bool isNt(char c) {
    switch(c) {
        case 'A': case 'C': case 'G': case 'T': return true;
        default: return false;
    }
}

// [[Rcpp::export]]
Rcpp::CharacterVector rcppGenerateUniqueNtMotifs(int k) {
    long int numKmers = 1 << (k + k);
    Rcpp::CharacterVector motifs (numKmers);
    for (long int i = 0; i < numKmers; i++) {
        motifs[i] = toNtKmer(i, k);
    }
    return motifs;
}

inline void updateSkip(int& skip, const char c, const int k) {
    if (!isNt(c)) {
        skip = k;
    } else if (skip > 0) {
        skip--;
    }
}

inline bool updateSkipAndReturnIfShouldntSkip(int& skip, const char c, const int k) {
    updateSkip(skip, c, k);
    return skip == 0;
}

// actual kmer counter - doesnt handle _NA_ for k = 1
inline void kmerCount(std::vector<double>& bins, const unsigned int mask, const std::string& seq, const int k) {
    
    int n = (int) seq.size();
    if (n < k) {
        return;
    }

    int skip = 0;
    unsigned long int kmer = 0;

    for (int i = 0; i < (k - 1); i++) { // this segment to initialize the kmer should be deletable if skip is adjusted?
        kmer = (kmer << 2) | toNtIndex(seq[i]);
        updateSkip(skip, seq[i], k);
    }

    for (int i = k - 1; i < n; i++) {
        kmer = ((kmer << 2) & mask) | toNtIndex(seq[i]);
        if (updateSkipAndReturnIfShouldntSkip(skip, seq[i], k)) {
            bins[kmer]++;
        }
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector rcppGetNtKmerPercent(const std::vector<std::string>& seqs, const int k) {
    const unsigned int mask = (1 << (k + k)) - 1;
    int numKmers = mask + 1;
    std::vector<double> bins (numKmers, 0);

    for (std::string seq : seqs) {
        kmerCount(bins, mask, seq, k);
    }

    long double binSum = scRepHelper::sum(bins);
    if (binSum == 0.0) { // pretty sure this can only happen if there arent valid seqs?
        return Rcpp::NumericVector (numKmers, R_NaReal);
    }
    
    double scaleFactor = 1 / binSum;
    for (int i = 0; i < numKmers; i++) {
        bins[i] *= scaleFactor;
    }

    return scRepHelper::convertZerosToNA(bins, numKmers);
}
