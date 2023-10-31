// bit-based nucleotide kmer counting

#include <Rcpp.h>
#include <vector>
#include <string>

inline const unsigned short int toIndex(const char x) {
    switch(x) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default: return 3;
    }
}

inline const char lastNt(unsigned int index) {
    switch(index & 3) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        default: return 'T';
    }
}

inline const std::string toKmer(unsigned int index, int k) {
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
    unsigned int numKmers = 1 << (k + k);
    Rcpp::CharacterVector motifs (numKmers);
    for (unsigned int i = 0; i < numKmers; i++) {
        motifs[i] = toKmer(i, k);
    }
    return motifs;
}

inline void updateSkip(int& skip, char c, int k) {
    if (!isNt(c)) {
        skip = k;
    } else if (skip > 0) {
        skip--;
    }
}

// doesnt handle _NA_ for k = 1
inline void kmerCount(std::vector<long double>& bins, const unsigned int mask, const std::string& seq, int k) {
    int skip = 0;
    unsigned int kmer = 0;

    for (int i = 0; i < k - 1; i++) {
        kmer = (kmer << 2) | toIndex(seq[i]);
        updateSkip(skip, seq[i], k);
    }

    for (int i = k - 1; i < (int) seq.size(); i++) {
        kmer = ((kmer << 2) & mask) | toIndex(seq[i]);

        updateSkip(skip, seq[i], k);
        if (skip == 0) {bins[kmer]++;}
    }
}

const long double sum(std::vector<long double>& v) {
    long double n = 0;
    for (long double num : v) {
        n += num;
    }
    return n;
}

Rcpp::NumericVector convertZerosToNA(std::vector<long double>& v, int numKmers) {
    Rcpp::NumericVector converted (numKmers, R_NaReal);
    for (int i = 0; i < numKmers; i++) {
        if (v[i] > 0) {
            converted[i] = v[i];
        }
    }
    return converted;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcppGetNtKmerPercent(const std::vector<std::string>& seqs, const int k) {
    const unsigned int mask = (1 << (k + k)) - 1;
    int numKmers = mask + 1;
    std::vector<long double> bins (numKmers, 0);

    for (std::string seq : seqs) {
        kmerCount(bins, mask, seq, k);
    }

    long double binSum = sum(bins);
    if (binSum == 0.0) { // pretty sure this can only happen if there arent valid seqs?
        return Rcpp::NumericVector (numKmers, R_NaReal);
    }
    
    const double scaleFactor = 1 / binSum;
    for (int i = 0; i < numKmers; i++) {
        bins[i] *= scaleFactor;
    }

    return convertZerosToNA(bins, numKmers);
}
