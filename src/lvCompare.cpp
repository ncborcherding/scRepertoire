#include <Rcpp.h>
#include <string>
#include <vector>

template <typename T>
T min(T a, T b, T c) {
    return std::min(a, std::min(b, c));
}

double editDist(const std::string& s1, const std::string& s2) {

    const int n = s2.size();
    const int m = s1.size();

    if (m == 0) return static_cast<double>(n);
    if (n == 0) return static_cast<double>(m);

    std::vector<int> prev(n + 1), curr(n + 1);

    for (int j = 0; j <= n; j++) {
        prev[j] = j;
    }

    for (int i = 1; i <= m; i++) {

        curr[0] = i;

        for (int j = 1; j <= n; j++) {
            curr[j] = min(
                curr[j - 1] + 1,
                prev[j] + 1,
                prev[j - 1] + ((s1[i - 1] == s2[j - 1]) ? 0 : 1)
            );
        }

        std::swap(prev, curr);
    }

    return static_cast<double>(prev[n]);
}

bool lenDiffWithinThreshold(const int len1, const int len2, const double threshold) {
    double lenDiff = static_cast<double>(std::abs(len1 - len2));
    double maxLen = static_cast<double>(std::max(len1, len2));
    return lenDiff <= (maxLen * (1 - threshold));
}

// [[Rcpp::export]]
Rcpp::DataFrame rcppGetSigSequenceEditDistEdgeListDf(
    const std::vector<std::string> sequences, const double threshold
) {

    std::vector<std::string> from, to;
    std::vector<double> distances;

    for (size_t i = 0; i < (sequences.size() - 1); i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {

            int len1 = sequences[i].size();
            int len2 = sequences[j].size();

            if (!lenDiffWithinThreshold(len1, len2, threshold)) {
                continue;
            }

            double meanLen = static_cast<double>(len1 + len2) * 0.5;
            double normalizedDistance = 1 - (editDist(sequences[i], sequences[j]) / meanLen);

            if (normalizedDistance < threshold) {
                continue;
            }

            from.push_back(sequences[i]);
            to.push_back(sequences[j]);
            distances.push_back(normalizedDistance);
        
        }
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("from") = from,
        Rcpp::Named("to") = to,
        Rcpp::Named("distance") = distances
    );

}
