#include <Rcpp.h>
#include <string>
#include <vector>

double editDist(const std::string& s1, const std::string& s2) {

    const auto n = s2.size();

    std::vector<int64_t> v0;
    v0.resize(n + 1);
    std::iota(v0.begin(), v0.end(), 0);

    auto v1 = v0;

    auto s1char = s1[0];
    for (size_t j = 0; j < n; j++) {
        auto delCost = v0[j + 1] + 1;
        auto insCost = v1[j] + 1;
        int substCost = s1char != s2[j];
        v1[j + 1] = std::min({v0[j] + substCost, delCost, insCost});
    }

    return static_cast<double>(v0[n]);
}

// [[Rcpp::export]]
Rcpp::DataFrame rcppGetSigSequenceEditDistEdgeListDf(
    std::vector<std::string> sequences, double threshold
) {

    std::vector<std::string> from, to;
    std::vector<double> distances;

    for (size_t i = 0; i < sequences.size() - 1; i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {

            int len1 = sequences[i].size();
            int len2 = sequences[j].size();

            if (std::abs(len1 - len2 > std::max(len1, len2) * (1 - threshold))) {
                continue;
            }

            double distance = editDist(sequences[i], sequences[j]);
            double normalizedDistance = 1 - (2 * distance / (len1 + len2));

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
