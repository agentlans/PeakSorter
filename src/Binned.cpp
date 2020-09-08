#include <Rcpp.h>

#include <queue>
#include <unordered_map>
#include <vector>
#include <set>
#include <algorithm>

using namespace Rcpp;
using std::vector;

// [[Rcpp::plugins(cpp11)]]

// LC-MS peak
struct Peak {
  int peak_id;
  int bin_id;
  double intensity;
  
  bool operator<(const Peak& peak) const {
    return intensity < peak.intensity;
  }
};

/* Returns the IDs of peaks (0, 1, 2, ...) ordered according to the binned method.
 bin_id[i] and intensity[i] represent the bin assignment and the intensity of peak[i] respectively. */
vector<int> binned_method_cpp(const vector<int>& bin_id, const vector<double>& intensity) {
  // Assume bin_id.size() == intensity.size()
  int n = bin_id.size();
  // Load the peaks into data structure
  std::unordered_map<int, std::priority_queue<Peak>> bins;
  std::set<int> distinct_bins;
  for (int i = 0; i < n; ++i) {
    // The bin of peak i
    int ib = bin_id[i];
    Peak p = {i, ib, intensity[i]};
    // Push peak into bin
    bins[ib].push(p);
    // Note that this bin exists
    distinct_bins.insert(ib);
  }
  // Take the peaks out of data structure in order
  vector<int> sorted_peaks;
  sorted_peaks.reserve(n);
  while (sorted_peaks.size() < n) {
    vector<Peak> top_peaks;
    top_peaks.reserve(distinct_bins.size());
    
    // For each bin, choose the top peak
    for (int bin : distinct_bins) {
      if (!bins[bin].empty()) {
        Peak top_peak = bins[bin].top();
        top_peaks.push_back(top_peak);
        bins[bin].pop();
      }
    }
    // The top peaks from highest to lowest intensity
    std::stable_sort(top_peaks.begin(), top_peaks.end(), 
              [](const Peak& p1, const Peak& p2) {
                return !(p1 < p2);
              });
    // Top peaks from highest to lowest intensity
    for (const auto& p : top_peaks) {
      sorted_peaks.push_back(p.peak_id);
    }
  }
  return sorted_peaks;
}

// Wrapper for the above function using Rcpp data structures

//' Internal function called as part of the binned algorithm. Not meant for the
//' end user.
//' 
//' In the following, i is 0, 1, 2, ...
//' @param bin_id A vector where bin_id[i] represents the bin number of peak i.
//' @param intensity A vector where intensity[i] represents the intensity of peak i.
//' @return Vector of i indicating the order of peaks as sorted by the binned method
// [[Rcpp::export]]
NumericVector binned_method_cpp_wrap(NumericVector bin_id, NumericVector intensity) {
  auto bi = as<vector<int>>(bin_id);
  auto inten = as<vector<double>>(intensity);
  return wrap(binned_method_cpp(bi, inten));
}

