#include "latency_distribution_table.hpp"

#include <print>
#include <sstream>

#include "absl/flags/flag.h"

namespace principia {
namespace nanobenchmarks {
namespace _latency_distribution_table {
namespace internal {

LatencyDistributionTable::LatencyDistributionTable(
    std::vector<double> const& quantiles)
    : quantiles_(quantiles) {}

double LatencyDistributionTable::min() const {
  return min_;
}

std::vector<double> const& LatencyDistributionTable::quantiles() const {
  return quantiles_;
}

void LatencyDistributionTable::SetSamples(std::vector<double> const& samples) {
  min_ = samples[0];
  for (double const q : quantiles_) {
    measures_.push_back(samples[(samples.size() - 1) * q]);
  }
}

std::string LatencyDistributionTable::Heading() {
  std::stringstream& out = *new std::stringstream;
  std::print(out, "{:>8}", "min");
  for (double const q : quantiles_) {
    if (q < 1e-3) {
      std::print(out, "{:>7}‱", 10'000 * q);
    } else if (q < 1e-2) {
      std::print(out, "{:>7}‰", 1000 * q);
    } else {
      std::print(out, "{:>7}%", 100 * q);
    }
  }
  return out.str();
}

inline std::string LatencyDistributionTable::Row() const {
  std::stringstream out;
  std::print(out, "{:8.2f}", min_);
  for (double const measure : measures_) {
    std::print(out, "{:+8.2f}", measure - min_);
  }
  return out.str();
}

LatencyDistributionTable operator*(double const a,
                                   LatencyDistributionTable const& x) {
  LatencyDistributionTable result(x.quantiles());
  result.min_ = a * x.min_;
  for (double const measure : x.measures_) {
    result.measures_.push_back(a * measure);
  }
  return result;
}

LatencyDistributionTable operator+(LatencyDistributionTable const& x,
                                   double const b) {
  LatencyDistributionTable result(x.quantiles());
  result.min_ = x.min_ + b;
  for (double const measure : x.measures_) {
    result.measures_.push_back(measure + b);
  }
  return result;
}

}  // namespace internal
}  // namespace _latency_distribution_table
}  // namespace nanobenchmarks
}  // namespace principia
