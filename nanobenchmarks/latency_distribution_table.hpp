#pragma once

#include <string>
#include <vector>

#include "absl/flags/declare.h"

ABSL_DECLARE_FLAG(std::vector<double>, quantiles);

namespace principia {
namespace nanobenchmarks {
namespace _latency_distribution_table {
namespace internal {

class LatencyDistributionTable {
 public:
  static std::vector<double> const& Quantiles();

  static std::string const& Heading();

  double min() const;

  void SetSamples(std::vector<double> const& samples);

  std::string Row() const;

 private:
  double min_ = 0.0;
  std::vector<double> measures_;

  friend LatencyDistributionTable operator*(double a,
                                            LatencyDistributionTable const& x);
  friend LatencyDistributionTable operator+(LatencyDistributionTable const& x,
                                            double b);
  friend LatencyDistributionTable operator-(LatencyDistributionTable const& x,
                                            double b);
};

LatencyDistributionTable operator*(double a, LatencyDistributionTable const& x);

LatencyDistributionTable operator+(LatencyDistributionTable const& x, double b);

LatencyDistributionTable operator-(LatencyDistributionTable const& x, double b);

}  // namespace internal

using internal::LatencyDistributionTable;

}  // namespace _latency_distribution_table
}  // namespace nanobenchmarks
}  // namespace principia
