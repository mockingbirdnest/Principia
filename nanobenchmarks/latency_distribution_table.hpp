#pragma once

#include <string>
#include <vector>

namespace principia {
namespace nanobenchmarks {
namespace _latency_distribution_table {
namespace internal {

class LatencyDistributionTable {
 public:
  explicit LatencyDistributionTable(std::vector<double> const& quantiles);

  double min() const;
  std::vector<double> const& quantiles() const;

  void SetSamples(std::vector<double> const& samples);

  std::string Heading();

  std::string Row() const;

 private:
  std::vector<double> const quantiles_;
  double min_;
  std::vector<double> measures_;

  friend LatencyDistributionTable operator*(double a,
                                            LatencyDistributionTable const& x);
  friend LatencyDistributionTable operator+(LatencyDistributionTable const& x,
                                            double b);
};

LatencyDistributionTable operator*(double a, LatencyDistributionTable const& x);

LatencyDistributionTable operator+(LatencyDistributionTable const& x, double b);

}  // namespace internal

using internal::LatencyDistributionTable;

}  // namespace _latency_distribution_table
}  // namespace nanobenchmarks
}  // namespace principia
