#include "clr_benchmarks_adapter/quantities.hpp"

#include <vector>

#include "benchmarks/quantities.hpp"
#include "quantities/named_quantities.hpp"

using principia::quantities::Momentum;

namespace principia {
namespace clr_benchmarks_adapter {

void QuantitiesCLRBenchmark::DimensionfulDiscreteCosineTransform() {
  std::vector<Momentum> output;
  principia::benchmarks::DimensionfulDiscreteCosineTransform(&output);
}

void QuantitiesCLRBenchmark::DoubleDiscreteCosineTransform() {
  std::vector<double> output;
  principia::benchmarks::DoubleDiscreteCosineTransform(&output);
}

}  // namespace clr_benchmarks_adapter
}  // namespace principia
