#pragma once

namespace principia {
namespace benchmarks {
namespace _metric {
namespace internal {

enum class Metric {
  Latency,
  Throughput
};

}  // namespace internal

using internal::Metric;

}  // namespace _metric
}  // namespace benchmarks
}  // namespace principia
