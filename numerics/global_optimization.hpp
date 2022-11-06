#pragma once

#include <functional>
#include <random>
#include <vector>

#include "geometry/hilbert.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Hilbert;
using quantities::Difference;
using quantities::Length;
using quantities::Product;
using quantities::Quotient;

// In this file |Argument| must be such that its difference belongs to a Hilbert
// space.

template<typename Scalar, typename Argument>
using Field = std::function<Scalar(Argument const&)>;

template<typename Scalar, typename Argument>
using Gradient =
    Product<Scalar,
            Quotient<Difference<Argument>,
                     typename Hilbert<Difference<Argument>>::Norm²Type>>;

template<typename Scalar, typename Argument>
class MultiLevelSingleLinkage {
 public:
  // A parallelepiped defined by its centre and the displacements of three
  // vertices.
  struct Box {
    Argument centre;
    std::array<Difference<Argument>, 3> vertices;
  };

  MultiLevelSingleLinkage(
      Box const& box,
      Field<Scalar, Argument> const& f,
      Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
      std::int64_t values_per_round);

 private:
  // The distribution must have bounds [-1, 1].  Returns a vector of size
  // |values_per_round|.
  static std::vector<Argument> GenerateArguments(
      Box const& box,
      std::int64_t values_per_round,
      std::uniform_real_distribution<>& distribution);
};

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia

#include "numerics/global_optimization_body.hpp"
