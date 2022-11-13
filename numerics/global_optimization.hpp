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
using quantities::Cube;
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

//TODO(phl): Could this be a self-standing function?
template<typename Scalar, typename Argument>
class MultiLevelSingleLinkage {
 public:
  // A parallelepiped defined by its centre and the displacements of three
  // vertices.  Random points are uniformly distributed in the box.
  struct Box {
    Argument centre;
    std::array<Difference<Argument>, 3> vertices;
  };

  MultiLevelSingleLinkage(
      Box const& box,
      Field<Scalar, Argument> const& f,
      Field<Gradient<Scalar, Argument>, Argument> const& grad_f);

  void FindGlobalMinimum(std::int64_t values_per_round,
                         std::int64_t number_of_rounds) const;

 private:
  // Returns a vector of size |values_per_round|.
  std::vector<Argument> GenerateArguments(Box const& box,
                                          std::int64_t values_per_round);

  // Returns the radius rₖ from [] specialized for 3 dimensions.
  static typename Hilbert<Difference<Argument>>::NormType rₖ(double σ,
                                                             std::int64_t kN);

  Box const box_;
  Cube<typename Hilbert<Difference<Argument>>::NormType> const box_measure_;
  Field<Scalar, Argument> const f_;
  Field<Gradient<Scalar, Argument>, Argument> const grad_f_;

  std::mt19937_64 random_;
  std::uniform_real_distribution<> distribution_;
};

}  // namespace internal_global_optimization

using internal_global_optimization::MultiLevelSingleLinkage;

}  // namespace numerics
}  // namespace principia

#include "numerics/global_optimization_body.hpp"
