#pragma once

#include "numerics/global_optimization.hpp"

#include <utility>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/cbrt.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Barycentre;

template<typename Scalar, typename Argument>
MultiLevelSingleLinkage<Scalar, Argument>::MultiLevelSingleLinkage(
    Box const& box,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f)
    : box_(box),
      box_measure_(8 * Wedge(box_.vertices[0],
                             Wedge(box_.vertices[1], box_.vertices[2])).Norm()),
      f_(f),
      grad_f_(grad_f),
      random_(42),
      distribution_(-1.0, 1.0) {}

template<typename Scalar, typename Argument>
void MultiLevelSingleLinkage<Scalar, Argument>::FindGlobalMinimum(
    std::int64_t const values_per_round,
    std::int64_t const number_of_rounds) const {
  std::vector<Argument> arguments;
  for (std::int64_t k = 0; k < number_of_rounds; ++k) {
    std::vector<Argument> const arguments_k =
        GenerateArguments(box_, values_per_round);
  }
}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::GenerateArguments(
    Box const& box,
    std::int64_t const values_per_round) {
  std::vector<Argument> result;
  for (std::int64_t i = 0; i < values_per_round; ++i) {
    Argument argument = box.centre;
    for (const auto& vertex : box.vertices) {
      argument += distribution_(random_) * vertex;
    }
    result.push_back(argument);
  }
  return result;
}

template<typename Scalar, typename Argument>
typename Hilbert<Difference<Argument>>::NormType
MultiLevelSingleLinkage<Scalar, Argument>::rₖ(double const σ,
                                              std::int64_t const kN) {
  return Cbrt(3.0 * box_measure_ * std::log(kN) / (4.0 * π * kN));
}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
