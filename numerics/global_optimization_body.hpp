#pragma once

#include "numerics/global_optimization.hpp"

#include <utility>

#include "geometry/barycentre_calculator.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Barycentre;

template<typename Scalar, typename Argument>
MultiLevelSingleLinkage<Scalar, Argument>::MultiLevelSingleLinkage(
    Box const& box,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f,
    std::int64_t const values_per_round) {}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::GenerateArguments(
    Box const& box,
    std::int64_t const values_per_round,
    std::uniform_real_distribution<>& distribution) {
  std::vector<Argument> result;
  for (std::int64_t i = 0; i < values_per_round; ++i) {
    Argument argument = box.centre;
    for (const auto& vertex : box.vertices) {
      argument += distribution() * vertex;
    }
    result.push_back(argument);
  }
  return result;
}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
