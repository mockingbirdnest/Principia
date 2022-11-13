#pragma once

#include "numerics/global_optimization.hpp"

#include <utility>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/gradient_descent.hpp"
#include "numerics/nearest_neighbour.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Barycentre;
using geometry::Wedge;
using quantities::Cbrt;

template<typename Scalar, typename Argument>
Cube<typename Hilbert<Difference<Argument>>::NormType>
MultiLevelSingleLinkage<Scalar, Argument>::Box::Measure() const {
  return 8 * Wedge(vertices[0], Wedge(vertices[1], vertices[2])).Norm();
}

template<typename Scalar, typename Argument>
MultiLevelSingleLinkage<Scalar, Argument>::MultiLevelSingleLinkage(
    Box const& box,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f)
    : box_(box),
      box_measure_(box_.Measure()),
      f_(f),
      grad_f_(grad_f),
      random_(42),
      distribution_(-1.0, 1.0) {}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::FindGlobalMinima(
    std::int64_t const values_per_round,
    std::int64_t const number_of_rounds,
    NormType const local_search_tolerance) {
  const std::int64_t N = values_per_round;

  std::vector<Argument> stationary_points;

  // The PCP tree used for nearest neighbour computation.  It gets updated as
  // new points are generated.
  //TODO(phl):parameters
  PrincipalComponentPartitioningTree<Argument> pcp_tree(
      /*values=*/{},
      /*max_values_per_cell=*/10);

  // Make sure that pointers don't get invalidated as new arguments are added.
  std::vector<Argument> arguments;
  arguments.reserve(N * number_of_rounds);

  // TODO(phl): This is quadratic.  Make the algorithm linear once we believe
  // that it is correct.
  for (std::int64_t k = 0; k < number_of_rounds; ++k) {
    auto const rₖ = CriticalRadius(/*σ=*/4, /*kN=*/k * N);

    // Generate N new random points and add them to the PCP tree.
    std::vector<Argument> argumentsₖ = RandomArguments(box_, N);
    for (auto& argumentₖ : argumentsₖ) {
      arguments.push_back(std::move(argumentₖ));
      pcp_tree.Add(&arguments.back());
    }

    for (auto& xᵢ : arguments) {
      auto* const xⱼ = pcp_tree.FindNearestNeighbour(
          xᵢ, [this, f_xᵢ = f_(xᵢ)](Argument const* const xⱼ) {
            return f_(*xⱼ) < f_xᵢ;
          });
      if (xⱼ == nullptr || (xᵢ - *xⱼ).Norm() > rₖ) {
        stationary_points.push_back(BroydenFletcherGoldfarbShanno(
            xᵢ, f_, grad_f_, local_search_tolerance));
      }
    }
  }

  return stationary_points;
}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::RandomArguments(
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
MultiLevelSingleLinkage<Scalar, Argument>::CriticalRadius(
    double const σ,
    std::int64_t const kN) {
  return Cbrt(3.0 * box_measure_ * std::log(kN) / (4.0 * π * kN));
}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
