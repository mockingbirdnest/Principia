#pragma once

#include "numerics/global_optimization.hpp"

#include <utility>
#include <vector>

#include "geometry/barycentre_calculator.hpp"
#include "numerics/gradient_descent.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using geometry::Barycentre;
using geometry::Wedge;
using quantities::Cbrt;
using quantities::Infinity;

// TODO(phl): Provide a way to parameterize the PCP trees?
constexpr int64_t pcp_tree_max_values_per_cell = 10;

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
      distribution_(-1.0, 1.0) {
  CHECK_LT(Cube<typename Hilbert<Difference<Argument>>::NormType>{},
           box_measure_);
}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::FindGlobalMinima(
    std::int64_t const points_per_round,
    std::int64_t const number_of_rounds,
    NormType const local_search_tolerance) {
  const std::int64_t N = points_per_round;

  // This is the set X* from [RT87b].  It contains |unique_ptr|s because we need
  // pointer stability for the PCP tree but we cannot precompute the vector
  // size.
  std::vector<std::unique_ptr<Argument>> stationary_points;

  // The PCP tree used for detecting proximity of the stationary points.  It
  // gets updated as new stationary points are found.
  PrincipalComponentPartitioningTree<Argument> stationary_point_neighbourhoods(
      /*values=*/{},
      pcp_tree_max_values_per_cell);

  // The sample points.  Make sure that pointers don't get invalidated as new
  // points are added.
  std::vector<Argument> points;
  points.reserve(N * number_of_rounds);

  // The PCP tree used for detecting proximity of the sample points .  It gets
  // updated as new points are generated.
  PrincipalComponentPartitioningTree<Argument> point_neighbourhoods(
      /*values=*/{},
      pcp_tree_max_values_per_cell);

  int number_of_local_searches = 0;

  // This structure corresponds to the list T in [RT87b].  Points are ordered
  // based on their distance to their nearest neighbour that has a lower value
  // of |f_|.  When new points are added to the map, they are initially given
  // the key |Infinity|, which is then refined when the nearest neighbour search
  // happens.  Each time through the outer loop, only the points in the upper
  // half of the map (distance greater than rₖ) are considered, and they are
  // moved down if a "too close" neighbour is found.
  std::multimap<NormType, Argument const*> schedule;

  for (std::int64_t k = 1; k <= number_of_rounds; ++k) {
    // Generate N new random points and add them to the PCP tree and to the
    // |schedule| map.  Note that while [RT87b] tells us in the description of
    // algorithm A that "it is no longer necessary to reduce the sample" they
    // also tell us in section 4 that they reduce the sample using γ = 0.1.
    // Also, [KS05] reduce the sample, but they don't seem to agree on whether
    // rₖ depends on γ.
    // Anyway, reducing the sample would be annoying with our data structures,
    // so let's go there, 'tis a silly place.
    std::vector<Argument> pointsₖ = RandomArguments(N);
    for (auto& pointₖ : pointsₖ) {
      points.push_back(std::move(pointₖ));
      Argument const* const pointₖ_pointer = &points.back();
      point_neighbourhoods.Add(pointₖ_pointer);
      schedule.emplace_hint(schedule.end(), Infinity<NormType>, pointₖ_pointer);
    }

    // Compute the radius below which we won't do a local search in this
    // iteration.
    auto const rₖ = CriticalRadius(/*σ=*/4, /*kN=*/k * N);

    // Process the points whose nearest neighbour is "sufficiently far" (or
    // unknown).
    for (auto it = schedule.upper_bound(rₖ); it != schedule.end();) {
      Argument const& xᵢ = *it->second;
      auto* const xⱼ = point_neighbourhoods.FindNearestNeighbour(
          xᵢ, [this, f_xᵢ = f_(xᵢ), rₖ, xᵢ](Argument const* const xⱼ) {
            return (xᵢ - *xⱼ).Norm() <= rₖ && f_(*xⱼ) < f_xᵢ;
          });

      if (xⱼ == nullptr) {
        // We must do a local search as the point xᵢ couldn't be added to an
        // existing cluster.
        ++number_of_local_searches;
        auto const stationary_point = BroydenFletcherGoldfarbShanno(
            xᵢ, f_, grad_f_, local_search_tolerance);

        // If the new stationary point is sufficiently far from the ones we
        // already know, record it.
        if (IsNewStationaryPoint(stationary_point,
                                 stationary_point_neighbourhoods,
                                 local_search_tolerance)) {
          stationary_points.push_back(
              std::make_unique<Argument>(stationary_point));
          stationary_point_neighbourhoods.Add(stationary_points.back().get());
        }
        // A local search has been started from xᵢ, so no point in considering
        // it again.
        it = schedule.erase(it);
      } else {
        // Move the point xᵢ "down" in the |schedule| map, based on the distance
        // to its nearest neighbour.
        auto const distance_to_xⱼ = (xᵢ - *xⱼ).Norm();
        DCHECK_LE(distance_to_xⱼ, rₖ);
        it = schedule.erase(it);
        // This insertion take places below |schedule.upper_bound(rₖ)|, so it
        // doesn't affect the current iteration.
        schedule.emplace(distance_to_xⱼ, &xᵢ);
      }
    }
  }

  DLOG(ERROR) << "Number of local searches: " << number_of_local_searches;

  std::vector<Argument> result;
  result.reserve(stationary_points.size());
  for (const auto& stationary_point : stationary_points) {
    result.push_back(*stationary_point);
  }

  return result;
}

template<typename Scalar, typename Argument>
bool MultiLevelSingleLinkage<Scalar, Argument>::IsNewStationaryPoint(
    Argument const& stationary_point,
    PrincipalComponentPartitioningTree<Argument> const&
        stationary_point_neighbourhoods,
    NormType const local_search_tolerance) {
  auto* const stationary_point_nearest_neighbour =
      stationary_point_neighbourhoods.FindNearestNeighbour(stationary_point);

  // The factor 2 below is because we are looking at the distance between two
  // approximate stationary points, not between an approximate stationary point
  // and the real one.
  return stationary_point_nearest_neighbour == nullptr ||
         (stationary_point - *stationary_point_nearest_neighbour).Norm() >
             2 * local_search_tolerance;
}

template<typename Scalar, typename Argument>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument>::RandomArguments(
    std::int64_t const values_per_round) {
  std::vector<Argument> result;
  for (std::int64_t i = 0; i < values_per_round; ++i) {
    Argument argument = box_.centre;
    for (const auto& vertex : box_.vertices) {
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
  return Cbrt(3.0 * box_measure_ * σ * std::log(kN) / (4.0 * π * kN));
}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
