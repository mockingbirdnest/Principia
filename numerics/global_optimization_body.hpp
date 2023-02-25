#pragma once

#include "numerics/global_optimization.hpp"

#include <algorithm>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/macros.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "numerics/gradient_descent.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using base::noreturn;
using geometry::Barycentre;
using geometry::Wedge;
using quantities::Cbrt;
using quantities::Infinity;
using quantities::Pow;
using namespace principia::base::_not_null;

// Values that are too small cause extra activity in |Add| and deeper recursion
// in |FindNearestNeighbour|.  Values that are too large cause extra function
// evaluations in the leaf search.
constexpr int64_t pcp_tree_max_values_per_cell = 20;

template<typename Scalar, typename Argument, int dimensions>
typename MultiLevelSingleLinkage<Scalar, Argument, dimensions>::Box::Measure
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::Box::measure() const {
  if constexpr (dimensions == 1) {
    return 2 * vertices[0].Norm();
  } else if constexpr (dimensions == 2) {
    return 4 * Wedge(vertices[0], vertices[1]).Norm();
  } else if constexpr (dimensions == 3) {
    return 8 * Wedge(vertices[0], Wedge(vertices[1], vertices[2])).Norm();
  }
  noreturn();
}

template<typename Scalar, typename Argument, int dimensions>
typename MultiLevelSingleLinkage<Scalar, Argument, dimensions>::NormType
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::Box::diametre() const {
  NormType diameter;
  for (auto const& vertex : vertices) {
    diameter = std::max(diameter, 2 * vertex.Norm());
  }
  return diameter;
}

template<typename Scalar, typename Argument, int dimensions>
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::MultiLevelSingleLinkage(
    Box const& box,
    Field<Scalar, Argument> f,
    Field<Gradient<Scalar, Argument>, Argument> grad_f)
    : box_(box),
      box_diametre_(box_.diametre()),
      box_measure_(box_.measure()),
      f_(std::move(f)),
      grad_f_(std::move(grad_f)),
      random_(42),
      distribution_(-1.0, 1.0) {
  for (auto const& vertex : box_.vertices) {
    CHECK_NE(vertex, Difference<Argument>{});
  }
}

template<typename Scalar, typename Argument, int dimensions>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::FindGlobalMaxima(
    std::int64_t const points_per_round,
    std::optional<std::int64_t> const number_of_rounds,
    NormType const local_search_tolerance) {
  Field<Scalar, Argument> minus_f = [this](Argument const& x) {
    return -f_(x);
  };
  Field<Gradient<Scalar, Argument>, Argument> minus_grad_f =
      [this](Argument const& x) {
    return -grad_f_(x);
  };

  return FindGlobalMinima(points_per_round,
                          number_of_rounds,
                          local_search_tolerance,
                          minus_f,
                          minus_grad_f);
}

template<typename Scalar, typename Argument, int dimensions>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::FindGlobalMinima(
    std::int64_t const points_per_round,
    std::optional<std::int64_t> const number_of_rounds,
    NormType const local_search_tolerance) {
  return FindGlobalMinima(points_per_round,
                          number_of_rounds,
                          local_search_tolerance,
                          f_,
                          grad_f_);
}

template<typename Scalar, typename Argument, int dimensions>
std::vector<Argument>
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::FindGlobalMinima(
    std::int64_t const points_per_round,
    std::optional<std::int64_t> const number_of_rounds,
    NormType const local_search_tolerance,
    Field<Scalar, Argument> const& f,
    Field<Gradient<Scalar, Argument>, Argument> const& grad_f) {
  // This is the set X* from [RT87b].
  Arguments stationary_points;

  // The sample points.  We'll have at least |points_per_round| points.
  Arguments points;
  points.reserve(points_per_round);

  // The PCP tree used for detecting proximity of the stationary points.  It
  // gets updated as new stationary points are found.
  PrincipalComponentPartitioningTree<Argument> stationary_point_neighbourhoods(
      /*values=*/{},
      pcp_tree_max_values_per_cell);

  // The PCP tree used for detecting proximity of the sample points .  It gets
  // updated as new points are generated.
  PrincipalComponentPartitioningTree<Argument> point_neighbourhoods(
      /*values=*/{},
      pcp_tree_max_values_per_cell);

  int number_of_local_searches = 0;

  // This structure corresponds to the list T in [RT87b].  Points are ordered
  // based on their distance to their nearest neighbour that has a lower value
  // of |f|.  When new points are added to the map, they are initially given the
  // key |Infinity|, which is then refined when the nearest neighbour search
  // happens.  Each time through the outer loop, only the points in the upper
  // half of the map (distance greater than rₖ) are considered, and they are
  // moved down if a "too close" neighbour is found.
  // We modify this map while iterating so we need iterator stability.
  std::multimap<Norm²Type, Argument const*> schedule;

  // There are two ways to iterate through this loop, corresponding to different
  // termination conditions.  We can either do a fixed number of rounds with a
  // fixed number of points each, as described in the bulk of [RT87a] and
  // [RT87b]; or we can use the optimal Bayesian stopping rule described in
  // proposition 2 of [RT87a] and equations (3) and (4) of [RT87b].  For
  // consistency with most of the description in [RT87a] and [RT87b], we call
  // |N| the number of points added in a particular iteration and |kN| the total
  // number of sampling points so far, even if |N| varies from one iteration to
  // the next and |kN| is not necessarily |k * N|.
  std::int64_t number_of_points_all_rounds;
  std::int64_t number_of_points_this_round;
  std::int64_t round = 1;
  auto& kN = number_of_points_all_rounds;
  auto& N = number_of_points_this_round;
  for (auto& k = round;; ++k) {
    if (k == 1) {
      N = points_per_round;
      kN = N;
    } else if (number_of_rounds.has_value()) {
      if (k > number_of_rounds.value()) {
        break;
      } else {
        N = points_per_round;
        kN += N;
      }
    } else {
      std::int64_t const w = stationary_points.size();
      // Note that |w| may be 0 here if the stationary points were all outside
      // the local search radius.
      DCHECK_GT(kN, w + 2);
      // [RT87b] equation 3.
      if (kN > w + 2 &&
          static_cast<double>(w * (kN - 1)) / static_cast<double>(kN - w - 2) <=
          w + 0.5) {
        break;
      } else {
        // [RT87b] equation 4.
        std::int64_t const M = (2 * w + 3) * w + 2;
        N = M - kN;
        DCHECK_LT(0, N);
        kN = M;
      }
      DLOG(ERROR) << "Round #" << k << " with " << N << " points";
    }

    // Generate N new random points and add them to the PCP tree and to the
    // |schedule| map.  Note that while [RT87b] tells us in the description of
    // algorithm A that "it is no longer necessary to reduce the sample" they
    // also tell us in section 4 that they reduce the sample using γ = 0.1.
    // Also, [KS05] reduce the sample, but they don't seem to agree on whether
    // rₖ depends on γ.
    // Anyway, reducing the sample would be annoying with our data structures,
    // so let's not go there, 'tis a silly place.
    Arguments pointsₖ = RandomArguments(N);
    for (auto& pointₖ : pointsₖ) {
      points.push_back(std::move(pointₖ));
      Argument const* const pointₖ_pointer = points.back().get();
      point_neighbourhoods.Add(pointₖ_pointer);
      schedule.emplace_hint(
          schedule.end(), Infinity<Norm²Type>, pointₖ_pointer);
    }

    // Compute the radius below which we won't do a local search in this
    // iteration.
    Norm²Type const rₖ² = CriticalRadius²(/*σ=*/4, kN);

    // Process the points whose nearest neighbour is "sufficiently far" (or
    // unknown).
    for (auto it = schedule.upper_bound(rₖ²); it != schedule.end();) {
      Argument const& xᵢ = *it->second;
      auto* const xⱼ = point_neighbourhoods.FindNearestNeighbour(
          xᵢ, [f_xᵢ = f(xᵢ), rₖ², f, xᵢ](Argument const* const xⱼ) {
            return (xᵢ - *xⱼ).Norm²() <= rₖ² && f(*xⱼ) < f_xᵢ;
          });

      if (xⱼ == nullptr) {
        // We must do a local search as xᵢ couldn't be added to an existing
        // cluster.  Note that the radius of the search has to be the diametre
        // of the box: it's possible that xᵢ would be near one vertex of the box
        // and the stationary point near the opposite vertex.
        ++number_of_local_searches;
        auto const stationary_point =
            BroydenFletcherGoldfarbShanno(xᵢ,
                                          f,
                                          grad_f,
                                          local_search_tolerance,
                                          box_diametre_);

        // If the new stationary point is sufficiently far from the ones we
        // already know, record it.
        if (stationary_point.has_value() &&
            IsNewStationaryPoint(stationary_point.value(),
                                 stationary_point_neighbourhoods,
                                 local_search_tolerance)) {
          stationary_points.push_back(
              std::make_unique<Argument>(stationary_point.value()));
          stationary_point_neighbourhoods.Add(stationary_points.back().get());
        }
        // A local search has been started from xᵢ, so no point in considering
        // it again.
        it = schedule.erase(it);
      } else {
        // Move the point xᵢ "down" in the |schedule| map, based on the distance
        // to its nearest neighbour.
        auto const distance²_to_xⱼ = (xᵢ - *xⱼ).Norm²();
        DCHECK_LE(distance²_to_xⱼ, rₖ²);
        it = schedule.erase(it);
        // This insertion take places below |schedule.upper_bound(rₖ²)|, so it
        // doesn't affect the current iteration.
        schedule.emplace(distance²_to_xⱼ, &xᵢ);
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

template<typename Scalar, typename Argument, int dimensions>
bool MultiLevelSingleLinkage<Scalar, Argument, dimensions>::
IsNewStationaryPoint(Argument const& stationary_point,
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

template<typename Scalar, typename Argument, int dimensions>
std::vector<not_null<std::unique_ptr<Argument>>>
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::RandomArguments(
    std::int64_t const values_per_round) {
  Arguments arguments;
  for (std::int64_t i = 0; i < values_per_round; ++i) {
    auto argument = make_not_null_unique<Argument>(box_.centre);
    for (const auto& vertex : box_.vertices) {
      *argument += distribution_(random_) * vertex;
    }
    arguments.push_back(std::move(argument));
  }
  return arguments;
}

template<typename Scalar, typename Argument, int dimensions>
typename MultiLevelSingleLinkage<Scalar, Argument, dimensions>::Norm²Type
MultiLevelSingleLinkage<Scalar, Argument, dimensions>::CriticalRadius²(
    double const σ,
    std::int64_t const kN) {
  if constexpr (dimensions == 1) {
    return Pow<2>(box_measure_ * σ * std::log(kN) / (2.0 * kN));
  } else if constexpr (dimensions == 2) {
    return box_measure_ * σ * std::log(kN) / (π * kN);
  } else if constexpr (dimensions == 3) {
    return Cbrt(Pow<2>(3.0 * box_measure_ * σ * std::log(kN) / (4.0 * π * kN)));
  }
  noreturn();
}

}  // namespace internal_global_optimization
}  // namespace numerics
}  // namespace principia
