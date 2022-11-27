#pragma once

#include <functional>
#include <optional>
#include <random>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/nearest_neighbour.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_global_optimization {

using base::not_null;
using geometry::Hilbert;
using quantities::Difference;
using quantities::Exponentiation;
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

// NOTE(phl): This could nearly be a self-standing function (it doesn't have
// much state) but having the type |Box| floating around would be unpleasant.
// Plus, that would be too many parameters in that function.
// The parameter |dimensions| may be used to specify the dimension of the
// problem.  It it is 1 or 2, the box is 1- or 2-dimensional and the computation
// of rₖ is adjusted accordingly.  In all cases, the dimensions of the box must
// be nonzero.
template<typename Scalar, typename Argument, int dimensions = 3>
class MultiLevelSingleLinkage {
 public:
  using NormType = typename Hilbert<Difference<Argument>>::NormType;

  // A parallelepiped defined by its centre and the displacements of three
  // vertices.  Random points are uniformly distributed in the box.
  struct Box {
    Argument centre;
    std::array<Difference<Argument>, dimensions> vertices;

    using Measure = Exponentiation<NormType, dimensions>;
    Measure measure() const;

    bool contains(Argument const& point) const;
  };

  MultiLevelSingleLinkage(
      Box const& box,
      Field<Scalar, Argument> const& f,
      Field<Gradient<Scalar, Argument>, Argument> const& grad_f);

  // If |number_of_rounds| is given, the algorithm does |number_of_rounds|
  // iterations, each time adding |points_per_round| to the sample.
  // If |number_of_rounds| is omitted, the first iteration uses
  // |points_per_round| points, and subsequent iterations adjust the number of
  // points (or the decision to terminate) based on the optimal Bayesian
  // stopping rule.
  // Beware!  The Bayesian stopping rule is typically more efficient, but it is
  // only technically correct (the best kind of correct) if the relative sizes
  // of the regions of attraction follow a uniform distribution.
  std::vector<Argument> FindGlobalMinima(
      std::int64_t points_per_round,
      std::optional<std::int64_t> number_of_rounds,
      NormType local_search_tolerance);

 private:
  using Norm²Type = typename Hilbert<Difference<Argument>>::Norm²Type;

  // We need pointer stability for the arguments as we store pointers, e.g., in
  // PCP trees.  We generally cannot |reserve| because we don't know the final
  // size of the vector, hence the |unique_ptr|s.
  using Arguments = std::vector<not_null<std::unique_ptr<Argument>>>;

  // Returns true iff the given |stationary_point| is sufficiently far from the
  // ones already in |stationary_point_neighbourhoods|.
  static bool IsNewStationaryPoint(
      Argument const& stationary_point,
      PrincipalComponentPartitioningTree<Argument> const&
          stationary_point_neighbourhoods,
      NormType local_search_tolerance);

  // Returns a vector of size |values_per_round|.  The points are in |box_|.
  Arguments RandomArguments(std::int64_t values_per_round);

  // Returns the square of the radius rₖ from [RT87a], eqn. 35, specialized for
  // |dimensions|.
  Norm²Type CriticalRadius²(double σ, std::int64_t kN);

  Box const box_;
  Box::Measure const box_measure_;
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
