#pragma once

#include "numerics/nearest_neighbour.hpp"

#include <algorithm>
#include <numeric>
#include <type_traits>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_nearest_neighbour {

using geometry::Barycentre;
using geometry::Bivector;
using geometry::Frame;
using geometry::SymmetricBilinearForm;
using geometry::SymmetricProduct;
using geometry::Wedge;
using quantities::Difference;

template<typename T>
using SymmetricBilinearFormFor =
    decltype(SymmetricProduct(std::declval<T>(), std::declval<T>()));

template<typename T>
SymmetricBilinearFormFor<T> ComputePrincipalComponentForm(
    std::vector<T> const& ts) {
  SymmetricBilinearFormFor<T> result;
  for (const auto& t : ts) {
    result += SymmetricProduct(t, t);
  }
  return result;
}

template<typename Value_>
PrincipalComponentPartitioningTree<Value_>::PrincipalComponentPartitioningTree(
    std::vector<Value> const& values)
    : centroid_(Barycentre<Value, double, std::vector>(
          values,
          std::vector<double>(values.size(), 1))) {
  // The displacements from the centroid.
  std::vector<Difference<Value>> displacements;
  displacements.reserve(values.size());
  for (auto const& value : values) {
    displacements.push_back(value - centroid_);
  }

  // These are indices in |displacements|, that we'll use to refer to actual
  // values.  We cannot use iterators because they would be invalidated by
  // |Add|.  We use
  std::vector<std::int32_t> indices(displacements.size());
  std::iota(indices.begin(), indices.end(), 0);

  // Compute the "inertia" of the |vectors| and diagonalize it.  This gives us
  // the direction of the largest eigenvalue, which defines the splitting plane.
  using Eigenframe = Frame<enum class EigenframeTag>;

  auto const form = ComputePrincipalComponentForm(displacements);
  auto const eigensystem = form.Diagonalize<Eigenframe>();
  // The first eigenvalue is the largest one.
  //TODO(phl):Bivector?
  auto const principal_axis =
      eigensystem.rotation(Bivector<double, Eigenframe>({1, 0, 0}));

  // Project the |vectors| on the principal axis.
  using Projection =
      decltype(Wedge(principal_axis, std::declval<Difference<Value>>()));
  std::vector<Projection> projections;

  auto const projection_less = [&projections](std::int32_t const left,
                                              std::int32_t const right) {
    return projections[left].coordinates() < projections[right].coordinates();
  };

  projections.reserve(displacements.size());
  for (auto const& displacement : displacements) {
    projections.push_back(Wedge(principal_axis, displacement));
  }

  // Find the median of the projections on the principal axis.  Cost: 3.4 * N.
  auto const mid_upper = indices.begin() + indices.size() / 2;
  std::nth_element(indices.begin(), mid_upper, indices.end(), projection_less);

  // Here |mid_upper| denotes an index that denotes the point that is "in the
  // middle" of the projections on the principal axis.  We do not wish to anchor
  // our separator plane on a point, because that would make it more likely that
  // queries have to dig into two branches of the tree.  Instead, we look for
  // the maximum projection in the "lower" half space, and anchor our plane
  // halfway between these points.  Cost: 0.25 * N.
  auto const mid_lower =
      std::max_element(indices.begin(), mid_upper, projection_less);

  auto const anchor = Barycentre<Difference<Value>, double>(
      {displacements[*mid_lower], displacements[*mid_upper]},
      {1, 1});
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(Value const& value) {}

template<typename Value_>
Value_ const& PrincipalComponentPartitioningTree<Value_>::FindNearest(
    Value const& value) const {}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
