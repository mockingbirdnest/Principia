#pragma once

#include "numerics/nearest_neighbour.hpp"

#include <algorithm>
#include <numeric>
#include <type_traits>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace numerics {
namespace internal_nearest_neighbour {

using base::make_not_null_unique;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::SymmetricBilinearForm;
using geometry::Wedge;


template<typename Value_>
PrincipalComponentPartitioningTree<Value_>::PrincipalComponentPartitioningTree(
    std::vector<Value> const& values)
    : centroid_(Barycentre<Value, double, std::vector>(
          values,
          std::vector<double>(values.size(), 1))) {
  displacements_.reserve(values.size());
  for (auto const& value : values) {
    displacements_.push_back(value - centroid_);
  }

  Indices indices(values.size());
  std::iota(indices.begin(), indices.end(), 0);  // Yes!  I did it!

  root_ = BuildTree(indices.begin(), indices.end(), indices.size());
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(Value const& value) {}

template<typename Value_>
Value_ const& PrincipalComponentPartitioningTree<Value_>::FindNearest(
    Value const& value) const {}

template<typename Value_>
not_null<std::unique_ptr<
typename PrincipalComponentPartitioningTree<Value_>::Node>>
PrincipalComponentPartitioningTree<Value_>::BuildTree(
    Indices::iterator const begin,
    Indices::iterator const end,
    std::int64_t const size) const {
  // Compute the "inertia" of the selected displacements and diagonalize it.
  // This gives us the direction of the largest eigenvalue, which defines the splitting plane.
  auto const form = ComputePrincipalComponentForm(begin, end);
  auto const eigensystem =
      form.template Diagonalize<PrincipalComponentsFrame>();
  // The first eigenvalue is the largest one.
  auto const principal_axis = eigensystem.rotation(
      Bivector<double, PrincipalComponentsFrame>({1, 0, 0}));

  // Project the |vectors| on the principal axis.
  using Projection =
      decltype(Wedge(principal_axis, std::declval<Displacement>()));
  std::vector<Projection> projections;

  auto const projection_less = [&projections](std::int32_t const left,
                                              std::int32_t const right) {
    return projections[left].coordinates() < projections[right].coordinates();
  };

  projections.reserve(size);
  for (auto it = begin; it != end; ++it) {
    auto const& displacement = displacements_[*it];
    projections.push_back(Wedge(principal_axis, displacement));
  }

  // Find the median of the projections on the principal axis.  Cost: 3.4 * N.
  auto const mid_upper = begin + size / 2;
  std::nth_element(begin, mid_upper, end, projection_less);

  // Here |mid_upper| denotes an index that denotes the point that is "in the
  // middle" of the projections on the principal axis.  We do not wish to anchor
  // our separator plane on a point, because that would make it more likely that
  // queries have to dig into two branches of the tree.  Instead, we look for
  // the maximum projection in the "lower" half space, and anchor our plane
  // halfway between these points.  Cost: 0.25 * N.
  auto const mid_lower = std::max_element(begin, mid_upper, projection_less);

  auto const anchor = Barycentre<Displacement, double>(
      {displacements_[*mid_lower], displacements_[*mid_upper]}, {1, 1});

  auto first_child = BuildTree(begin, mid_upper, size / 2);
  auto second_child = BuildTree(mid_upper, end, size - size / 2);

  return make_not_null_unique<Node>(
      Node{.principal_axis = principal_axis,
           .anchor = anchor,
           .tree = Children(std::move(first_child), std::move(second_child))});
}

template<typename Value_>
typename PrincipalComponentPartitioningTree<Value_>::
DisplacementSymmetricBilinearForm
PrincipalComponentPartitioningTree<Value_>::ComputePrincipalComponentForm(
    Indices::iterator const begin,
    Indices::iterator const end) const {
  DisplacementSymmetricBilinearForm result;
  for (auto it = begin; it != end; ++it) {
    auto const& displacement = displacements_[*it];
    result += SymmetricProduct(displacement, displacement);
  }
  return result;
}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
