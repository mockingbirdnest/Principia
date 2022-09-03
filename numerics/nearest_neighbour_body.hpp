#pragma once

#include "numerics/nearest_neighbour.hpp"

#include <algorithm>
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
    std::vector<Value> const& values,
    std::int64_t const max_values_per_cell)
    : max_values_per_cell_(max_values_per_cell),
      centroid_(Barycentre<Value, double, std::vector>(
          values,
          std::vector<double>(values.size(), 1))) {
  displacements_.reserve(values.size());
  for (auto const& value : values) {
    displacements_.push_back(value - centroid_);
  }

  Indices indices;
  indices.reserve(displacements_.size());
  for (int i = 0; i < displacements_.size(); ++i) {
    indices[i] = {.index = i, .projection = Norm{}};
  }

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
  if (size <= max_values_per_cell_) {
    // We are done subdividing, return a leaf.
    Leaf leaf;
    for (auto it = begin; it != end; ++it) {
      leaf.push_back(it->index);
    }
    return make_not_null_unique<Node>(leaf);
  }

  // Compute the "inertia" of the selected displacements and diagonalize it.
  // This gives us the direction of the largest eigenvalue, which defines the
  // splitting plane.
  auto const form = ComputePrincipalComponentForm(begin, end);
  auto const eigensystem =
      form.template Diagonalize<PrincipalComponentsFrame>();
  // The first eigenvalue is the largest one.
  auto const principal_axis =
      eigensystem.rotation(PrincipalComponentsAxis({1, 0, 0}));

  // Project the |vectors| on the principal axis.
  for (auto it = begin; it != end; ++it) {
    auto const& displacement = displacements_[it->index];
    it->projection = InnerProduct(principal_axis, displacement);
  }

  // A comparator for the projections.
  auto const projection_less = [](Index const& left, Index const& right) {
    return left.projection < right.projection;
  };

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
      {displacements_[mid_lower->index], displacements_[mid_upper->index]},
      {1, 1});

  auto first_child = BuildTree(begin, mid_upper, size / 2);
  auto second_child = BuildTree(mid_upper, end, size - size / 2);

  return make_not_null_unique<Node>(
      Internal{.principal_axis = principal_axis,
               .anchor = anchor,
               .children = {std::move(first_child), std::move(second_child)}});
}

template<typename Value_>
typename PrincipalComponentPartitioningTree<Value_>::
DisplacementSymmetricBilinearForm
PrincipalComponentPartitioningTree<Value_>::ComputePrincipalComponentForm(
    Indices::iterator const begin,
    Indices::iterator const end) const {
  DisplacementSymmetricBilinearForm result;
  for (auto it = begin; it != end; ++it) {
    auto const& displacement = displacements_[it->index];
    result += SymmetricProduct(displacement, displacement);
  }
  return result;
}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
