#pragma once

#include "numerics/nearest_neighbour.hpp"

#include <algorithm>
#include <limits>
#include <vector>
#include <type_traits>

#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_nearest_neighbour {

using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;

constexpr std::int32_t no_min_index = -1;

template<typename Value_>
PrincipalComponentPartitioningTree<Value_>::PrincipalComponentPartitioningTree(
    std::vector<not_null<Value const*>> const& values,
    std::int64_t const max_values_per_cell)
    : values_(values),
      max_values_per_cell_(max_values_per_cell) {
  CHECK_LE(values_.size(), std::numeric_limits<std::int32_t>::max());
  if (!values_.empty()) {
    Initialize();
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(
    not_null<Value const*> const value) {
  values_.push_back(value);
  if (values_.size() == 1) {
    Initialize();
  } else {
    // There is no good way to rebalance a PCP tree (or a kd tree for that
    // matter).  Bkd trees provide some kind of solution, but their lookup is
    // expensive.  We are assuming that the tree won't get so unbalanced that
    // lookups get costly.
    auto const displacement = *value - centroid_;
    std::int32_t const index = displacements_.size();
    displacements_.push_back(displacement);
    Add(index, *root_);
  }
}

template<typename Value_>
Value_ const* PrincipalComponentPartitioningTree<Value_>::FindNearestNeighbour(
    Value const& value,
    Filter const& filter) const {
  if (displacements_.empty()) {
    return nullptr;
  }
  Norm² min_distance²;
  std::int32_t min_index;
  Find(value - centroid_,
       filter,
       /*parent=*/nullptr,
       *root_,
       min_distance², min_index,
       /*must_check_other_side=*/nullptr);

  // In the end, this is why we retain the values: we want to return a pointer
  // that the client gave us.
  return min_index == no_min_index
             ? nullptr
             : static_cast<Value const*>(values_[min_index]);
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Initialize() {
  // Compute the centroid of the values.
  std::vector<Value> values_for_barycentre;
  values_for_barycentre.reserve(values_.size());
  for (auto const value : values_) {
    values_for_barycentre.push_back(*value);
  }
  centroid_ = Barycentre<Value, double, std::vector>(
      values_for_barycentre,
      std::vector<double>(values_for_barycentre.size(), 1));

  // Compute the displacements from the centroid.
  displacements_.reserve(values_.size());
  for (auto const value : values_) {
    displacements_.push_back(*value - centroid_);
  }

  // Fill indices to denote the displacements.
  Indices indices;
  indices.reserve(displacements_.size());
  for (int i = 0; i < displacements_.size(); ++i) {
    indices.push_back({.index = i, .projection = Norm{}});
  }

  // Finally, build the tree.
  root_ = BuildTree(indices.begin(), indices.end(), indices.size());
}

template<typename Value_>
not_null<std::unique_ptr<
typename PrincipalComponentPartitioningTree<Value_>::Node>>
PrincipalComponentPartitioningTree<Value_>::BuildTree(
    typename Indices::iterator const begin,
    typename Indices::iterator const end,
    std::int64_t const size) const {
  if (size <= max_values_per_cell_) {
    // We are done subdividing, return a leaf.
    Leaf leaf;
    leaf.reserve(size);
    for (auto it = begin; it != end; ++it) {
      leaf.push_back(it->index);
    }
    return make_not_null_unique<Node>(leaf);
  }

  // Compute the "inertia" of the selected displacements and diagonalize it.
  // This gives us the direction of the largest eigenvalue, which defines the
  // separator plane.
  auto const form = ComputePrincipalComponentForm(begin, end);
  auto const eigensystem =
      form.template Diagonalize<PrincipalComponentsFrame>();
  // The last eigenvalue is the largest one.
  auto const principal_axis =
      eigensystem.rotation(PrincipalComponentsAxis({0, 0, 1}));

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
    typename Indices::iterator const begin,
    typename Indices::iterator const end) const {
  DisplacementSymmetricBilinearForm result;
  for (auto it = begin; it != end; ++it) {
    result += SymmetricSquare(displacements_[it->index]);
  }
  return result;
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(std::int32_t const index,
                                                     Node& node) {
  if (std::holds_alternative<Internal>(node)) {
    Add(index, std::get<Internal>(node));
  } else if (std::holds_alternative<Leaf>(node)) {
    Add(index, std::get<Leaf>(node), node);
  } else {
    LOG(FATAL) << "Unexpected node";
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(std::int32_t const index,
                                                     Internal const& internal) {
  Norm const projection = InnerProduct(internal.principal_axis,
                                       displacements_[index] - internal.anchor);
  if (projection < Norm{}) {
    Add(index, *internal.children.first);
  } else {
    Add(index, *internal.children.second);
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(std::int32_t const index,
                                                     Leaf& leaf,
                                                     Node& node) {
  leaf.push_back(index);
  if (leaf.size() > max_values_per_cell_) {
    // The leaf is full, we need to split it by building a (sub)tree based on
    // it.
    Indices indices;
    indices.reserve(leaf.size());
    for (const std::int32_t index : leaf) {
      indices.push_back({.index = index, .projection = Norm{}});
    }
    auto const subtree =
        BuildTree(indices.begin(), indices.end(), indices.size());
    CHECK(std::holds_alternative<Internal>(*subtree));
    node = std::move(*subtree);
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Find(
    Displacement const& displacement,
    Filter const& filter,
    Internal const* const parent,
    Node const& node,
    Norm²& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  if (std::holds_alternative<Internal>(node)) {
    Find(displacement,
         filter,
         parent,
         std::get<Internal>(node),
         min_distance², min_index,
         must_check_other_side);
  } else if (std::holds_alternative<Leaf>(node)) {
    Find(displacement,
         filter,
         parent,
         std::get<Leaf>(node),
         min_distance², min_index,
         must_check_other_side);
  } else {
    LOG(FATAL) << "Unexpected node";
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Find(
    Displacement const& displacement,
    Filter const& filter,
    Internal const* parent,
    Internal const& internal,
    Norm²& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  Norm const projection =
      InnerProduct(internal.principal_axis, displacement - internal.anchor);

  // We first search on the "preferred" side of the separator plane, i.e., the
  // one where |displacement| is located.  We may have to look at the "other"
  // side, if |displacement| is too close to that plane.
  Node const* preferred_side;
  if (projection < Norm{}) {
    preferred_side = internal.children.first.get();
  } else {
    preferred_side = internal.children.second.get();
  }

  std::int32_t preferred_min_index;
  Norm² preferred_min_distance²;
  bool preferred_must_check_other_side;
  Find(displacement,
       filter,
       &internal,
       *preferred_side,
       preferred_min_distance², preferred_min_index,
       &preferred_must_check_other_side);

  if (preferred_must_check_other_side) {
    Node const* other_side;
    if (projection < Norm{}) {
      other_side = internal.children.second.get();
    } else {
      other_side = internal.children.first.get();
    }

    // We omit |must_check_other_side| because there is no point in checking the
    // other side again.
    // « Maintenant, rien qu'en traversant, j'ai fait qu'en face soit en face.
    // L'autre côté a changé de côté! » [GT71]
    std::int32_t other_min_index;
    Norm² other_min_distance²;
    Find(displacement,
         filter,
         parent,
         *other_side,
         other_min_distance², other_min_index,
         /*must_check_other_side=*/nullptr);

    if (other_min_distance² < preferred_min_distance²) {
      min_distance² = other_min_distance²;
      min_index = other_min_index;
    } else {
      min_distance² = preferred_min_distance²;
      min_index = preferred_min_index;
    }
  } else {
    min_distance² = preferred_min_distance²;
    min_index = preferred_min_index;
  }

  // Finally, check if we are close to the separator plane of the |parent|, in
  // which case *it* will need to check the other side of its plane.
  if (parent != nullptr && must_check_other_side != nullptr) {
    Norm² const projection² = Pow<2>(
        InnerProduct(parent->principal_axis, displacement - parent->anchor));
    *must_check_other_side = projection² < min_distance²;
  }
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Find(
    Displacement const& displacement,
    Filter const& filter,
    Internal const* const parent,
    Leaf const& leaf,
    Norm²& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  CHECK(!leaf.empty());

  // Find the point in this leaf which is the closest to |displacement|.
  min_distance² = Infinity<Norm²>;
  min_index = no_min_index;
  for (auto const index : leaf) {
    auto const distance² = (displacements_[index] - displacement).Norm²();
    // Skip the values that are filtered out.  Note that *all* the values may be
    // filtered out.
    if (distance² < min_distance² &&
        (filter == nullptr || filter(values_[index]))) {
      min_distance² = distance²;
      min_index = index;
    }
  }

  // If there is a parent node, and |displacement| is very close to the
  // separator plane, we'll need to check the other side of the plane.
  if (parent != nullptr && must_check_other_side != nullptr) {
    Norm² const projection² = Pow<2>(
        InnerProduct(parent->principal_axis, displacement - parent->anchor));
    *must_check_other_side = projection² < min_distance²;
  }
}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
