#pragma once

#include "numerics/nearest_neighbour.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <type_traits>

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "geometry/barycentre_calculator.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _nearest_neighbour {
namespace internal {

using namespace principia::geometry::_barycentre_calculator;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;

constexpr std::int32_t no_min_index = -1;

template<affine Point_, typename Payload_>
PrincipalComponentPartitioningTree<Point_, Payload_>::
PrincipalComponentPartitioningTree(
    std::vector<not_null<Point const*>> const& points,
    std::int64_t const max_points_per_cell)
    : points_(points),
      max_points_per_cell_(max_points_per_cell) {
  CHECK_LE(points_.size(), std::numeric_limits<std::int32_t>::max());
  if constexpr (has_payload) {
    payloads_.resize(points_.size());
  }
  if (!points_.empty()) {
    Initialize();
  }
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Add(
    not_null<Point const*> const point) {
  if constexpr (has_payload) {
    payloads_.push_back(std::nullopt);
  }
  AddPoint(point);
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Add(
    not_null<Point const*> const point,
    Payload const& payload)
  requires has_payload {
  payloads_.push_back(payload);
  AddPoint(point);
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Optimize() {
  // Build a permutation of the points in the order they appear in the leaves,
  // and update the leaf indices accordingly.
  std::vector<std::int32_t> permutation;
  permutation.resize(points_.size());
  std::int32_t new_index = 0;
  for (auto& node : nodes_) {
    if (auto* const old_indices = std::get_if<Leaf>(&node);
        old_indices != nullptr) {
      for (auto& old_index : *old_indices) {
        permutation[new_index] = old_index;
        old_index = new_index;
        ++new_index;
      }
    }
  }

  // Apply the permutation to our vectors.  I stole this code from
  // https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095/.
  for (std::int32_t i = 0; i < permutation.size(); ++i) {
    auto current = i;
    while (i != permutation[current]) {
      auto const next = permutation[current];
      std::swap(points_[current], points_[next]);
      std::swap(displacements_[current], displacements_[next]);
      if constexpr (has_payload) {
        std::swap(payloads_[current], payloads_[next]);
      }
      permutation[current] = current;
      current = next;
    }
    permutation[current] = current;
  }
}
template<affine Point_, typename Payload_>
Point_ const*
PrincipalComponentPartitioningTree<Point_, Payload_>::FindNearestNeighbour(
    Point const& point,
    NormType const& radius,
    Filter const& filter) const {
  if (displacements_.empty()) {
    return nullptr;
  }
  Norm²Type const radius² = Pow<2>(radius);
  Norm²Type min_distance²;
  std::int32_t min_index;
  Find(point - centroid_,
       radius²,
       filter,
       /*parent=*/nullptr,
       *root_,
       min_distance², min_index,
       /*must_check_other_side=*/nullptr);

  // In the end, this is why we retain the points: we want to return a pointer
  // that the client gave us.
  return min_index == no_min_index
             ? nullptr
             : static_cast<Point const*>(points_[min_index]);
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Initialize() {
  // Compute the centroid of the points.
  BarycentreCalculator<Point, double> centroid;
  for (auto const point : points_) {
    centroid.Add(*point, 1);
  }
  centroid_ = centroid.Get();

  // Compute the displacements from the centroid.
  displacements_.reserve(points_.size());
  for (auto const point : points_) {
    displacements_.push_back(*point - centroid_);
  }

  // Fill indices to denote the displacements.
  Indices indices;
  indices.reserve(displacements_.size());
  for (int i = 0; i < displacements_.size(); ++i) {
    indices.push_back({.index = i, .projection = NormType{}});
  }

  // Finally, build the tree.
  root_ = BuildTree(indices.begin(), indices.end(), indices.size());
}

template<affine Point_, typename Payload_>
not_null<typename PrincipalComponentPartitioningTree<Point_, Payload_>::Node*>
PrincipalComponentPartitioningTree<Point_, Payload_>::BuildTree(
    typename Indices::iterator const begin,
    typename Indices::iterator const end,
    std::int64_t const size,
    Node* const node) {
  if (size <= max_points_per_cell_) {
    // We are done subdividing, return a leaf.  About half the leaves end up
    // being split, with the unsplit leaves being about 50% full, so this wastes
    // about 25% of the space to save allocations.
    Leaf leaf;
    leaf.reserve(max_points_per_cell_);
    for (auto it = begin; it != end; ++it) {
      leaf.push_back(it->index);
    }
    nodes_.push_back(std::move(leaf));
    return &nodes_.back();
  }

  // Compute the "inertia" of the selected displacements and diagonalize it.
  // This gives us the direction of the largest eigenpoint, which defines the
  // separator plane.
  auto const form = ComputePrincipalComponentForm(begin, end);
  auto const eigensystem =
      form.template Diagonalize<PrincipalComponentsFrame>();
  // The last eigenpoint is the largest one.
  auto const principal_axis =
      eigensystem.rotation(PrincipalComponentsAxis({0, 0, 1}));

  // Project the `vectors` on the principal axis.
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

  // Here `mid_upper` denotes an index that denotes the point that is "in the
  // middle" of the projections on the principal axis.  We do not wish to anchor
  // our separator plane on a point, because that would make it more likely that
  // queries have to dig into two branches of the tree.  Instead, we look for
  // the maximum projection in the "lower" half space, and anchor our plane
  // halfway between these points.  Cost: 0.25 * N.
  auto const mid_lower = std::max_element(begin, mid_upper, projection_less);

  auto const anchor = Barycentre(
      {displacements_[mid_lower->index], displacements_[mid_upper->index]});

  auto const first_child = BuildTree(begin, mid_upper, size / 2);
  auto const second_child = BuildTree(mid_upper, end, size - size / 2);

  Internal internal{.principal_axis = principal_axis,
                    .anchor = anchor,
                    .children = {first_child, second_child}};
  if (node == nullptr) {
    nodes_.push_back(std::move(internal));
    return &nodes_.back();
  } else {
    *node = std::move(internal);
    return node;
  }
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::AddPoint(
    not_null<Point const*> const point) {
  points_.push_back(point);
  if (points_.size() == 1) {
    Initialize();
  } else {
    // There is no good way to rebalance a PCP tree (or a kd tree for that
    // matter).  Bkd trees provide some kind of solution, but their lookup is
    // expensive.  We are assuming that the tree won't get so unbalanced that
    // lookups get costly.
    auto const displacement = *point - centroid_;
    std::int32_t const index = displacements_.size();
    displacements_.push_back(displacement);
    Add(index, *root_);
  }
}

template<affine Point_, typename Payload_>
typename PrincipalComponentPartitioningTree<Point_, Payload_>::
DisplacementSymmetricBilinearForm
PrincipalComponentPartitioningTree<Point_, Payload_>::
ComputePrincipalComponentForm(
    typename Indices::iterator const begin,
    typename Indices::iterator const end) const {
  DisplacementSymmetricBilinearForm result;
  for (auto it = begin; it != end; ++it) {
    result += SymmetricSquare(displacements_[it->index]);
  }
  return result;
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Add(
    std::int32_t const index,
    Node& node) {
  std::visit(
      [this, index, &node](auto&& v) {
        using V = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<V, Internal>) {
          Add(index, v);
        } else if constexpr (std::is_same_v<V, Leaf>) {
          Add(index, v, &node);
        } else {
          static_assert(false, "Unexpected node type");
        }
      },
      node);
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Add(
    std::int32_t const index,
    Internal const& internal) {
  NormType const projection = InnerProduct(
      internal.principal_axis, displacements_[index] - internal.anchor);
  if (projection < NormType{}) {
    Add(index, *internal.children.first);
  } else {
    Add(index, *internal.children.second);
  }
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Add(
    std::int32_t const index,
    Leaf& leaf,
    not_null<Node*> const node) {
  if (leaf.size() < max_points_per_cell_) {
    leaf.push_back(index);
  } else {
    // The leaf is full, we need to split it by building a (sub)tree based on
    // it.  This code is written is such a way that the `leaf` never has more
    // than `max_points_per_cell_` elements, the number that we preallocated.
    Indices indices;
    indices.reserve(max_points_per_cell_ + 1);
    indices.push_back({.index = index, .projection = NormType{}});
    for (const std::int32_t index : leaf) {
      indices.push_back({.index = index, .projection = NormType{}});
    }
    auto const subtree = BuildTree(indices.begin(),
                                   indices.end(),
                                   indices.size(),
                                   /*node=*/node);
    CHECK(std::holds_alternative<Internal>(*subtree));
  }
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Find(
    Displacement const& displacement,
    Norm²Type const& radius²,
    Filter const& filter,
    Internal const* const parent,
    Node const& node,
    Norm²Type& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  std::visit(
      [&](auto&& v) {
        using V = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<V, Internal>) {
          Find(displacement,
               radius²,
               filter,
               parent,
               v,
               min_distance²,
               min_index,
               must_check_other_side);
        } else if constexpr (std::is_same_v<V, Leaf>) {
          Find(displacement,
               radius²,
               filter,
               parent,
               v,
               min_distance²,
               min_index,
               must_check_other_side);
        } else {
          static_assert(false, "Unexpected node type");
        }
      },
      node);
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Find(
    Displacement const& displacement,
    Norm²Type const& radius²,
    Filter const& filter,
    Internal const* parent,
    Internal const& internal,
    Norm²Type& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  NormType const projection =
      InnerProduct(internal.principal_axis, displacement - internal.anchor);

  // We first search on the "preferred" side of the separator plane, i.e., the
  // one where `displacement` is located.  We may have to look at the "other"
  // side, if `displacement` is too close to that plane.
  not_null<Node const*> const preferred_side = projection < NormType{}
                                                   ? internal.children.first
                                                   : internal.children.second;

  std::int32_t preferred_min_index;
  Norm²Type preferred_min_distance²;
  bool preferred_must_check_other_side;
  Find(displacement,
       radius²,
       filter,
       &internal,
       *preferred_side,
       preferred_min_distance², preferred_min_index,
       &preferred_must_check_other_side);

  if (preferred_must_check_other_side) {
    not_null<Node const*> const other_side = projection < NormType{}
                                                 ? internal.children.second
                                                 : internal.children.first;

    // We omit `must_check_other_side` because there is no point in checking the
    // other side again.
    // « Maintenant, rien qu'en traversant, j'ai fait qu'en face soit en face.
    // L'autre côté a changé de côté! » [GT71]
    std::int32_t other_min_index;
    Norm²Type other_min_distance²;
    Find(displacement,
         radius²,
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

  // Finally, check if we are close to the separator plane of the `parent`, in
  // which case *it* will need to check the other side of its plane.
  if (parent != nullptr && must_check_other_side != nullptr) {
    Norm²Type const projection² = Pow<2>(
        InnerProduct(parent->principal_axis, displacement - parent->anchor));
    *must_check_other_side =
        projection² < min_distance² && projection² <= radius²;
  }
}

template<affine Point_, typename Payload_>
void PrincipalComponentPartitioningTree<Point_, Payload_>::Find(
    Displacement const& displacement,
    Norm²Type const& radius²,
    Filter const& filter,
    Internal const* const parent,
    Leaf const& leaf,
    Norm²Type& min_distance²,
    std::int32_t& min_index,
    bool* const must_check_other_side) const {
  CHECK(!leaf.empty());

  // Find the point in this leaf which is the closest to `displacement`.
  min_distance² = Infinity<Norm²Type>;
  min_index = no_min_index;
  // If `Optimize()` was called, this loop accesses sequential data in
  // `displacements_` and should be cache-friendly.
  for (auto const index : leaf) {
    auto const distance² = (displacements_[index] - displacement).Norm²();
    // Skip the points that are filtered out.  Note that *all* the points may be
    // filtered out.  The order of comparisons matters because `min_distance²`
    // is typically small.
    if (distance² < min_distance² && distance² <= radius²) {
      bool must_update;
      if (filter == nullptr) {
        must_update = true;
      } else if constexpr (has_payload) {
        must_update = filter(points_[index], payloads_[index]);
      } else {
        must_update = filter(points_[index]);
      }
      if (must_update) {
        min_distance² = distance²;
        min_index = index;
      }
    }
  }

  // If there is a parent node, and `displacement` is very close to the
  // separator plane, we'll need to check the other side of the plane.
  if (parent != nullptr && must_check_other_side != nullptr) {
    Norm²Type const projection² = Pow<2>(
        InnerProduct(parent->principal_axis, displacement - parent->anchor));
    *must_check_other_side =
        projection² < min_distance² && projection² <= radius²;
  }
}

}  // namespace internal
}  // namespace _nearest_neighbour
}  // namespace numerics
}  // namespace principia
