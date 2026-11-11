#pragma once

#include <functional>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

#include "base/algebra.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/symmetric_bilinear_form.hpp"

namespace principia {
namespace numerics {
namespace _nearest_neighbour {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_symmetric_bilinear_form;

// Principal component partitioning trees (PCP trees) are introduced by [WZ91]
// in the context of quantization.  Their use for nearest neighbour search was
// proposed by [ZJL02].  We choose to use this algorithm not because it is easy
// but because it is coordinate-free.
// The `Payload` may be used to store client-specific data that is useful for
// filtering.
template<euclidean Point_, typename Payload_ = std::nullopt_t>
class PrincipalComponentPartitioningTree {
  static constexpr bool has_payload = !std::is_same_v<Payload_, std::nullopt_t>;

 public:
  using Point = Point_;
  using Payload = Payload_;

  // A displacement from the centroid.
  using Displacement = Difference<Point>;

  // The type of the norm (and its square) of `Displacement`.
  using NormType = NormType<Displacement>;
  using Norm²Type = Norm²Type<Displacement>;

  using Filter = std::conditional_t<
      has_payload,
      std::function<bool(Point const*, std::optional<Payload> const&)>,
      std::function<bool(Point const*)>>;

  // We stop subdividing a cell when it contains `max_points_per_cell` or fewer
  // points.  This API takes (non-owning) pointers so that the client can relate
  // the points given here to the ones it gets from `FindNearestNeighbour`.  The
  // vector `points` may be empty, in which case the object is not usable until
  // the first point has been `Add`ed.
  PrincipalComponentPartitioningTree(
      std::vector<not_null<Point const*>> const& points,
      std::int64_t max_points_per_cell);

  // Adds a new point to the tree, restructuring it as needed.  Optionally
  // supports adding a payload, if the tree has one.  This is the only way to
  // add payloads.
  void Add(not_null<Point const*> point);
  void Add(not_null<Point const*> point, Payload const& payload)
    requires has_payload;

  // Finds the nearest neighbour of the given `point` at a distance less than or
  // equal to `radius`.  Returns nullptr if the tree is empty.  Only the points
  // for which `filter` returns true are considered.
  // Beware! if `filter` rejects many points and `radius` is large, the search
  // may end up looking at the entire tree.
  Point const* FindNearestNeighbour(Point const& point,
                                    NormType const& radius,
                                    Filter const& filter = nullptr) const;

 private:
  // A frame used to compute the principal components.
  using PrincipalComponentsFrame =
      Frame<struct PrincipalComponentsFrameTag>;

  // A unit vector corresponding to `Displacement`.
  using Axis = NormalizedType<Displacement>;

  // A form that operates on `Displacement`s.
  // NOTE(phl): We don't have SymmetricSquare for Bivector, so this effectively
  // means that Displacement is a Vector.
  using DisplacementSymmetricBilinearForm =
      decltype(SymmetricSquare(std::declval<Displacement>()));

  // The principal components (a.k.a. eigensystem) of the above form.
  using DisplacementPrincipalComponentsSystem =
      typename DisplacementSymmetricBilinearForm::template Eigensystem<
          PrincipalComponentsFrame>;

  // A unit vector in the principal components frame.
  using PrincipalComponentsAxis = decltype(
      std::declval<DisplacementPrincipalComponentsSystem>().rotation.Inverse()(
          std::declval<Axis>()));

  // The declarations of the tree structure.
  struct Internal;

  using Leaf = std::vector<std::int32_t>;  // Indices in `displacements_`.

  using Node = std::variant<Internal, Leaf>;

  using Children = std::pair<not_null<std::unique_ptr<Node>>,
                             not_null<std::unique_ptr<Node>>>;

  struct Internal {
    Axis principal_axis;
    Displacement anchor;
    Children children;
  };

  // The construction of the tree uses this type, which contains an index in the
  // `displacements_` array and storage for the projection of the corresponding
  // `Displacement` on the current principal axis.  We use 32-bit integers
  // because these objects will be swapped a lot, so the less memory we touch
  // the better.
  struct Index {
    std::int32_t index;
    NormType projection;
  };
  using Indices = std::vector<Index>;

  // Called when the first point is added to the tree to initialize the
  // `centroid_` and the `root_`.
  void Initialize();

  // Constructs a tree for the displacements given by the index range
  // [begin, end[.  `size` must be equal to `std::distance(begin, end)`, but is
  // passed by the caller for efficiency.
  not_null<std::unique_ptr<Node>> BuildTree(typename Indices::iterator begin,
                                            typename Indices::iterator end,
                                            std::int64_t size) const;

  // Returns the symmetric bilinear form that represents the "inertia" of the
  // displacements given by the index range [begin, end[.
  DisplacementSymmetricBilinearForm ComputePrincipalComponentForm(
      typename Indices::iterator begin,
      typename Indices::iterator end) const;

  // Adds `point` to the tree, initializing it if needed.  The caller must take
  // care of managing `payloads_` if needed.
  void AddPoint(not_null<Point const*> point);

  // Adds the point at the given `index` by first finding the location in the
  // subtree rooted at `node` where it would be located, and then adding it to
  // the leaf (if there is room) or splitting the leaf (if not).
  void Add(std::int32_t index, Node& node);

  // Specializations for internal nodes and leaves, respectively.
  void Add(std::int32_t index, Internal const& internal);
  void Add(std::int32_t index, Leaf& leaf, Node& node);

  // Finds the point closest to `displacement` in the `node` and its children,
  // and returns its index and its (squared) distance.  If `displacement` is
  // close to the separator plane of `parent`, sets `must_check_other_side` to
  // true.  That pointer may be null if the client doesn't want to check this
  // condition.  `parent` should be null for the root of the tree.
  void Find(Displacement const& displacement,
            Norm²Type const& radius²,
            Filter const& filter,
            Internal const* parent,
            Node const& node,
            Norm²Type& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;

  // Specializations for internal nodes and leaves, respectively.
  void Find(Displacement const& displacement,
            Norm²Type const& radius²,
            Filter const& filter,
            Internal const* parent,
            Internal const& internal,
            Norm²Type& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;
  void Find(Displacement const& displacement,
            Norm²Type const& radius²,
            Filter const& filter,
            Internal const* parent,
            Leaf const& leaf,
            Norm²Type& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;

  // Construction parameters.
  std::vector<not_null<Point const*>> points_;
  std::int64_t const max_points_per_cell_;

  // The centroid of the points passed at construction, or the first point
  // passed to `Add` if no point is passed at construction.
  Point centroid_;

  // The displacements from the centroid.  The indices are the same as for
  // `points_`.
  std::vector<Displacement> displacements_;

  // Payloads.  Note that the not-a-payload type, `std::nullopt_t` cannot be
  // stored in an `optional`, so it helps detect mistakes.
  std::conditional_t<has_payload,
                     std::vector<std::optional<Payload>>,
                     std::monostate>
      payloads_;

  std::unique_ptr<Node> root_;
};

}  // namespace internal

using internal::PrincipalComponentPartitioningTree;

}  // namespace _nearest_neighbour
}  // namespace numerics
}  // namespace principia

#include "numerics/nearest_neighbour_body.hpp"
