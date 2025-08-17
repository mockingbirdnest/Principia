#pragma once

#include <functional>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _nearest_neighbour {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::quantities::_arithmetic;

// Principal component partitioning trees (PCP trees) are introduced by [WZ91]
// in the context of quantization.  Their use for nearest neighbour search was
// proposed by [ZJL02].  We choose to use this algorithm not because it is easy
// but because it is coordinate-free.
template<typename Value_>
class PrincipalComponentPartitioningTree {
 public:
  using Value = Value_;

  using Filter = std::function<bool(Value const*)>;

  // We stop subdividing a cell when it contains `max_values_per_cell` or fewer
  // values.  This API takes (non-owning) pointers so that the client can relate
  // the values given here to the ones it gets from `FindNearestNeighbour`.  The
  // vector `values` may be empty, in which case the object is not usable until
  // the first value has been `Add`ed.
  PrincipalComponentPartitioningTree(
      std::vector<not_null<Value const*>> const& values,
      std::int64_t max_values_per_cell);

  // Adds a new value to the tree, restructuring it as needed.
  void Add(not_null<Value const*> value);

  // Finds the nearest neighbour of the given `value`.  Returns nullptr if the
  // tree is empty.  Only the values for which `filter` returns true are
  // considered.
  Value const* FindNearestNeighbour(Value const& value,
                                    Filter const& filter = nullptr) const;

 private:
  // A frame used to compute the principal components.
  using PrincipalComponentsFrame =
      Frame<struct PrincipalComponentsFrameTag>;

  // A displacement from the centroid.
  using Displacement = Difference<Value>;

  // The type of the norm (and its square) of `Displacement`.
  using Norm = typename Hilbert<Displacement>::NormType;
  using Norm² = typename Hilbert<Displacement>::Norm²Type;

  // A unit vector corresponding to `Displacement`.
  using Axis = typename Hilbert<Displacement>::NormalizedType;

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
    Norm projection;
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

  // Adds the value at the given `index` by first finding the location in the
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
            Filter const& filter,
            Internal const* parent,
            Node const& node,
            Norm²& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;

  // Specializations for internal nodes and leaves, respectively.
  void Find(Displacement const& displacement,
            Filter const& filter,
            Internal const* parent,
            Internal const& internal,
            Norm²& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;
  void Find(Displacement const& displacement,
            Filter const& filter,
            Internal const* parent,
            Leaf const& leaf,
            Norm²& min_distance²,
            std::int32_t& min_index,
            bool* must_check_other_side) const;

  // Construction parameters.
  std::vector<not_null<Value const*>> values_;
  std::int64_t const max_values_per_cell_;

  // The centroid of the values passed at construction, or the first value
  // passed to `Add` if no value is passed at construction.
  Value centroid_;

  // The displacements from the centroid.  The indices are the same as for
  // `values_`.
  std::vector<Displacement> displacements_;

  std::unique_ptr<Node> root_;
};

}  // namespace internal

using internal::PrincipalComponentPartitioningTree;

}  // namespace _nearest_neighbour
}  // namespace numerics
}  // namespace principia

#include "numerics/nearest_neighbour_body.hpp"
