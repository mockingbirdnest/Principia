#pragma once

#include <utility>
#include <variant>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_nearest_neighbour {

using base::not_null;
using geometry::Bivector;
using geometry::Frame;
using geometry::SymmetricProduct;
using quantities::Difference;

template<typename Value_>
class PrincipalComponentPartitioningTree {
 public:
  using Value = Value_;

  PrincipalComponentPartitioningTree(std::vector<Value> const& values);

  void Add(Value const& value);

  Value const& FindNearest(Value const& value) const;

 private:
  // A frame used to compute the principal components.
  using PrincipalComponentsFrame =
      Frame<enum class PrincipalComponentsFrameTag>;

  // A displacement from the centroid.
  using Displacement = Difference<Value>;

  // A number of useful types for building data structures and parameters/
  // results.
  using DisplacementSymmetricBilinearForm =
      decltype(SymmetricProduct(std::declval<Displacement>(),
                                std::declval<Displacement>()));
  using DisplacementPrincipalComponentsSystem =
      typename DisplacementSymmetricBilinearForm::template Eigensystem<
          PrincipalComponentsFrame>;
  // TODO(phl): Unclear if the Bivector thingy is a good idea.
  using PrincipalAxis =
      decltype(std::declval<DisplacementPrincipalComponentsSystem>().rotation(
          std::declval<Bivector<double, PrincipalComponentsFrame>>()));

  // These are indices in |displacements|, that we'll use to refer to actual
  // values.  We cannot use iterators because they would be invalidated by
  // |Add|.  We use 32-bit integers because these will be swapped a lot, so we
  // don't want to touch too much memory.
  using Indices = std::vector<std::int32_t>;

  struct Node;

  using Children = std::pair<not_null<std::unique_ptr<Node>>,
                             not_null<std::unique_ptr<Node>>>;

  using Leaf = std::vector<Displacement>;

  struct Node {
    PrincipalAxis principal_axis;
    Displacement anchor;
    std::variant<Children, Leaf> tree;
  };

  not_null<std::unique_ptr<Node>> BuildTree(Indices::iterator begin,
                                            Indices::iterator end,
                                            std::int64_t size) const;

   DisplacementSymmetricBilinearForm ComputePrincipalComponentForm(
      Indices::iterator begin,
      Indices::iterator end) const;

  // The centroid of the values passed at construction.
  Value centroid_;

  // The displacements from the centroid.
  std::vector<Displacement> displacements_;

  std::unique_ptr<Node> root_;
};

}  // namespace internal_nearest_neighbour

using internal_nearest_neighbour::PrincipalComponentPartitioningTree;

}  // namespace numerics
}  // namespace principia

#include "numerics/nearest_neighbour_body.hpp"
