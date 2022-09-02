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

  struct Node;

  using Children = std::pair<not_null<std::unique_ptr<Node>>,
                             not_null<std::unique_ptr<Node>>>;

  using Leaf = std::vector<Displacement>;

  struct Node {
    PrincipalAxis principal_axis;
    Displacement anchor;
    std::variant<Children, Leaf> tree;
  };

  static not_null<std::unique_ptr<Node>> BuildTree();

  static DisplacementSymmetricBilinearForm ComputePrincipalComponentForm(
      std::vector<Displacement> const& displacements);

  Value centroid_;
};

}  // namespace internal_nearest_neighbour

using internal_nearest_neighbour::PrincipalComponentPartitioningTree;

}  // namespace numerics
}  // namespace principia

#include "numerics/nearest_neighbour_body.hpp"
