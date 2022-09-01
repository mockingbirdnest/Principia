#pragma once

#include "numerics/nearest_neighbour.hpp"

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
  std::vector<Difference<Value>> differences;
  differences.reserve(values.size());
  for (auto const& value : values) {
    differences.push_back(value - centroid_);
  }

  using Eigenframe = Frame<enum class EigenframeTag>;

  auto const form = ComputePrincipalComponentForm(differences);
  auto const eigensystem = form.Diagonalize<Eigenframe>();
  // The first eigenvalue is the largest one.
  auto const principal_plane =
      eigensystem.rotation(Bivector<double, Eigenframe>({1, 0, 0}));
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(Value const& value) {}

template<typename Value_>
Value_ const& PrincipalComponentPartitioningTree<Value_>::FindNearest(
    Value const& value) const {}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
