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
  // These are members of a vector space, not |Vector| per se.
  std::vector<Difference<Value>> vectors;
  vectors.reserve(values.size());
  for (auto const& value : values) {
    vectors.push_back(value - centroid_);
  }

  // Compute the "inertia" of the |vectors| and diagonalize it.  This gives us
  // the direction of the largest eigenvalue, which defines the splitting plane.
  using Eigenframe = Frame<enum class EigenframeTag>;

  auto const form = ComputePrincipalComponentForm(vectors);
  auto const eigensystem = form.Diagonalize<Eigenframe>();
  // The first eigenvalue is the largest one.
  auto const principal_plane =
      eigensystem.rotation(Bivector<double, Eigenframe>({1, 0, 0}));

  // Project the |vectors| on the principal axis.
  std::vector<Difference<Value>> projections;
  projections.reserve(vectors.size());
  for (auto const& v : vectors) {
    auto const projection_on_plane = Wedge(principal_plane, v);
    projections.push_back(principal_plane * projection_on_plane)
  }

  // Find the median in linear time.  This is where our plane will lie.  Note
  // that |mid| is orthogonal to the plane.
  auto const mid = projections.begin() + projections.size() / 2;
  std::nth_element(projections.begin(), mid, projections.end());
  auto const anchor = *mid;
}

template<typename Value_>
void PrincipalComponentPartitioningTree<Value_>::Add(Value const& value) {}

template<typename Value_>
Value_ const& PrincipalComponentPartitioningTree<Value_>::FindNearest(
    Value const& value) const {}

}  // namespace internal_nearest_neighbour
}  // namespace numerics
}  // namespace principia
