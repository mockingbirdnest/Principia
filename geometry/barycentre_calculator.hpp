#pragma once

#include <optional>
#include <vector>
#include <utility>

#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _barycentre_calculator {
namespace internal {

using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
class BarycentreCalculator final {
 public:
  BarycentreCalculator() = default;

  void Add(Point const& point, Weight const& weight);
  Point Get() const;

  // The sum of the weights added so far.
  Weight const& weight() const;

 private:
  bool empty_ = true;
  Product<Difference<Point>, Weight> weighted_sum_;
  Weight weight_;

  // We need reference values to convert points into vectors, if needed.  We
  // pick default-constructed objects as they don't introduce any inaccuracies
  // in the computations.
  // If we have an additive group, Point the same as Difference<Point>, which is
  // a vector, so no reference is needed.
  static std::conditional_t<additive_group<Point>, std::nullopt_t, Point> const
      reference_;
};

template<affine Point, homogeneous_field Weight, std::size_t size>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Point Barycentre(Point const (&points)[size], Weight const (&weights)[size]);
template<real_affine_space Point, std::size_t size>
Point Barycentre(Point const (&points)[size], double const (&weights)[size]);
template<real_affine_space Point, std::size_t size>
Point Barycentre(Point const (&points)[size]);

}  // namespace internal

using internal::Barycentre;
using internal::BarycentreCalculator;

}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia

#include "geometry/barycentre_calculator_body.hpp"
