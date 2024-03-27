#pragma once

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

  Weight const& weight() const;

 private:
  bool empty_ = true;
  Product<Difference<Point>, Weight> weighted_sum_;
  Weight weight_;

  // We need reference values to convert points into vectors, if needed.  We
  // pick default-constructed objects as they don't introduce any inaccuracies
  // in the computations.
  static Point const reference_;
};
/*
template<typename Vector, typename Weight>
  requires homogeneous_vector_space<Vector, Weight>
class BarycentreCalculator<Vector, Weight> final {
 public:
  BarycentreCalculator() = default;

  void Add(Vector const& vector, Weight const& weight);
  Vector Get() const;

  // The sum of the weights added so far.
  Weight const& weight() const;

 private:
  bool empty_ = true;
  Product<Vector, Weight> weighted_sum_;
  Weight weight_;
};
*/
// |T| is anything for which a specialization of BarycentreCalculator exists.
template<typename T, typename Weight>
T Barycentre(std::pair<T, T> const& ts,
             std::pair<Weight, Weight> const& weights);
template<typename T, typename Weight, template<typename...> class Container>
T Barycentre(Container<T> const& ts, Container<Weight> const& weights);

}  // namespace internal

using internal::Barycentre;
using internal::BarycentreCalculator;

}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia

#include "geometry/barycentre_calculator_body.hpp"
