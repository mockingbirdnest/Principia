#pragma once

#include "geometry/barycentre_calculator.hpp"

#include <utility>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace geometry {
namespace _barycentre_calculator {
namespace internal {

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
void BarycentreCalculator<Point, Weight>::Add(Point const& point,
                                              Weight const& weight) {
  auto const weighted_sum_diff = (point - reference_) * weight;
  if (empty_) {
    weighted_sum_ = weighted_sum_diff;
    weight_ = weight;
    empty_ = false;
  } else {
    weighted_sum_ += weighted_sum_diff;
    weight_ += weight;
  }
}

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Point BarycentreCalculator<Point, Weight>::Get() const {
  CHECK(!empty_) << "Empty BarycentreCalculator";
  return reference_ + weighted_sum_ / weight_;
}

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Weight const& BarycentreCalculator<Point, Weight>::weight() const {
  return weight_;
}

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Point const BarycentreCalculator<Point, Weight>::reference_;

template<typename T, typename Scalar>
T Barycentre(std::pair<T, T> const& ts,
             std::pair<Scalar, Scalar> const& weights) {
  BarycentreCalculator<T, Scalar> calculator;
  calculator.Add(ts.first, weights.first);
  calculator.Add(ts.second, weights.second);
  return calculator.Get();
}

template<typename T, typename Scalar, template<typename...> class Container>
T Barycentre(Container<T> const& ts, Container<Scalar> const& weights) {
  CHECK_EQ(ts.size(), weights.size()) << "Ts and weights of unequal sizes";
  CHECK(!ts.empty()) << "Empty input";
  BarycentreCalculator<T, Scalar> calculator;
  auto ts_it = ts.begin();
  auto weights_it = weights.begin();
  for (; ts_it != ts.end() && weights_it != weights.end();
       ++ts_it, ++weights_it) {
    calculator.Add(*ts_it, *weights_it);
  }
  return calculator.Get();
}

}  // namespace internal
}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia
