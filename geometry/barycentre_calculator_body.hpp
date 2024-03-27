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

}  // namespace internal
}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia
