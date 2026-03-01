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
  auto const weighted_sum_diff = [&]() {
    if constexpr (additive_group<Point>) {
      return point * weight;
    } else {
      return (point - reference_) * weight;
    }
  }();
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
  if constexpr (additive_group<Point>) {
    return weighted_sum_ / weight_;
  } else {
    return reference_ + weighted_sum_ / weight_;
  }
}

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Weight const& BarycentreCalculator<Point, Weight>::weight() const {
  return weight_;
}

template<affine Point, homogeneous_field Weight>
  requires homogeneous_vector_space<Difference<Point>, Weight>
std::conditional_t<additive_group<Point>, std::nullopt_t, Point> const
    BarycentreCalculator<Point, Weight>::reference_;

template<affine Point, homogeneous_field Weight, std::size_t size>
  requires homogeneous_vector_space<Difference<Point>, Weight>
Point Barycentre(Point const (&points)[size], Weight const (&weights)[size]) {
  static_assert(size != 0);
  BarycentreCalculator<Point, Weight> calculator;
  for (int i = 0; i < size; ++i) {
    calculator.Add(points[i], weights[i]);
  }
  return calculator.Get();
}

template<real_affine_space Point, std::size_t size>
Point Barycentre(Point const (&points)[size], double const (&weights)[size]) {
  return Barycentre<Point, double>(points, weights);
}

template<real_affine_space Point, std::size_t size>
Point Barycentre(Point const (&points)[size]) {
  static_assert(size != 0);
  Difference<Point> total{};
  static const Point origin{};
  for (int i = 0; i < size; ++i) {
    if constexpr (additive_group<Point>) {
      total += points[i];
    } else {
      total += points[i] - origin;
    }
  }
  if constexpr (additive_group<Point>) {
    return total / size;
  } else {
    return origin + total / size;
  }
}

}  // namespace internal
}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia
