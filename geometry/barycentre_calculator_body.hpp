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
  for (int i = 0; i < size; ++i) {
    if constexpr (additive_group<Point>) {
      total += points[i];
    } else {
      total += points[i] - Point{};
    }
  }
  if constexpr (additive_group<Point>) {
    return total / size;
  } else {
    return Point{} + total / size;
  }
}

template<typename T, typename Scalar, template<typename...> class Container>
T Barycentre(Container<T> const& ts, Container<Scalar> const& weights) {
  CHECK_EQ(ts.size(), weights.size()) << "Ts and weights of unequal sizes";
  CHECK(!ts.empty()) << "Empty input";
  BarycentreCalculator<T, Scalar> calculator;
  auto ts_it = ts.begin();
  auto weights_it = weights.begin();
  for (;
       ts_it != ts.end() && weights_it != weights.end();
       ++ts_it, ++weights_it) {
    calculator.Add(*ts_it, *weights_it);
  }
  return calculator.Get();
}

}  // namespace internal
}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia
