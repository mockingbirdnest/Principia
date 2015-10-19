#pragma once

#include "geometry/barycentre_calculator.hpp"

#include <vector>

#include "glog/logging.h"

namespace principia {
namespace geometry {

template<typename Vector, typename Scalar>
void BarycentreCalculator<Vector, Scalar>::Add(Vector const& vector,
                                               Scalar const& weight) {
  if (empty_) {
    weighted_sum_ = vector * weight;
    weight_ = weight;
    empty_ = false;
  } else {
    weighted_sum_ += vector * weight;
    weight_ += weight;
  }
}

template<typename Vector, typename Scalar>
Vector BarycentreCalculator<Vector, Scalar>::Get() const {
  CHECK(!empty_) << "Empty BarycentreCalculator";
  return Vector(weighted_sum_ / weight_);
}

template<typename T, typename Scalar>
T Barycentre(std::vector<T> const& ts, std::vector<Scalar> const& weights) {
  CHECK_EQ(ts.size(), weights.size()) << "Ts and weights of unequal sizes";
  CHECK(!ts.empty()) << "Empty input";
  BarycentreCalculator<T, Scalar> calculator;
  for (size_t i = 0; i < ts.size(); ++i) {
    calculator.Add(ts[i], weights[i]);
  }
  return calculator.Get();
}

}  // namespace geometry
}  // namespace principia
