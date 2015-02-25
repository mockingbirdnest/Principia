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

}  // namespace geometry
}  // namespace principia
