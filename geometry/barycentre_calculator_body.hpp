#pragma once

#include "geometry/barycentre_calculator.hpp"

#include <utility>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace geometry {
namespace internal_barycentre_calculator {

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

template<typename Vector, typename Scalar>
Scalar const& BarycentreCalculator<Vector, Scalar>::weight() const {
  return weight_;
}

template<typename T, typename Scalar>
T Barycentre(std::pair<T, T> const & ts,
             std::pair<Scalar, Scalar> const & weights) {
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
  for (;
       ts_it != ts.end() && weights_it != weights.end();
       ++ts_it, ++weights_it) {
    calculator.Add(*ts_it, *weights_it);
  }
  return calculator.Get();
}

}  // namespace internal_barycentre_calculator
}  // namespace geometry
}  // namespace principia
