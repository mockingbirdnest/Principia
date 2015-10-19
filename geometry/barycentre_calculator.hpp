#pragma once

#include <vector>
#include <utility>

namespace principia {
namespace geometry {

// |Vector| must be a vector space over the field |Scalar|.
template<typename Vector, typename Scalar>
class BarycentreCalculator {
 public:
  BarycentreCalculator() = default;
  ~BarycentreCalculator() = default;

  void Add(Vector const& vector, Scalar const& weight);
  Vector Get() const;

 private:
  bool empty_ = true;
  decltype(std::declval<Vector>() * std::declval<Scalar>()) weighted_sum_;
  Scalar weight_;
};

// |T| is anything for which a specialization of BarycentreCalculator exists.
template<typename T, typename Scalar>
T Barycentre(std::vector<T> const& ts, std::vector<Scalar> const& weights);

}  // namespace geometry
}  // namespace principia

#include "geometry/barycentre_calculator_body.hpp"
