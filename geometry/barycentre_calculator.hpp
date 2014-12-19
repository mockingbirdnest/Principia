#pragma once

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
  Vector const Get() const;

 private:
  bool empty_ = true;
  decltype(std::declval<Vector>() * std::declval<Scalar>()) weighted_sum_;
  Scalar weight_;
};

}  // namespace geometry
}  // namespace principia

#include "geometry/barycentre_calculator_body.hpp"
