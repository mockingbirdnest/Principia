#pragma once

#include "Grassmann.hpp"

namespace principia {
namespace geometry {

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
class LinearMap {
 public:
  LinearMap() = default;
  virtual ~LinearMap() = default;

  virtual Multivector<Scalar, ToFrame, Rank> ActOn(
      Multivector<Scalar, FromFrame, Rank> const& multivector) = 0;
};

}  // namespace geometry
}  // namespace principia
