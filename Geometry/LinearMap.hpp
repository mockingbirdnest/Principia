#pragma once

#include "Grassmann.hpp"

namespace Principia {
namespace Geometry {

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
class LinearMap {
  public:
    LinearMap() {}
    virtual ~LinearMap() {}

    virtual Multivector<Scalar, ToFrame, Rank> ActOn(
        Multivector<Scalar, FromFrame, Rank> multivector) = 0;
};

}  // namespace Geometry
}  // namespace Principia