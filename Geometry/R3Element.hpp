#pragma once

#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

namespace principia {
namespace geometry {

// R3Element is an element of a 3-dimensional dimensionful vector space on the
// field R, represented by Dimensionless. It is the underlying data type for
// the more advanced strongly typed structures of the Grassmann algebras and
// affine spaces.
template<typename Scalar>
struct R3Element {
 public:
  R3Element() = default;
  R3Element(Scalar const& x, Scalar const& y, Scalar const& z);

  Scalar&       operator[](int const index);
  Scalar const& operator[](int const index) const;

  Scalar x;
  Scalar y;
  Scalar z;
};

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator*(quantities::Dimensionless const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            quantities::Dimensionless const& right);
template<typename Scalar>
R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                            quantities::Dimensionless const& right);

template<typename LDimension, typename RScalar>
R3Element<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3Element<RScalar> const& right);
template<typename LScalar, typename RDimension>
R3Element<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right);
template<typename LScalar, typename RDimension>
R3Element<quantities::Quotient<LScalar, quantities::Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right);

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right);

template<typename Scalar>
void operator+=(R3Element<Scalar>& left,
                R3Element<Scalar> const& right);
template<typename Scalar>
void operator-=(R3Element<Scalar>& left,
                R3Element<Scalar> const& right);

template<typename Scalar>
void operator*=(R3Element<Scalar>& left,
                quantities::Dimensionless const& right);
template<typename Scalar>
void operator/=(R3Element<Scalar>& left,
                quantities::Dimensionless const& right);

template<typename LScalar, typename RScalar>
R3Element<quantities::Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right);

template<typename LScalar, typename RScalar>
quantities::Product<LScalar, RScalar> Dot(R3Element<LScalar> const& left,
                                          R3Element<RScalar> const& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/R3Element-body.hpp"
