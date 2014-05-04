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

template<typename DLeft, typename Right>
R3Element<quantities::Product<quantities::Quantity<DLeft>, Right>> operator*(
    quantities::Quantity<DLeft> const& left,
    R3Element<Right> const& right);
template<typename Left, typename DRight>
R3Element<quantities::Product<Left, quantities::Quantity<DRight>>> operator*(
    R3Element<Left> const& left,
    quantities::Quantity<DRight> const& right);
template<typename Left, typename DRight>
R3Element<quantities::Quotient<Left, quantities::Quantity<DRight>>> operator/(
    R3Element<Left> const& left,
    quantities::Quantity<DRight> const& right);

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

template<typename Left, typename Right>
R3Element<quantities::Product<Left, Right>> Cross(
    R3Element<Left> const& left,
    R3Element<Right> const& right);

template<typename Left, typename Right>
quantities::Product<Left, Right> Dot(R3Element<Left> const& left,
                                     R3Element<Right> const& right);

}  // namespace geometry
}  // namespace principia

#include "R3Element-body.hpp"
