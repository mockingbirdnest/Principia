#pragma once

#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

namespace Principia {
namespace Geometry {
// R3Element is an element of a 3-dimensional dimensionful vector space on the
// field R, represented by Dimensionless. It is the underlying data type for
// the more advanced strongly typed structures of the Grassmann algebras and
// affine spaces.
template<typename Scalar>
struct R3Element {
 public:
  R3Element(Scalar const& x, Scalar const& y, Scalar const& z);

  Scalar&       operator[] (int const index);
  Scalar const& operator[] (int const index) const;

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
R3Element<Scalar> operator*(Quantities::Dimensionless const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            Quantities::Dimensionless const& right);
template<typename Scalar>
R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                            Quantities::Dimensionless const& right);

template<typename DLeft, typename Right>
R3Element<Quantities::Product<Quantities::Quantity<DLeft>, Right>> operator*(
    Quantities::Quantity<DLeft> const& left,
    R3Element<Right> const& right);
template<typename Left, typename DRight>
R3Element<Quantities::Product<Left, Quantities::Quantity<DRight>>> operator*(
    R3Element<Left> const& left,
    Quantities::Quantity<DRight> const& right);
template<typename Left, typename DRight>
R3Element<Quantities::Quotient<Left, Quantities::Quantity<DRight>>> operator/(
    R3Element<Left> const& left,
    Quantities::Quantity<DRight> const& right);

template<typename Scalar>
void operator+=(R3Element<Scalar>& left,
                R3Element<Scalar> const& right);
template<typename Scalar>
void operator-=(R3Element<Scalar>& left,
                R3Element<Scalar> const& right);

template<typename Scalar>
void operator*=(R3Element<Scalar>& left,
                Quantities::Dimensionless const& right);
template<typename Scalar>
void operator/=(R3Element<Scalar>& left,
                Quantities::Dimensionless const& right);

template<typename Left, typename Right>
R3Element<Quantities::Product<Left, Right>> Cross(
    R3Element<Left> const& left,
    R3Element<Right> const& right);

template<typename Left, typename Right>
Quantities::Product<Left, Right> Dot(R3Element<Left> const& left,
                                     R3Element<Right> const& right);

}  // namespace Geometry
}  // namespace Principia

#include "R3Element-body.hpp"
