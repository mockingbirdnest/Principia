#pragma once

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"
#include <string>

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

template<typename D>
R3Element<Quantities::Quantity<D>> operator+(
    R3Element<Quantities::Quantity<D>> const& right);
template<typename D>
R3Element<Quantities::Quantity<D>> operator-(
    R3Element<Quantities::Quantity<D>> const& right);

template<typename D>
R3Element<Quantities::Quantity<D>> operator+(
    R3Element<Quantities::Quantity<D>> const& left,
    R3Element<Quantities::Quantity<D>> const& right);
template<typename D>
R3Element<Quantities::Quantity<D>> operator-(
    R3Element<Quantities::Quantity<D>> const& left,
    R3Element<Quantities::Quantity<D>> const& right);

template<typename D>
R3Element<Quantities::Quantity<D>> operator*(
    Quantities::Dimensionless const& left,
    R3Element<Quantities::Quantity<D>> const& right);
template<typename D>
R3Element<Quantities::Quantity<D>> operator*(
    R3Element<Quantities::Quantity<D>> const& left,
    Quantities::Dimensionless const& right);
template<typename D>
R3Element<Quantities::Quantity<D>> operator/(
    R3Element<Quantities::Quantity<D>> const& left,
    Quantities::Dimensionless const& right);

template<typename DLeft, typename DRight>
R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                              Quantities::Quantity<DRight>>> operator*(
    Quantities::Quantity<DRight> const& left,
    R3Element<Quantities::Quantity<DRight>> const& right);
template<typename DLeft, typename DRight>
R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                              Quantities::Quantity<DRight>>> operator*(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    Quantities::Quantity<DRight> const& right);
template<typename DLeft, typename DRight>
R3Element<Quantities::Quotient<Quantities::Quantity<DLeft>,
                               Quantities::Quantity<DRight>>> operator/(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    Quantities::Quantity<DRight> const& right);

template<typename D>
void operator+=(R3Element<Quantities::Quantity<D>>& left,
                R3Element<Quantities::Quantity<D>> const& right);
template<typename D>
void operator-=(R3Element<Quantities::Quantity<D>>& left,
                R3Element<Quantities::Quantity<D>> const& right);

template<typename D>
void operator*=(R3Element<Quantities::Quantity<D>>& left,
                Quantities::Dimensionless const& right);
template<typename D>
void operator/=(R3Element<Quantities::Quantity<D>>& left,
                Quantities::Dimensionless const& right);

template<typename Left, typename Right>
R3Element<Quantities::Product<Left, Right>> Cross(
    R3Element<Left> const& left,
    R3Element<Right> const& right);

template<typename DLeft, typename DRight>
Quantities::Product<Quantities::Quantity<DLeft>,
                    Quantities::Quantity<DRight>> Dot(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    R3Element<Quantities::Quantity<DRight>> const& right);

}  // namespace Geometry
}  // namespace Principia

#include "R3Element-body.hpp"
