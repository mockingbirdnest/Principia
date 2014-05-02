#pragma once

#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

#include "R3Element.hpp"

namespace Principia {
namespace Geometry {
template<typename Scalar, typename Frame, unsigned int Rank>
struct Multivector;

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 1> {
  Multivector(R3Element<Scalar> coordinates) : coordinates(coordinates) {};
  R3Element<Scalar> coordinates;
};

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 2> {
  Multivector(R3Element<Scalar> coordinates) : coordinates(coordinates) {};
  R3Element<Scalar> coordinates;
};

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 3> {
  Multivector(Scalar coordinates) : coordinates(coordinates) {};
  Scalar coordinates;
};

template<typename Scalar, typename Frame>
using Vector = Multivector<Scalar, Frame, 1>;
template<typename Scalar, typename Frame>
using Bivector = Multivector<Scalar, Frame, 2>;
template<typename Scalar, typename Frame>
using Trivector = Multivector<Scalar, Frame, 3>;

template<typename LeftScalar, typename RightScalar, typename Frame>
Quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Vector<LeftScalar, Frame> const& left,
    Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Bivector<LeftScalar, Frame> const& left,
    Bivector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Trivector<LeftScalar, Frame> const& left,
    Trivector<RightScalar, Frame> const& right);

template<typename LeftScalar, typename RightScalar, typename Frame>
Bivector<Quantities::Product<LeftScalar, RightScalar>,
         Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                      Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Trivector<Quantities::Product<LeftScalar, RightScalar>,
          Frame> Wedge(Bivector<LeftScalar, Frame> const& left,
                       Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Trivector<Quantities::Product<LeftScalar, RightScalar>,
          Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                       Bivector<RightScalar, Frame> const& right);

// Lie bracket on V ^ V = so(V).
template<typename LeftScalar, typename RightScalar, typename Frame>
Bivector<Quantities::Product<LeftScalar, RightScalar>,
         Frame> Commutator(Bivector<LeftScalar, Frame> const& left,
                           Bivector<RightScalar, Frame> const& right);

// Left action of V ^ V = so(V) on V.
template<typename LeftScalar, typename RightScalar, typename Frame>
Vector<Quantities::Product<LeftScalar, RightScalar>,
       Frame> operator*(Bivector<LeftScalar, Frame> const& left,
                        Vector<RightScalar, Frame> const& right);

// Right action of V ^ V = so(V) on V* = V.
template<typename LeftScalar, typename RightScalar, typename Frame>
Vector<Quantities::Product<LeftScalar, RightScalar>,
       Frame> operator*(Vector<LeftScalar, Frame> const& left,
                        Bivector<RightScalar, Frame> const& right);

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator+(Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator-(Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>

Multivector<T, Frame,
            Rank> operator+(Multivector<T, Frame, Rank> const& left,
                            Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator-(Multivector<T, Frame, Rank> const& left,
                            Multivector<T, Frame, Rank> const& right);

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator*(Quantities::Dimensionless const& left,
                            Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator*(Multivector<T, Frame, Rank> const& left,
                            Quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator/(Multivector<T, Frame, Rank> const& left,
                            Quantities::Dimensionless const& right);

template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<U, T> , Frame,
            Rank> operator*(U const& left,
                            Multivector<T, Frame, Rank> const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<T, U>, Frame,
            Rank> operator*(Multivector<T, Frame, Rank> const& left,
                            Quantities::Dimensionless const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Quotient<T, U>, Frame,
            Rank> operator/(Multivector<T, Frame, Rank> const& left,
                            Quantities::Dimensionless const& right);

template<typename T, typename Frame, unsigned int Rank>
void operator+=(Multivector<T, Frame, Rank>& left,
                Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator-=(Multivector<T, Frame, Rank>& left,
                Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator*=(Multivector<T, Frame, Rank>& left,
                Quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator/=(Multivector<T, Frame, Rank>& left,
                Quantities::Dimensionless const& right);

}  // namespace Geometry
}  // namespace Principia

#include "Grassmann-body.hpp"
