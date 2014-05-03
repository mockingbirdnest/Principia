#pragma once

#include "Geometry/R3Element.hpp"
#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

namespace principia {
namespace geometry {
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
quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Vector<LeftScalar, Frame> const& left,
    Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Bivector<LeftScalar, Frame> const& left,
    Bivector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Trivector<LeftScalar, Frame> const& left,
    Trivector<RightScalar, Frame> const& right);

template<typename LeftScalar, typename RightScalar, typename Frame>
Bivector<quantities::Product<LeftScalar, RightScalar>,
         Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                      Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Trivector<quantities::Product<LeftScalar, RightScalar>,
          Frame> Wedge(Bivector<LeftScalar, Frame> const& left,
                       Vector<RightScalar, Frame> const& right);
template<typename LeftScalar, typename RightScalar, typename Frame>
Trivector<quantities::Product<LeftScalar, RightScalar>,
          Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                       Bivector<RightScalar, Frame> const& right);

// Lie bracket on V ^ V = so(V).
template<typename LeftScalar, typename RightScalar, typename Frame>
Bivector<quantities::Product<LeftScalar, RightScalar>,
         Frame> Commutator(Bivector<LeftScalar, Frame> const& left,
                           Bivector<RightScalar, Frame> const& right);

// Left action of V ^ V = so(V) on V.
template<typename LeftScalar, typename RightScalar, typename Frame>
Vector<quantities::Product<LeftScalar, RightScalar>,
       Frame> operator*(Bivector<LeftScalar, Frame> const& left,
                        Vector<RightScalar, Frame> const& right);

// Right action of V ^ V = so(V) on V* = V.
template<typename LeftScalar, typename RightScalar, typename Frame>
Vector<quantities::Product<LeftScalar, RightScalar>,
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
            Rank> operator*(quantities::Dimensionless const& left,
                            Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator*(Multivector<T, Frame, Rank> const& left,
                            quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator/(Multivector<T, Frame, Rank> const& left,
                            quantities::Dimensionless const& right);

template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<quantities::Product<U, T> , Frame,
            Rank> operator*(U const& left,
                            Multivector<T, Frame, Rank> const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<quantities::Product<T, U>, Frame,
            Rank> operator*(Multivector<T, Frame, Rank> const& left,
                            quantities::Dimensionless const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<quantities::Quotient<T, U>, Frame,
            Rank> operator/(Multivector<T, Frame, Rank> const& left,
                            quantities::Dimensionless const& right);

template<typename T, typename Frame, unsigned int Rank>
void operator+=(Multivector<T, Frame, Rank>& left,
                Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator-=(Multivector<T, Frame, Rank>& left,
                Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator*=(Multivector<T, Frame, Rank>& left,
                quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
void operator/=(Multivector<T, Frame, Rank>& left,
                quantities::Dimensionless const& right);

}  // namespace geometry
}  // namespace principia

#include "Grassmann-body.hpp"
