#pragma once

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"
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
  R3Element<Scalar> Coordinates;
};

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 3> {
  Multivector(Scalar coordinates) : coordinates(coordinates) {};
  Scalar Coordinates;
};

template<typename Scalar, typename Frame>
using Vector = Multivector<Scalar, Frame, 1>;

template<typename Scalar, typename Frame>
using Bivector = Multivector<Scalar, Frame, 2>;

template<typename Scalar, typename Frame>
using Trivector = Multivector<Scalar, Frame, 3>;

template<typename Scalar, typename Frame>
Scalar InnerProduct(Vector<Scalar, Frame> const& left,
                    Vector<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
Scalar InnerProduct(Bivector<Scalar, Frame> const& left,
                    Bivector<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
Scalar InnerProduct(Trivector<Scalar, Frame> const& left,
                    Trivector<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
Bivector<Scalar, Frame> Wedge(Vector<Scalar, Frame> const& Left,
                              Vector<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
Trivector<Scalar, Frame> Wedge(Bivector<Scalar, Frame> const& Left,
                               Vector<Scalar, Frame> const& right);

template<typename Scalar, typename Frame>
Trivector<Scalar, Frame> Wedge(Vector<Scalar, Frame> const& Left,
                               Bivector<Scalar, Frame> const& right);

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator+ (Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator- (Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>

Multivector<T, Frame,
            Rank> operator+ (Multivector<T, Frame, Rank> const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator- (Multivector<T, Frame, Rank> const& left,
                             Multivector<T, Frame, Rank> const& right);

template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator* (Quantities::Dimensionless const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator* (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
template<typename T, typename Frame, unsigned int Rank>
Multivector<T, Frame,
            Rank> operator/ (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);

template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<U, T> , Frame,
            Rank> operator* (U const& left,
                             Multivector<T, Frame, Rank> const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Product<T, U>, Frame,
            Rank> operator* (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);
template<typename T, typename U, typename Frame, unsigned int Rank>
Multivector<Quantities::Quotient<T, U>, Frame,
            Rank> operator/ (Multivector<T, Frame, Rank> const& left,
                             Quantities::Dimensionless const& right);

}
}

#include "Grassmann-body.hpp"
