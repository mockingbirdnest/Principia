#pragma once

namespace principia {
namespace Geometry {

template<typename LeftScalar, typename RightScalar, typename Frame>
inline quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Vector<LeftScalar, Frame> const& left,
    Vector<RightScalar, Frame> const& right) {
  return Dot(left.coordinates, right.coordinates);
}
template<typename LeftScalar, typename RightScalar, typename Frame>
inline quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Bivector<LeftScalar, Frame> const& left,
    Bivector<RightScalar, Frame> const& right) {
  return Dot(left.coordinates, right.coordinates);
}
template<typename LeftScalar, typename RightScalar, typename Frame>
inline quantities::Product<LeftScalar, RightScalar> InnerProduct(
    Trivector<LeftScalar, Frame> const& left,
    Trivector<RightScalar, Frame> const& right) {
  return left.coordinates * right.coordinates;
}

template<typename LeftScalar, typename RightScalar, typename Frame>
inline Bivector<quantities::Product<LeftScalar, RightScalar>,
                Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                             Vector<RightScalar, Frame> const& right) {
  return Bivector<quantities::Product<LeftScalar, RightScalar>,
                   Frame>(Cross(left.coordinates, right.coordinates));
}
template<typename LeftScalar, typename RightScalar, typename Frame>
inline Trivector<quantities::Product<LeftScalar, RightScalar>,
                 Frame> Wedge(Bivector<LeftScalar, Frame> const& left,
                              Vector<RightScalar, Frame> const& right) {
  return Trivector<quantities::Product<LeftScalar, RightScalar>,
                   Frame>(Dot(left.coordinates, right.coordinates));
}
template<typename LeftScalar, typename RightScalar, typename Frame>
inline Trivector<quantities::Product<LeftScalar, RightScalar>,
                 Frame> Wedge(Vector<LeftScalar, Frame> const& left,
                              Bivector<RightScalar, Frame> const& right) {
  return Trivector<quantities::Product<LeftScalar, RightScalar>,
                   Frame>(Dot(left.coordinates, right.coordinates));
}

// Lie bracket on V ^ V = so(V).
template<typename LeftScalar, typename RightScalar, typename Frame>
inline Bivector<quantities::Product<LeftScalar, RightScalar>,
                Frame> Commutator(Bivector<LeftScalar, Frame> const& left,
                                  Bivector<RightScalar, Frame> const& right) {
  return Bivector<quantities::Product<LeftScalar, RightScalar>,
                  Frame>(Cross(left.coordinates, right.coordinates));
}

// Left action of V ^ V = so(V) on V.
template<typename LeftScalar, typename RightScalar, typename Frame>
inline Vector<quantities::Product<LeftScalar, RightScalar>,
              Frame> operator*(Bivector<LeftScalar, Frame> const& left,
                               Vector<RightScalar, Frame> const& right) {
  return Vector<quantities::Product<LeftScalar, RightScalar>,
                Frame>(Cross(left.coordinates, right.coordinates));
}

// Right action of V ^ V = so(V) on V* = V.
template<typename LeftScalar, typename RightScalar, typename Frame>
inline Vector<quantities::Product<LeftScalar, RightScalar>,
              Frame> operator*(Vector<LeftScalar, Frame> const& left,
                               Bivector<RightScalar, Frame> const& right) {
  return Vector<quantities::Product<LeftScalar, RightScalar>,
                Frame>(Cross(left.coordinates, right.coordinates));
}

template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator+(Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(+right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator-(Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(-right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>

inline Multivector<T, Frame,
                   Rank> operator+(Multivector<T, Frame, Rank> const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates + right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator-(Multivector<T, Frame, Rank> const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates - right.coordinates);
}

template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator*(quantities::Dimensionless const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left * right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator*(Multivector<T, Frame, Rank> const& left,
                                   quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates * right);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator/(Multivector<T, Frame, Rank> const& left,
                                   quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates / right);
}

template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<quantities::Product<U, T> , Frame,
                   Rank> operator*(U const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left * right.coordinates);
}
template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<quantities::Product<T, U>, Frame,
                   Rank> operator*(Multivector<T, Frame, Rank> const& left,
                                   quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates * right);
}
template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<quantities::Quotient<T, U>, Frame,
                   Rank> operator/(Multivector<T, Frame, Rank> const& left,
                                   quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates / right);
}

template<typename T, typename Frame, unsigned int Rank>
inline void operator+=(Multivector<T, Frame, Rank>& left,
                       Multivector<T, Frame, Rank> const& right) {
  left.coordinates += right.coordinates;
}
template<typename T, typename Frame, unsigned int Rank>
inline void operator-=(Multivector<T, Frame, Rank>& left,
                       Multivector<T, Frame, Rank> const& right) {
  left.coordinates -= right.coordinates;
}
template<typename T, typename Frame, unsigned int Rank>
inline void operator*=(Multivector<T, Frame, Rank>& left,
                       quantities::Dimensionless const& right) {
  left.coordinates *= right;
}
template<typename T, typename Frame, unsigned int Rank>
inline void operator/=(Multivector<T, Frame, Rank>& left,
                       quantities::Dimensionless const& right) {
  left.coordinates /= right;
}

}  // namespace Geometry
}  // namespace principia
