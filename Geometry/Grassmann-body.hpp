#pragma once

namespace Principia {
namespace Geometry {

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
  return Multivector<T, Frame, Rank>(left.coordinates + left.coordinates);
}

template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator*(Quantities::Dimensionless const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left * right.coordinates);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator*(Multivector<T, Frame, Rank> const& left,
                                   Quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates * right);
}
template<typename T, typename Frame, unsigned int Rank>
inline Multivector<T, Frame,
                   Rank> operator/(Multivector<T, Frame, Rank> const& left,
                                   Quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates / right);
}

template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<Quantities::Product<U, T> , Frame,
                   Rank> operator*(U const& left,
                                   Multivector<T, Frame, Rank> const& right) {
  return Multivector<T, Frame, Rank>(left * right.coordinates);
}
template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<Quantities::Product<T, U>, Frame,
                   Rank> operator*(Multivector<T, Frame, Rank> const& left,
                                   Quantities::Dimensionless const& right) {
  return Multivector<T, Frame, Rank>(left.coordinates * right);
}
template<typename T, typename U, typename Frame, unsigned int Rank>
inline Multivector<Quantities::Quotient<T, U>, Frame,
                   Rank> operator/(Multivector<T, Frame, Rank> const& left,
                                   Quantities::Dimensionless const& right) {
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
                       Quantities::Dimensionless const& right) {
  left.coordinates *= right;
}


template<typename T, typename Frame, unsigned int Rank>
inline void operator/=(Multivector<T, Frame, Rank>& left,
                       Quantities::Dimensionless const& right) {
  left.coordinates /= right;
}


}
}
