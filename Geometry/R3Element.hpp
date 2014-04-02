#pragma once

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"
#include <string>

namespace Principia {
namespace Geometry {
  template<typename Scalar>
  struct R3Element {
   public:
    R3Element(Scalar const& x, Scalar const& y, Scalar const& z) :
      x(x), y(y), z(z) {};

    Scalar&       operator[] (int const index);
    Scalar& const operator[] (int const index) const;

    Scalar x;
    Scalar y;
    Scalar z;
  };
template<typename T>

R3Element<T> operator+ (R3Element<T> const& right);
template<typename T>
R3Element<T> operator- (R3Element<T> const& right);
template<typename T>

R3Element<T> operator+ (R3Element<T> const& left, R3Element<T> const& right);
template<typename T>
R3Element<T> operator- (R3Element<T> const& left, R3Element<T> const& right);
template<typename T>

R3Element<T> operator* (Quantities::Dimensionless const& left,
                        R3Element<T> const& right);
template<typename T>
R3Element<T> operator* (R3Element<T> const& left,
                        Quantities::Dimensionless const& right);
template<typename T>
R3Element<T> operator/ (R3Element<T> const& left,
                        Quantities::Dimensionless const& right);

template<typename T, typename U>
R3Element<Quantities::Product<U, T>> operator* (U const& left,
                                                R3Element<T> const& right);
template<typename T, typename U>
R3Element<Quantities::Product<T, U>> operator* (R3Element<T> const& left,
                                                U const& right);
template<typename T, typename U>
R3Element<Quantities::Quotient<T, U>> operator/ (R3Element<T> const& left,
                                                 U const& right);

template<typename T>
R3Element<T> operator+= (R3Element<T>& left, R3Element<T> const& right);
template<typename T>
R3Element<T> operator-= (R3Element<T>& left, R3Element<T> const& right);

template<typename T>
R3Element<T> operator*= (R3Element<T>& left,
                         Quantities::Dimensionless const& right);
template<typename T>
R3Element<T> operator/= (R3Element<T>& left,
                         Quantities::Dimensionless const& right);

template<typename T, typename U>
R3Element<Quantities::Product<T, U>> Cross(R3Element<T> const& left,
                                           R3Element<U> const& right);
template<typename T, typename U>
Quantities::Product<T, U> Dot(R3Element<T> const& left,
                              R3Element<U> const& right);
}
}
