#pragma once

#include "..\Quantities\Dimensionless.hpp"
#include "..\Quantities\Quantities.hpp"
#include <string>

namespace Principia {
namespace Geometry {
  template<typename T>
  class R3Element {
   public:
    R3Element(T x, T y, T z) : x_(x), y_(y), z_(z) {};
    T X() const;
    void X(T x);
    T Y() const;
    void Y(T y);
    T Z() const;
    void Z(T z);
   private:
    T x_;
    T y_;
    T z_;
  };
template<typename T>

R3Element<T> operator+ (R3Element<T> right);
template<typename T>
R3Element<T> operator- (R3Element<T> right);
template<typename T>

R3Element<T> operator+ (R3Element<T> left, R3Element<T> right);
template<typename T>
R3Element<T> operator- (R3Element<T> left, R3Element<T> right);
template<typename T>

R3Element<T> operator* (Quantities::Dimensionless left, R3Element<T> right);
template<typename T>
R3Element<T> operator* (R3Element<T> left, Quantities::Dimensionless right);
template<typename T>
R3Element<T> operator/ (R3Element<T> left, Quantities::Dimensionless right);

template<typename T, typename Scalar>
R3Element<Quantities::Product<Scalar, T>> operator* (Scalar left,
                                                     R3Element<T> right);
template<typename T, typename Scalar>
R3Element<Quantities::Product<T, Scalar>> operator* (R3Element<T> left,
                                                     Scalar right);
template<typename T, typename Scalar>
R3Element<Quantities::Quotient<T, Scalar>> operator/ (R3Element<T> left,
                                                     Scalar right);
}
}
