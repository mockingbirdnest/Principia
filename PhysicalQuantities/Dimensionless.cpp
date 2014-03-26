// Dimensionless.cpp

#include "stdafx.h"

#include "Dimensionless.h"

namespace PhysicalQuantities {

inline Dimensionless::Dimensionless(double value) : value_(value) {}
inline double Dimensionless::Value() const { return value_; }

inline Dimensionless operator+(Dimensionless const& right) { return right; }
inline Dimensionless operator-(Dimensionless const& right) {
  return -right.Value();
}
inline Dimensionless operator+(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.Value() + right.Value();
}
inline Dimensionless operator-(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.Value() - right.Value();
}
inline Dimensionless operator*(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.Value() * right.Value();
}
inline Dimensionless operator/(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.Value() / right.Value();
}

inline void operator+=(Dimensionless& left, Dimensionless const& right) {
  left = left + right;
}
inline void operator-=(Dimensionless& left, Dimensionless const& right) {
  left = left - right;
}
inline void operator*=(Dimensionless& left, Dimensionless const& right) {
  left = left * right;
}
inline void operator/=(Dimensionless& left, Dimensionless const& right) {
  left = left / right;
}

inline bool operator>(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() > right.Value();
}
inline bool operator<(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() < right.Value();
}
inline bool operator>=(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() >= right.Value();
}
inline bool operator<=(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() <= right.Value();
}

}