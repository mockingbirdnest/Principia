// Dimensionless.ipp

#pragma once

namespace Quantities {
inline Dimensionless::Dimensionless(double value) : value_(value) {}
inline double Dimensionless::Value() const { return value_; }
// TODO(robin): This should not be inlined.
inline Dimensionless Exponentiate(Dimensionless const& base, 
                                  int const exponent) {
  if (exponent < 0) {
    return Exponentiate(1 / base, -exponent);
  } else if (exponent == 0) { 
    return 1;
  } else if (exponent == 1) {
    return base;
  } else if (exponent % 2 == 0) {
    return Exponentiate(base * base, exponent / 2);
  } else {
    return base * Exponentiate(base * base, (exponent - 1) / 2);
  }
}

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
inline bool operator==(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() == right.Value();
}
inline bool operator!=(Dimensionless const& left, Dimensionless const& right) {
  return left.Value() != right.Value();
}

inline Dimensionless Abs(Dimensionless const& number) {
  return std::abs(number.Value());
}

inline std::wstring ToString(Dimensionless const& number) {
  wchar_t result[50];
  std::swprintf(result, L"%.16e", number.Value());
  return result;
}
}