// Dimensionless.ipp

#pragma once

namespace principia {
namespace quantities {
inline Dimensionless::Dimensionless(double const value) : value_(value) {}
inline double Dimensionless::Value() const { return value_; }
template<int Exponent>
inline Dimensionless Dimensionless::Pow() const {
  return this->Pow(Exponent);
}
inline Dimensionless Dimensionless::Pow(int const exponent) const {
  if (exponent < 0) {
    return (1 / *this).Pow(-exponent);
  } else if (exponent == 0) { 
    return 1;
  } else if (exponent == 1) {
    return *this;
  } else if (exponent & 1) {
    return *this * (*this * *this).Pow((exponent - 1) / 2);
  } else {
    return (*this * *this).Pow(exponent / 2);
  }
}

inline Dimensionless operator+(Dimensionless const& right) {
  return +right.Value();
}
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

inline std::wstring ToString(Dimensionless const& number,
                             unsigned char const precision) {
  wchar_t result[50];
  std::swprintf(result, 49, (L"%."+ std::to_wstring(precision) + L"e").c_str(),
                number.Value());
  return result;
}
}  // namespace quantities
}  // namespace principia
