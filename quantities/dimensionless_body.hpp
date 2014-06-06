#pragma once

#include <algorithm>
#include <cmath>
#include <string>

namespace principia {
namespace quantities {

inline Dimensionless::Dimensionless() : value_(0) {}
inline Dimensionless::Dimensionless(double const value) : value_(value) {}

inline Dimensionless Dimensionless::SIUnit() {
  return 1;
}

inline double Dimensionless::value() const {
  return value_;
}

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
  return right;
}

inline Dimensionless operator-(Dimensionless const& right) {
  return -right.value_;
}

inline Dimensionless operator+(double const left,
                               Dimensionless const& right) {
  return left + right.value_;
}
inline Dimensionless operator+(Dimensionless const& left,
                               double const right) {
  return left.value_ + right;
}
inline Dimensionless operator+(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.value_ + right.value_;
}

inline Dimensionless operator-(double const left,
                               Dimensionless const& right) {
  return left - right.value_;
}
inline Dimensionless operator-(Dimensionless const& left,
                               double const right) {
  return left.value_ - right;
}
inline Dimensionless operator-(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.value_ - right.value_;
}

inline Dimensionless operator*(double const left,
                               Dimensionless const& right) {
  return left * right.value_;
}
inline Dimensionless operator*(Dimensionless const& left,
                               double const right) {
  return left.value_ * right;
}
inline Dimensionless operator*(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.value_ * right.value_;
}

inline Dimensionless operator/(double const left,
                               Dimensionless const& right) {
  return left / right.value_;
}
inline Dimensionless operator/(Dimensionless const& left,
                               double const right) {
  return left.value_ / right;
}
inline Dimensionless operator/(Dimensionless const& left,
                               Dimensionless const& right) {
  return left.value_ / right.value_;
}

inline void operator+=(Dimensionless& left,  // NOLINT(runtime/references)
                       double const right) {
  left.value_ += right;
}
inline void operator+=(Dimensionless& left,  // NOLINT(runtime/references)
                       Dimensionless const& right) {
  left.value_ += right.value_;
}

inline void operator-=(Dimensionless& left,  // NOLINT(runtime/references)
                       double const right) {
  left.value_ -= right;
}
inline void operator-=(Dimensionless& left,  // NOLINT(runtime/references)
                       Dimensionless const& right) {
  left.value_ -= right.value_;
}

inline void operator*=(Dimensionless& left,  // NOLINT(runtime/references)
                       double const right) {
  left.value_ *= right;
}
inline void operator*=(Dimensionless& left,  // NOLINT(runtime/references)
                       Dimensionless const& right) {
  left.value_ *= right.value_;
}

inline void operator/=(Dimensionless& left,  // NOLINT(runtime/references)
                       double const right) {
  left.value_ /= right;
}
inline void operator/=(Dimensionless& left,  // NOLINT(runtime/references)
                       Dimensionless const& right) {
  left.value_ /= right.value_;
}

inline bool operator<(double const left, Dimensionless const& right) {
  return left < right.value_;
}

inline bool operator<(Dimensionless const& left, double const right) {
  return left.value_ < right;
}

inline bool operator<(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ < right.value_;
}

inline bool operator>(double const left, Dimensionless const& right) {
  return left > right.value_;
}

inline bool operator>(Dimensionless const& left, double const right) {
  return left.value_ > right;
}

inline bool operator>(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ > right.value_;
}

inline bool operator<=(double const left, Dimensionless const& right) {
  return left <= right.value_;
}

inline bool operator<=(Dimensionless const& left, double const right) {
  return left.value_ <= right;
}

inline bool operator<=(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ <= right.value_;
}

inline bool operator>=(double const left, Dimensionless const& right) {
  return left >= right.value_;
}

inline bool operator>=(Dimensionless const& left, double const right) {
  return left.value_ >= right;
}

inline bool operator>=(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ >= right.value_;
}

inline bool operator==(double const left, Dimensionless const& right) {
  return left == right.value_;
}

inline bool operator==(Dimensionless const& left, double const right) {
  return left.value_ == right;
}

inline bool operator==(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ == right.value_;
}

inline bool operator!=(double const left, Dimensionless const& right) {
  return left != right.value_;
}

inline bool operator!=(Dimensionless const& left, double const right) {
  return left.value_ != right;
}

inline bool operator!=(Dimensionless const& left, Dimensionless const& right) {
  return left.value_ != right.value_;
}

inline Dimensionless Abs(Dimensionless const& number) {
  return std::abs(number.value_);
}

inline std::string ToString(Dimensionless const& number,
                            unsigned char const precision) {
  char result[50];
  sprintf_s(result, ("%."+ std::to_string(precision) + "e").c_str(),
            number.value());
  return result;
}

inline std::ostream& operator<<(std::ostream& out,
                                Dimensionless const& number) {
  return out << ToString(number);
}

}  // namespace quantities
}  // namespace principia
