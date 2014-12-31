#pragma once

#include "geometry/pair.hpp"

namespace principia {
namespace geometry {

template<typename T1, typename T2>
Pair<T1, T2>::Pair(T1 const& t1, T2 const& t2)
    : t1_(t1),
      t2_(t2) {}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator+(
    typename vector_of<Pair>::type const& right) const {
  return Pair<T1, T2>(t1_ + right.t1_, t2_ + right.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator-(
    typename vector_of<Pair>::type const& right) const {
  return Pair<T1, T2>(t1_ - right.t1_, t2_ - right.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator+=(
    typename vector_of<Pair>::type const& right) {
  t1_ += right.t1_;
  t2_ += right.t2_;
  return *this;
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator-=(
    typename vector_of<Pair>::type const& right) {
  t1_ -= right.t1_;
  t2_ -= right.t2_;
  return *this;
}

template<typename T1, typename T2>
bool Pair<T1, T2>::operator==(Pair const& right) const {
  return t1_ == right.t1_ && t2_ == right.t2_;
}

template<typename T1, typename T2>
bool Pair<T1, T2>::operator!=(Pair const& right) const {
  return t1_ != right.t1_ || t2_ != right.t2_;
}

template<typename T1, typename T2>
typename vector_of<Pair<T1, T2>>::type operator-(
    typename enable_if_affine<Pair<T1, T2>>::type const& left,
    Pair<T1, T2> const& right) {
  return typename vector_of<Pair<T1, T2>>::type(left.t1_ - right.t1_,
                                                left.t2_ - right.t2_);
}

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator+(
  Pair<T1, T2> const& right) {
  return right;
}

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator-(
  Pair<T1, T2> const& right) {
  return Pair<T1, T2>(-right.t1_, -right.t2_);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<decltype(std::declval<Scalar>() * std::declval<T1>()),
         decltype(std::declval<Scalar>() * std::declval<T2>())>>::type
operator*(Scalar const left, Pair<T1, T2> const& right) {
  return Pair<decltype(std::declval<Scalar>() * std::declval<T1>()),
              decltype(std::declval<Scalar>() * std::declval<T2>())>(
      left * right.t1_, left * right.t2_);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<decltype(std::declval<T1>() * std::declval<Scalar>()),
         decltype(std::declval<T2>() * std::declval<Scalar>())>>::type
operator*(Pair<T1, T2> const& left, Scalar const right) {
  return Pair<T1, T2>(left.t1_ * right, left.t2_ * right);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<decltype(std::declval<T1>() * std::declval<Scalar>()),
         decltype(std::declval<T2>() * std::declval<Scalar>())>>::type
operator/(Pair<T1, T2> const& left, Scalar const right) {
  return Pair<T1, T2>(left.t1_ / right, left.t2_ / right);
}

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type& operator*=(
    Pair<T1, T2>& left,  // NOLINT(runtime/references)
    double const right) {
  left.t1_ *= right;
  left.t2_ *= right;
  return left;
}

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type& operator/=(
    Pair<T1, T2>& left,  // NOLINT(runtime/references)
    double const right) {
  left.t1_ /= right;
  left.t2_ /= right;
  return left;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair) {
  out << "{" << pair.t1_ << ", " << pair.t2_ << "}";
  return out;
}

}  // namespace geometry
}  // namespace principia
