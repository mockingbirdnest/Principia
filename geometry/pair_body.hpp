#pragma once

#include "geometry/pair.hpp"

#include <string>

#include "geometry/serialization.hpp"

namespace principia {
namespace geometry {
namespace _pair {
namespace internal {

using namespace principia::geometry::_serialization;

template<typename T1, typename T2>
Pair<T1, T2>::Pair(T1 const& t1, T2 const& t2)
    : t1_(t1),
      t2_(t2) {}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator+(
    vector_of_t<Pair<T1, T2>> const& right) const {
  return Pair<T1, T2>(t1_ + right.t1_, t2_ + right.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator-(
    vector_of_t<Pair<T1, T2>> const& right) const {
  return Pair<T1, T2>(t1_ - right.t1_, t2_ - right.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator+=(vector_of_t<Pair<T1, T2>> const& right) {
  t1_ += right.t1_;
  t2_ += right.t2_;
  return *this;
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator-=(vector_of_t<Pair<T1, T2>> const& right) {
  t1_ -= right.t1_;
  t2_ -= right.t2_;
  return *this;
}

template<typename T1, typename T2>
template<typename U1, typename U2>
enable_if_vector_t<Pair<U1, U2>>& Pair<T1, T2>::operator*=(double const right) {
  t1_ *= right;
  t2_ *= right;
  return *this;
}

template<typename T1, typename T2>
template<typename U1, typename U2>
enable_if_vector_t<Pair<U1, U2>>& Pair<T1, T2>::operator/=(double const right) {
  t1_ /= right;
  t2_ /= right;
  return *this;
}

template<typename T1, typename T2>
T1 const& Pair<T1, T2>::position() const
  requires is_position_v<T1> {
  return t1_;
}

template<typename T1, typename T2>
T1 const& Pair<T1, T2>::displacement() const
  requires is_displacement_v<T1> {
  return t1_;
}

template<typename T1, typename T2>
T2 const& Pair<T1, T2>::velocity() const
  requires is_velocity_v<T2> {
  return t2_;
}

template<typename T1, typename T2>
void Pair<T1, T2>::WriteToMessage(
    not_null<serialization::Pair*> const message) const {
  PointOrMultivectorSerializer<T1, serialization::Pair::Element>::
      WriteToMessage(t1_, message->mutable_t1());
  PointOrMultivectorSerializer<T2, serialization::Pair::Element>::
      WriteToMessage(t2_, message->mutable_t2());
}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::ReadFromMessage(serialization::Pair const& message)
  requires serializable<T1> && serializable<T2> {
  T1 const t1 = PointOrMultivectorSerializer<T1, serialization::Pair::Element>::
                    ReadFromMessage(message.t1());
  T2 const t2 = PointOrMultivectorSerializer<T2, serialization::Pair::Element>::
                    ReadFromMessage(message.t2());
  return {t1, t2};
}

template<typename T1, typename T2>
vector_of_t<Pair<T1, T2>> operator-(
    enable_if_affine_t<Pair<T1, T2>> const& left,
    Pair<T1, T2> const& right) {
  return vector_of_t<Pair<T1, T2>>(left.t1_ - right.t1_,
                                                left.t2_ - right.t2_);
}

template<typename T1, typename T2>
enable_if_vector_t<Pair<T1, T2>> operator+(Pair<T1, T2> const& right) {
  return right;
}

template<typename T1, typename T2>
enable_if_vector_t<Pair<T1, T2>> operator-(Pair<T1, T2> const& right) {
  return Pair<T1, T2>(-right.t1_, -right.t2_);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Product<Scalar, T1>, Product<Scalar, T2>>>::type
operator*(Scalar const left, Pair<T1, T2> const& right) {
  return Pair<Product<Scalar, T1>, Product<Scalar, T2>>(left * right.t1_,
                                                        left * right.t2_);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Product<T1, Scalar>, Product<T2, Scalar>>>::type
operator*(Pair<T1, T2> const& left, Scalar const right) {
  return Pair<Product<T1, Scalar>, Product<T2, Scalar>>(left.t1_ * right,
                                                        left.t2_ * right);
}

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Quotient<T1, Scalar>, Quotient<T2, Scalar>>>::type
operator/(Pair<T1, T2> const& left, Scalar const right) {
  return Pair<Quotient<T1, Scalar>, Quotient<T2, Scalar>>(left.t1_ / right,
                                                          left.t2_ / right);
}

template<typename T1, typename T2>
std::string DebugString(Pair<T1, T2> const& pair) {
  return "{" + DebugString(pair.t1_) + ", " + DebugString(pair.t2_) + "}";
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair) {
  out << "{" << pair.t1_ << ", " << pair.t2_ << "}";
  return out;
}

}  // namespace internal
}  // namespace _pair
}  // namespace geometry

namespace base {
namespace _mappable {
namespace internal {

template<typename Functor, typename T1, typename T2>
typename Mappable<Functor,
                  Pair<T1, T2>,
                  enable_if_vector_t<Pair<T1, T2>, void>>::type
Mappable<Functor, Pair<T1, T2>, enable_if_vector_t<Pair<T1, T2>, void>>::Do(
    Functor const& functor,
    Pair<T1, T2> const& pair) {
  return type(functor(pair.t1_), functor(pair.t2_));
}

}  // namespace internal
}  // namespace _mappable
}  // namespace base

}  // namespace principia
