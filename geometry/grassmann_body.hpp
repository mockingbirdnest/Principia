#pragma once

#include "geometry/grassmann.hpp"

#include <string>

#include "base/not_constructible.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace _grassmann {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::geometry::_sign;
using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename Frame>
Multivector<Scalar, Frame, 1>::Multivector(R3Element<Scalar> const& coordinates)
    : coordinates_(coordinates) {}

template<typename Scalar, typename Frame>
Multivector<Scalar, Frame, 2>::Multivector(R3Element<Scalar> const& coordinates)
    : coordinates_(coordinates) {}

template<typename Scalar, typename Frame>
Multivector<Scalar, Frame, 3>::Multivector(Scalar const& coordinates)
    : coordinates_(coordinates) {}

template<typename Scalar, typename Frame>
R3Element<Scalar> const& Multivector<Scalar, Frame, 1>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
R3Element<Scalar> const& Multivector<Scalar, Frame, 2>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
Scalar const& Multivector<Scalar, Frame, 3>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
Scalar Multivector<Scalar, Frame, 1>::Norm() const {
  return coordinates_.Norm();
}

template<typename Scalar, typename Frame>
Scalar Multivector<Scalar, Frame, 2>::Norm() const {
  return coordinates_.Norm();
}

template<typename Scalar, typename Frame>
Scalar Multivector<Scalar, Frame, 3>::Norm() const {
  // When |Scalar| is double, ADL will not find |Abs|.
  return quantities::Abs(coordinates_);
}

template<typename Scalar, typename Frame>
Square<Scalar> Multivector<Scalar, Frame, 1>::Norm²() const {
  return coordinates_.Norm²();
}

template<typename Scalar, typename Frame>
Square<Scalar> Multivector<Scalar, Frame, 2>::Norm²() const {
  return coordinates_.Norm²();
}

template<typename Scalar, typename Frame>
Square<Scalar> Multivector<Scalar, Frame, 3>::Norm²() const {
  return coordinates_ * coordinates_;
}

template<typename Scalar, typename Frame>
template<typename S>
Multivector<Scalar, Frame, 1>
    Multivector<Scalar, Frame, 1>::OrthogonalizationAgainst(
        Multivector<S, Frame, 1> const& multivector) const {
  return Multivector(
      coordinates_.OrthogonalizationAgainst(multivector.coordinates_));
}

template<typename Scalar, typename Frame>
template<typename S>
Multivector<Scalar, Frame, 2>
    Multivector<Scalar, Frame, 2>::OrthogonalizationAgainst(
        Multivector<S, Frame, 2> const& multivector) const {
  return Multivector(
      coordinates_.OrthogonalizationAgainst(multivector.coordinates_));
}

template<typename Scalar, typename Frame>
void Multivector<Scalar, Frame, 1>::WriteToMessage(
      not_null<serialization::Multivector*> const message) const {
  Frame::WriteToMessage(message->mutable_frame());
  coordinates_.WriteToMessage(message->mutable_vector());
}

template<typename Scalar, typename Frame>
void Multivector<Scalar, Frame, 2>::WriteToMessage(
      not_null<serialization::Multivector*> const message) const {
  Frame::WriteToMessage(message->mutable_frame());
  coordinates_.WriteToMessage(message->mutable_bivector());
}

template<typename Scalar, typename Frame>
void Multivector<Scalar, Frame, 3>::WriteToMessage(
      not_null<serialization::Multivector*> const message) const {
  Frame::WriteToMessage(message->mutable_frame());
  coordinates_.WriteToMessage(message->mutable_trivector());
}

template<typename Scalar, typename Frame>
template<typename, typename>
Multivector<Scalar, Frame, 1> Multivector<Scalar, Frame, 1>::ReadFromMessage(
    serialization::Multivector const& message) {
  Frame::ReadFromMessage(message.frame());
  CHECK(message.has_vector());
  return Multivector(R3Element<Scalar>::ReadFromMessage(message.vector()));
}

template<typename Scalar, typename Frame>
template<typename, typename>
Multivector<Scalar, Frame, 2> Multivector<Scalar, Frame, 2>::ReadFromMessage(
    serialization::Multivector const& message) {
  Frame::ReadFromMessage(message.frame());
  CHECK(message.has_bivector());
  return Multivector(R3Element<Scalar>::ReadFromMessage(message.bivector()));
}

template<typename Scalar, typename Frame>
template<typename, typename>
Multivector<Scalar, Frame, 3> Multivector<Scalar, Frame, 3>::ReadFromMessage(
    serialization::Multivector const& message) {
  Frame::ReadFromMessage(message.frame());
  CHECK(message.has_trivector());
  return Multivector(Scalar::ReadFromMessage(message.trivector()));
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Vector<LScalar, Frame> const& left,
                                       Vector<RScalar, Frame> const& right) {
  return Dot(left.coordinates(), right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Bivector<LScalar, Frame> const& left,
                                       Bivector<RScalar, Frame> const& right) {
  return Dot(left.coordinates(), right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Trivector<LScalar, Frame> const& left,
                                       Trivector<RScalar, Frame> const& right) {
  return left.coordinates() * right.coordinates();
}

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
Trivector<Product<LScalar, RScalar>, Frame> Wedge(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Trivector<Product<LScalar, RScalar>, Frame>(
      Dot(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
Trivector<Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Trivector<Product<LScalar, RScalar>, Frame>(
      Dot(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Commutator(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 1> Normalize(
    Multivector<Scalar, Frame, 1> const& multivector) {
  return Multivector<double, Frame, 1>(Normalize(multivector.coordinates()));
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 2> Normalize(
    Multivector<Scalar, Frame, 2> const& multivector) {
  return Multivector<double, Frame, 2>(Normalize(multivector.coordinates()));
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 3> Normalize(
    Multivector<Scalar, Frame, 3> const& multivector) {
  Scalar const norm = multivector.Norm();
  CHECK_NE(Scalar(), norm);
  return multivector / norm;
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 1> NormalizeOrZero(
    Multivector<Scalar, Frame, 1> const& multivector) {
  return Multivector<double, Frame, 1>(
      NormalizeOrZero(multivector.coordinates()));
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 2> NormalizeOrZero(
    Multivector<Scalar, Frame, 2> const& multivector) {
  return Multivector<double, Frame, 2>(
      NormalizeOrZero(multivector.coordinates()));
}

template<typename Scalar, typename Frame>
Multivector<double, Frame, 3> NormalizeOrZero(
    Multivector<Scalar, Frame, 3> const& multivector) {
  Scalar const norm = multivector.Norm();
  if (norm == Scalar()) {
    static Multivector<double, Frame, 3> zero(0);
    return zero;
  } else {
    return multivector / norm;
  }
}

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

// Implementation from [Kah06], §12 "Mangled Angles", p. 47.
template<typename LScalar, typename RScalar, typename Frame>
Angle AngleBetween(Vector<LScalar, Frame> const& v,
                   Vector<RScalar, Frame> const& w) {
  auto const v_norm_w = v * w.Norm();
  auto const w_norm_v = w * v.Norm();
  return 2 * ArcTan((v_norm_w - w_norm_v).Norm(), (v_norm_w + w_norm_v).Norm());
}

template<typename LScalar, typename RScalar, typename Frame>
Angle AngleBetween(Bivector<LScalar, Frame> const& v,
                   Bivector<RScalar, Frame> const& w) {
  auto const v_norm_w = v * w.Norm();
  auto const w_norm_v = w * v.Norm();
  return 2 * ArcTan((v_norm_w - w_norm_v).Norm(), (v_norm_w + w_norm_v).Norm());
}

template<typename LScalar, typename RScalar, typename PScalar, typename Frame>
Angle OrientedAngleBetween(Vector<LScalar, Frame> const& v,
                           Vector<RScalar, Frame> const& w,
                           Bivector<PScalar, Frame> const& positive) {
  return Sign(InnerProduct(Wedge(v, w), positive)) * AngleBetween(v, w);
}

template<typename LScalar, typename RScalar, typename PScalar, typename Frame>
Angle OrientedAngleBetween(Bivector<LScalar, Frame> const& v,
                           Bivector<RScalar, Frame> const& w,
                           Bivector<PScalar, Frame> const& positive) {
  return Sign(InnerProduct(Commutator(v, w), positive)) * AngleBetween(v, w);
}

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Vector<Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Bivector<Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(+right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(-right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(
      left.coordinates() + right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(
      left.coordinates() - right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank>
operator*(LScalar const& left,
          Multivector<RScalar, Frame, rank> const& right) {
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      left * right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank>
operator*(Multivector<LScalar, Frame, rank> const& left,
          RScalar const& right) {
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      left.coordinates() * right);
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Quotient<LScalar, RScalar>, Frame, rank>
operator/(Multivector<LScalar, Frame, rank> const& left,
          RScalar const& right) {
  return Multivector<Quotient<LScalar, RScalar>, Frame, rank>(
      left.coordinates() / right);
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedMultiplyAdd(
    Multivector<LScalar, Frame, rank> const& a,
    RScalar const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedMultiplyAdd;
  using _r3_element::FusedMultiplyAdd;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedMultiplyAdd(a.coordinates(), b, c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedMultiplySubtract(
    Multivector<LScalar, Frame, rank> const& a,
    RScalar const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedMultiplySubtract;
  using _r3_element::FusedMultiplySubtract;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedMultiplySubtract(a.coordinates(), b, c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedNegatedMultiplyAdd(
    Multivector<LScalar, Frame, rank> const& a,
    RScalar const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedNegatedMultiplyAdd;
  using _r3_element::FusedNegatedMultiplyAdd;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedNegatedMultiplyAdd(a.coordinates(), b, c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank>
FusedNegatedMultiplySubtract(
    Multivector<LScalar, Frame, rank> const& a,
    RScalar const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedNegatedMultiplySubtract;
  using _r3_element::FusedNegatedMultiplySubtract;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedNegatedMultiplySubtract(a.coordinates(), b, c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedMultiplyAdd(
    LScalar const& a,
    Multivector<RScalar, Frame, rank> const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedMultiplyAdd;
  using _r3_element::FusedMultiplyAdd;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedMultiplyAdd(a, b.coordinates(), c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedMultiplySubtract(
    LScalar const& a,
    Multivector<RScalar, Frame, rank> const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedMultiplySubtract;
  using _r3_element::FusedMultiplySubtract;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedMultiplySubtract(a, b.coordinates(), c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank> FusedNegatedMultiplyAdd(
    LScalar const& a,
    Multivector<RScalar, Frame, rank> const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedNegatedMultiplyAdd;
  using _r3_element::FusedNegatedMultiplyAdd;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedNegatedMultiplyAdd(a, b.coordinates(), c.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame, int rank, typename>
Multivector<Product<LScalar, RScalar>, Frame, rank>
FusedNegatedMultiplySubtract(
    LScalar const& a,
    Multivector<RScalar, Frame, rank> const& b,
    Multivector<Product<LScalar, RScalar>, Frame, rank> const& c) {
  using quantities::_elementary_functions::FusedNegatedMultiplySubtract;
  using _r3_element::FusedNegatedMultiplySubtract;
  return Multivector<Product<LScalar, RScalar>, Frame, rank>(
      FusedNegatedMultiplySubtract(a, b.coordinates(), c.coordinates()));
}

template<typename Scalar, typename Frame, int rank>
bool operator==(Multivector<Scalar, Frame, rank> const& left,
                Multivector<Scalar, Frame, rank> const& right) {
  return left.coordinates() == right.coordinates();
}

template<typename Scalar, typename Frame, int rank>
bool operator!=(Multivector<Scalar, Frame, rank> const& left,
                Multivector<Scalar, Frame, rank> const& right) {
  return left.coordinates() != right.coordinates();
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator+=(
    Multivector<Scalar, Frame, rank>& left,
    Multivector<Scalar, Frame, rank> const& right) {
  left.coordinates_ += right.coordinates_;
  return left;
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator-=(
    Multivector<Scalar, Frame, rank>& left,
    Multivector<Scalar, Frame, rank> const& right) {
  left.coordinates_ -= right.coordinates_;
  return left;
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator*=(
    Multivector<Scalar, Frame, rank>& left,
    double const right) {
  left.coordinates_ *= right;
  return left;
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator/=(
    Multivector<Scalar, Frame, rank>& left,
    double const right) {
  left.coordinates_ /= right;
  return left;
}

template<typename Scalar, typename Frame, int rank>
std::string DebugString(Multivector<Scalar, Frame, rank> const& multivector) {
  // This |using| is required for the |Trivector|, whose |DebugString(Scalar)|
  // will not be found by ADL if |Scalar| is |double|.
  using quantities::DebugString;
  return DebugString(multivector.coordinates());
}

template<typename Scalar, typename Frame, int rank>
std::ostream& operator<<(std::ostream& out,
                         Multivector<Scalar, Frame, rank> const& multivector) {
  out << DebugString(multivector);
  return out;
}

}  // namespace internal
}  // namespace _grassmann
}  // namespace geometry
}  // namespace principia
