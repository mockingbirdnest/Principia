#pragma once

#include <string>

#include "geometry/rotation.hpp"

namespace principia {
namespace geometry {

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
inline R3Element<Scalar> const&
Multivector<Scalar, Frame, 1>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
inline R3Element<Scalar> const&
Multivector<Scalar, Frame, 2>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
inline Scalar const& Multivector<Scalar, Frame, 3>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
inline Scalar Multivector<Scalar, Frame, 1>::Norm() const {
  return coordinates_.Norm();
}

template<typename Scalar, typename Frame>
inline Scalar Multivector<Scalar, Frame, 2>::Norm() const {
  return coordinates_.Norm();
}

template<typename Scalar, typename Frame>
inline Scalar Multivector<Scalar, Frame, 3>::Norm() const {
  return quantities::Abs(coordinates_);
}

template<typename LScalar, typename RScalar, typename Frame>
inline quantities::Product<LScalar, RScalar> InnerProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Dot(left.coordinates(), right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
inline quantities::Product<LScalar, RScalar> InnerProduct(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Dot(left.coordinates(), right.coordinates());
}

template<typename LScalar, typename RScalar, typename Frame>
inline quantities::Product<LScalar, RScalar> InnerProduct(
    Trivector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right) {
  return left.coordinates() * right.coordinates();
}

template<typename LScalar, typename RScalar, typename Frame>
inline Bivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Bivector<quantities::Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
inline Trivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Trivector<quantities::Product<LScalar, RScalar>, Frame>(
      Dot(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
inline Trivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Trivector<quantities::Product<LScalar, RScalar>, Frame>(
      Dot(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
inline Bivector<quantities::Product<LScalar, RScalar>, Frame> Commutator(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Bivector<quantities::Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<double, Frame, rank> Normalize(
    Multivector<Scalar, Frame, rank> const& multivector) {
  return multivector / multivector.Norm();
}

template<typename LScalar, typename RScalar, typename Frame>
inline Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Vector<quantities::Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename LScalar, typename RScalar, typename Frame>
inline Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Vector<quantities::Product<LScalar, RScalar>, Frame>(
      Cross(left.coordinates(), right.coordinates()));
}

template<typename Frame>
Rotation<Frame, Frame> Exp(Bivector<quantities::Angle, Frame> const& exponent) {
  return Rotation<Frame, Frame>(exponent.Norm(), exponent);
}

template<typename LScalar, typename RScalar, typename Frame>
inline Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right) {
  return Vector<quantities::Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
inline Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right) {
  return Vector<quantities::Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
inline Bivector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right) {
  return Bivector<quantities::Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}
template<typename LScalar, typename RScalar, typename Frame>
inline Bivector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right) {
  return Bivector<quantities::Product<LScalar, RScalar>, Frame>(
      left.coordinates() * right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(+right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(-right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(
      left.coordinates() + right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(
      left.coordinates() - right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator*(
    double const left,
    Multivector<Scalar, Frame, rank> const& right) {
  return Multivector<Scalar, Frame, rank>(left * right.coordinates());
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator*(
    Multivector<Scalar, Frame, rank> const& left,
    double const right) {
  return Multivector<Scalar, Frame, rank>(left.coordinates() * right);
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank> operator/(
    Multivector<Scalar, Frame, rank> const& left,
    double const right) {
  return Multivector<Scalar, Frame, rank>(left.coordinates() / right);
}

template<typename LDimension, typename RScalar, typename Frame, int rank>
inline Multivector<
    quantities::Product<quantities::Quantity<LDimension>, RScalar>,
    Frame,
    rank>
operator*(quantities::Quantity<LDimension> const& left,
          Multivector<RScalar, Frame, rank> const& right) {
  return Multivector<
      quantities::Product<quantities::Quantity<LDimension>, RScalar>,
      Frame,
      rank>(left * right.coordinates());
}

template<typename LScalar, typename RDimension, typename Frame, int rank>
inline Multivector<
    quantities::Product<LScalar, quantities::Quantity<RDimension>>,
    Frame,
    rank>
operator*(Multivector<LScalar, Frame, rank> const& left,
          quantities::Quantity<RDimension> const& right) {
  return Multivector<
      quantities::Product<LScalar, quantities::Quantity<RDimension>>,
      Frame,
      rank>(left.coordinates() * right);
}

template<typename LScalar, typename RDimension, typename Frame, int rank>
inline Multivector<
    quantities::Quotient<LScalar, quantities::Quantity<RDimension>>,
    Frame,
    rank>
operator/(Multivector<LScalar, Frame, rank> const& left,
          quantities::Quantity<RDimension> const& right) {
  return Multivector<
      quantities::Quotient<LScalar, quantities::Quantity<RDimension>>,
      Frame,
      rank>(left.coordinates() / right);
}

template<typename Scalar, typename Frame, int rank>
inline bool operator==(Multivector<Scalar, Frame, rank> const& left,
                       Multivector<Scalar, Frame, rank> const& right) {
  return left.coordinates() == right.coordinates();
}

template<typename Scalar, typename Frame, int rank>
inline bool operator!=(Multivector<Scalar, Frame, rank> const& left,
                       Multivector<Scalar, Frame, rank> const& right) {
  return left.coordinates() != right.coordinates();
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank>& operator+=(
    Multivector<Scalar, Frame, rank>& left,  // NOLINT(runtime/references)
    Multivector<Scalar, Frame, rank> const& right) {
  return left = left + right;
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank>& operator-=(
    Multivector<Scalar, Frame, rank>& left,  // NOLINT(runtime/references)
    Multivector<Scalar, Frame, rank> const& right) {
  return left = left - right;
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank>& operator*=(
    Multivector<Scalar, Frame, rank>& left,  // NOLINT(runtime/references)
    double const right) {
  return left = left * right;
}

template<typename Scalar, typename Frame, int rank>
inline Multivector<Scalar, Frame, rank>& operator/=(
    Multivector<Scalar, Frame, rank>& left,  // NOLINT(runtime/references)
    double const right) {
  return left = left / right;
}

template<typename Scalar, typename Frame, int rank>
std::string DebugString(Multivector<Scalar, Frame, rank> const& multivector) {
  // This |using| is required for the |Trivector|, since we need an ambiguity
  // between |geometry::DebugString(R3Element<Scalar> const&)| and
  // |quantities::DebugString(Scalar const&)| in order for the template magic
  // to work out.
  using quantities::DebugString;
  return DebugString(multivector.coordinates());
}

template<typename Scalar, typename Frame, int rank>
std::ostream& operator<<(std::ostream& out,
                         Multivector<Scalar, Frame, rank> const& multivector) {
  out << DebugString(multivector);
  return out;
}

}  // namespace geometry
}  // namespace principia
