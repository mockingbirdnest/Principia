#pragma once

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
inline R3Element<Scalar> Multivector<Scalar, Frame, 1>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
inline R3Element<Scalar> Multivector<Scalar, Frame, 2>::coordinates() const {
  return coordinates_;
}

template<typename Scalar, typename Frame>
inline Scalar Multivector<Scalar, Frame, 3>::coordinates() const {
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

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator+(
    Multivector<Scalar, Frame, Rank> const& right) {
  return Multivector<Scalar, Frame, Rank>(+right.coordinates());
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator-(
    Multivector<Scalar, Frame, Rank> const& right) {
  return Multivector<Scalar, Frame, Rank>(-right.coordinates());
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator+(
    Multivector<Scalar, Frame, Rank> const& left,
    Multivector<Scalar, Frame, Rank> const& right) {
  return Multivector<Scalar, Frame, Rank>(
      left.coordinates() + right.coordinates());
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator-(
    Multivector<Scalar, Frame, Rank> const& left,
    Multivector<Scalar, Frame, Rank> const& right) {
  return Multivector<Scalar, Frame, Rank>(
      left.coordinates() - right.coordinates());
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator*(
    double const left,
    Multivector<Scalar, Frame, Rank> const& right) {
  return Multivector<Scalar, Frame, Rank>(left * right.coordinates());
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator*(
    Multivector<Scalar, Frame, Rank> const& left,
    double const right) {
  return Multivector<Scalar, Frame, Rank>(left.coordinates() * right);
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank> operator/(
    Multivector<Scalar, Frame, Rank> const& left,
    double const right) {
  return Multivector<Scalar, Frame, Rank>(left.coordinates() / right);
}

template<typename LDimension, typename RScalar, typename Frame,
         unsigned int Rank>
inline Multivector<
    quantities::Product<quantities::Quantity<LDimension>, RScalar>,
    Frame,
    Rank>
operator*(quantities::Quantity<LDimension> const& left,
          Multivector<RScalar, Frame, Rank> const& right) {
  return Multivector<
      quantities::Product<quantities::Quantity<LDimension>, RScalar>,
      Frame,
      Rank>(left * right.coordinates());
}

template<typename LScalar, typename RDimension, typename Frame,
         unsigned int Rank>
inline Multivector<
    quantities::Product<LScalar, quantities::Quantity<RDimension>>,
    Frame,
    Rank>
operator*(Multivector<LScalar, Frame, Rank> const& left,
          quantities::Quantity<RDimension> const& right) {
  return Multivector<
      quantities::Product<LScalar, quantities::Quantity<RDimension>>,
      Frame,
      Rank>(left.coordinates() * right);
}

template<typename LScalar, typename RDimension, typename Frame,
         unsigned int Rank>
inline Multivector<
    quantities::Quotient<LScalar, quantities::Quantity<RDimension>>,
    Frame,
    Rank>
operator/(Multivector<LScalar, Frame, Rank> const& left,
          quantities::Quantity<RDimension> const& right) {
  return Multivector<
      quantities::Quotient<LScalar, quantities::Quantity<RDimension>>,
      Frame,
      Rank>(left.coordinates() / right);
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline bool operator==(Multivector<Scalar, Frame, Rank> const& left,
                       Multivector<Scalar, Frame, Rank> const& right) {
  return left.coordinates() == right.coordinates();
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline bool operator!=(Multivector<Scalar, Frame, Rank> const& left,
                       Multivector<Scalar, Frame, Rank> const& right) {
  return left.coordinates() != right.coordinates();
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank>& operator+=(
    Multivector<Scalar, Frame, Rank>& left,  // NOLINT(runtime/references)
    Multivector<Scalar, Frame, Rank> const& right) {
  return left = left + right;
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank>& operator-=(
    Multivector<Scalar, Frame, Rank>& left,  // NOLINT(runtime/references)
    Multivector<Scalar, Frame, Rank> const& right) {
  return left = left - right;
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank>& operator*=(
    Multivector<Scalar, Frame, Rank>& left,  // NOLINT(runtime/references)
    double const right) {
  return left = left * right;
}

template<typename Scalar, typename Frame, unsigned int Rank>
inline Multivector<Scalar, Frame, Rank>& operator/=(
    Multivector<Scalar, Frame, Rank>& left,  // NOLINT(runtime/references)
    double const right) {
  return left = left / right;
}

template<typename Scalar, typename Frame, unsigned int Rank>
std::ostream& operator<<(std::ostream& out,
                         Multivector<Scalar, Frame, Rank> const& multivector) {
  out << multivector.coordinates();
  return out;
}

}  // namespace geometry
}  // namespace principia
