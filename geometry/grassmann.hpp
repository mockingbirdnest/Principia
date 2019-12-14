
#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>

#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/quantities.hpp"
#include "quantities/traits.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

FORWARD_DECLARE_FROM(rotation,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Rotation);

namespace internal_grassmann {

using base::not_null;
using quantities::Angle;
using quantities::is_quantity;
using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Square;

// A multivector of rank |rank| on a three-dimensional real inner product
// space bearing the dimensionality of |Scalar|, i.e., an element of
// ⋀ⁿ Scalar³. Do not use this type for |rank == 0| (scalars), use the |Scalar|
// type directly instead.
// |Frame| represents a reference frame together with an orthonormal basis.
template<typename Scalar, typename Frame, int rank>
class Multivector;

template<typename Scalar, typename Frame>
class Multivector<Scalar, Frame, 1> final {
 public:
  constexpr Multivector();
  explicit Multivector(R3Element<Scalar> const& coordinates);

  R3Element<Scalar> const& coordinates() const;
  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  template<typename S>
  Multivector OrthogonalizationAgainst(
      Multivector<S, Frame, 1> const& multivector) const;

  void WriteToMessage(
      not_null<serialization::Multivector*> message) const;
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static Multivector ReadFromMessage(serialization::Multivector const& message);

 private:
  R3Element<Scalar> coordinates_;

  template<typename S, typename F, int r>
  friend class Multivector;

  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator+=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator-=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator*=(Multivector<S, F, r>& left,
                                          double right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator/=(Multivector<S, F, r>& left,
                                          double right);
};

template<typename Scalar, typename Frame>
class Multivector<Scalar, Frame, 2> final {
 public:
  constexpr Multivector();
  explicit Multivector(R3Element<Scalar> const& coordinates);

  R3Element<Scalar> const& coordinates() const;
  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  template<typename S>
  Multivector OrthogonalizationAgainst(
      Multivector<S, Frame, 2> const& multivector) const;

  void WriteToMessage(not_null<serialization::Multivector*> message) const;
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static Multivector ReadFromMessage(serialization::Multivector const& message);

 private:
  R3Element<Scalar> coordinates_;

  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator+=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator-=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator*=(Multivector<S, F, r>& left,
                                          double right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator/=(Multivector<S, F, r>& left,
                                          double right);
};

template<typename Scalar, typename Frame>
class Multivector<Scalar, Frame, 3> final {
 public:
  constexpr Multivector();
  explicit Multivector(Scalar const& coordinates);

  Scalar const& coordinates() const;
  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  void WriteToMessage(not_null<serialization::Multivector*> message) const;
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  static Multivector ReadFromMessage(serialization::Multivector const& message);

 private:
  Scalar coordinates_;

  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator+=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator-=(Multivector<S, F, r>& left,
                                          Multivector<S, F, r> const& right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator*=(Multivector<S, F, r>& left,
                                          double right);
  template<typename S, typename F, int r>
  friend Multivector<S, F, r>& operator/=(Multivector<S, F, r>& left,
                                          double right);
};

template<typename Scalar, typename Frame>
using Vector = Multivector<Scalar, Frame, 1>;
template<typename Scalar, typename Frame>
using Bivector = Multivector<Scalar, Frame, 2>;
template<typename Scalar, typename Frame>
using Trivector = Multivector<Scalar, Frame, 3>;

template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Vector<LScalar, Frame> const& left,
                                       Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Bivector<LScalar, Frame> const& left,
                                       Bivector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Product<LScalar, RScalar> InnerProduct(Trivector<LScalar, Frame> const& left,
                                       Trivector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Trivector<Product<LScalar, RScalar>, Frame> Wedge(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Trivector<Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

// Lie bracket on 𝑉 ∧ 𝑉 ≅ 𝖘𝔬(𝑉).
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> Commutator(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

// Returns multivector / ‖multivector‖.
template<typename Scalar, typename Frame, int rank>
Multivector<double, Frame, rank> Normalize(
    Multivector<Scalar, Frame, rank> const& multivector);

// Same as above, but returns the zero multivector if its argument is zero.
template<typename Scalar, typename Frame, int rank>
Multivector<double, Frame, rank> NormalizeOrZero(
    Multivector<Scalar, Frame, rank> const& multivector);

// Left action of 𝑉 ∧ 𝑉 ≅ 𝖘𝔬(𝑉) on 𝑉.
template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

// Right action of 𝑉 ∧ 𝑉 ≅ 𝖘𝔬(𝑉) on 𝑉* ≅ 𝑉.
template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

// Exponential map 𝑉 ∧ 𝑉 ≅ 𝖘𝔬(𝑉) → SO(𝑉).
template<typename Frame>
Rotation<Frame, Frame> Exp(Bivector<Angle, Frame> const& exponent);

// The result is in [0, π]; the function is commutative.
template<typename LScalar, typename RScalar, typename Frame>
Angle AngleBetween(Vector<LScalar, Frame> const& v,
                   Vector<RScalar, Frame> const& w);

// The result is in [0, π]; the function is commutative.
template<typename LScalar, typename RScalar, typename Frame>
Angle AngleBetween(Bivector<LScalar, Frame> const& v,
                   Bivector<RScalar, Frame> const& w);

// The result is in [-π, π]; the function is anticommutative, the result is in
// [0, π] if |InnerProduct(Wedge(v, w), positive) >= 0|.
template<typename LScalar, typename RScalar, typename PScalar, typename Frame>
Angle OrientedAngleBetween(Vector<LScalar, Frame> const& v,
                           Vector<RScalar, Frame> const& w,
                           Bivector<PScalar, Frame> const& positive);

// The result is in [-π, π]; the function is anticommutative, the result is in
// [0, π] if |InnerProduct(Commutator(v, w), positive) >= 0|.
template<typename LScalar, typename RScalar, typename PScalar, typename Frame>
Angle OrientedAngleBetween(Bivector<LScalar, Frame> const& v,
                           Bivector<RScalar, Frame> const& w,
                           Bivector<PScalar, Frame> const& positive);

template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Vector<Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Bivector<Product<LScalar, RScalar>, Frame> operator*(
    Trivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& right);
template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& right);

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator+(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right);
template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank> operator-(
    Multivector<Scalar, Frame, rank> const& left,
    Multivector<Scalar, Frame, rank> const& right);

template<typename LScalar, typename RScalar, typename Frame, int rank,
         typename = std::enable_if_t<is_quantity<LScalar>::value>>
Multivector<Product<LScalar, RScalar>, Frame, rank>
operator*(LScalar const& left,
          Multivector<RScalar, Frame, rank> const& right);

template<typename LScalar, typename RScalar, typename Frame, int rank,
         typename = std::enable_if_t<is_quantity<RScalar>::value>>
Multivector<Product<LScalar, RScalar>, Frame, rank>
operator*(Multivector<LScalar, Frame, rank> const& left,
          RScalar const& right);

template<typename LScalar, typename RScalar, typename Frame, int rank,
         typename = std::enable_if_t<is_quantity<RScalar>::value>>
Multivector<Quotient<LScalar, RScalar>, Frame, rank>
operator/(Multivector<LScalar, Frame, rank> const& left,
          RScalar const& right);

template<typename Scalar, typename Frame, int rank>
bool operator==(Multivector<Scalar, Frame, rank> const& left,
                Multivector<Scalar, Frame, rank> const& right);
template<typename Scalar, typename Frame, int rank>
bool operator!=(Multivector<Scalar, Frame, rank> const& left,
                Multivector<Scalar, Frame, rank> const& right);

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator+=(
    Multivector<Scalar, Frame, rank>& left,
    Multivector<Scalar, Frame, rank> const& right);
template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator-=(
    Multivector<Scalar, Frame, rank>& left,
    Multivector<Scalar, Frame, rank> const& right);
template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator*=(
    Multivector<Scalar, Frame, rank>& left,
    double right);
template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>& operator/=(
    Multivector<Scalar, Frame, rank>& left,
    double right);

template<typename Scalar, typename Frame, int rank>
std::string DebugString(Multivector<Scalar, Frame, rank> const& multivector);

template<typename Scalar, typename Frame, int rank>
std::ostream& operator<<(std::ostream& out,
                         Multivector<Scalar, Frame, rank> const& multivector);

}  // namespace internal_grassmann

using internal_grassmann::AngleBetween;
using internal_grassmann::Bivector;
using internal_grassmann::Commutator;
using internal_grassmann::Exp;
using internal_grassmann::InnerProduct;
using internal_grassmann::Multivector;
using internal_grassmann::Normalize;
using internal_grassmann::NormalizeOrZero;
using internal_grassmann::OrientedAngleBetween;
using internal_grassmann::Trivector;
using internal_grassmann::Vector;
using internal_grassmann::Wedge;

}  // namespace geometry
}  // namespace principia

#include "geometry/grassmann_body.hpp"
