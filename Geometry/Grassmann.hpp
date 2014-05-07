#pragma once

#include "Geometry/R3Element.hpp"
#include "Quantities/Dimensionless.hpp"
#include "Quantities/Quantities.hpp"

namespace principia {
namespace geometry {

// A multivector of rank Rank on a three-dimensional real inner product space
// bearing the dimensionality of Scalar.
// Do not use this type for Rank = 0 (scalar), use the Scalar type directly
// instead.
template<typename Scalar, typename Frame, unsigned int Rank>
struct Multivector;

template<typename Scalar, typename Frame>
class Multivector<Scalar, Frame, 1> {
 public:
  explicit Multivector(R3Element<Scalar> const& coordinates);
  ~Multivector() = default;

  R3Element<Scalar> coordinates() const;

 private:
  R3Element<Scalar> coordinates_;
};

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 2> {
  explicit Multivector(R3Element<Scalar> const& coordinates);
  ~Multivector() = default;

  R3Element<Scalar> coordinates() const;

 private:
  R3Element<Scalar> coordinates_;
};

template<typename Scalar, typename Frame>
struct Multivector<Scalar, Frame, 3> {
  explicit Multivector(Scalar const& coordinates);
  ~Multivector() = default;

  Scalar coordinates() const;

 private:
  Scalar coordinates_;
};

template<typename Scalar, typename Frame>
using Vector = Multivector<Scalar, Frame, 1>;
template<typename Scalar, typename Frame>
using Bivector = Multivector<Scalar, Frame, 2>;
template<typename Scalar, typename Frame>
using Trivector = Multivector<Scalar, Frame, 3>;

template<typename LScalar, typename RScalar, typename Frame>
quantities::Product<LScalar, RScalar> InnerProduct(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
quantities::Product<LScalar, RScalar> InnerProduct(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
quantities::Product<LScalar, RScalar> InnerProduct(
    Trivector<LScalar, Frame> const& left,
    Trivector<RScalar, Frame> const& right);

template<typename LScalar, typename RScalar, typename Frame>
Bivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Trivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);
template<typename LScalar, typename RScalar, typename Frame>
Trivector<quantities::Product<LScalar, RScalar>, Frame> Wedge(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

// Lie bracket on V ^ V = so(V).
template<typename LScalar, typename RScalar, typename Frame>
Bivector<quantities::Product<LScalar, RScalar>, Frame> Commutator(
    Bivector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

// Left action of V ^ V = so(V) on V.
template<typename LScalar, typename RScalar, typename Frame>
Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Bivector<LScalar, Frame> const& left,
    Vector<RScalar, Frame> const& right);

// Right action of V ^ V = so(V) on V* = V.
template<typename LScalar, typename RScalar, typename Frame>
Vector<quantities::Product<LScalar, RScalar>, Frame> operator*(
    Vector<LScalar, Frame> const& left,
    Bivector<RScalar, Frame> const& right);

template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator+(
    Multivector<Scalar, Frame, Rank> const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator-(
    Multivector<Scalar, Frame, Rank> const& right);

template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator+(
    Multivector<Scalar, Frame, Rank> const& left,
    Multivector<Scalar, Frame, Rank> const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator-(
    Multivector<Scalar, Frame, Rank> const& left,
    Multivector<Scalar, Frame, Rank> const& right);

template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator*(
    quantities::Dimensionless const& left,
    Multivector<Scalar, Frame, Rank> const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator*(
    Multivector<Scalar, Frame, Rank> const& left,
    quantities::Dimensionless const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
Multivector<Scalar, Frame, Rank> operator/(
    Multivector<Scalar, Frame, Rank> const& left,
    quantities::Dimensionless const& right);

template<typename LDimension, typename RScalar, typename Frame,
         unsigned int Rank>
Multivector<quantities::Product<quantities::Quantity<LDimension>, RScalar>,
            Frame, Rank>
operator*(quantities::Quantity<LDimension> const& left,
          Multivector<RScalar, Frame, Rank> const& right);

template<typename LScalar, typename RDimension, typename Frame,
         unsigned int Rank>
Multivector<quantities::Product<LScalar, quantities::Quantity<RDimension>>,
            Frame, Rank>
operator*(Multivector<LScalar, Frame, Rank> const& left,
          quantities::Quantity<RDimension> const& right);

template<typename LScalar, typename RDimension, typename Frame,
         unsigned int Rank>
Multivector<quantities::Quotient<LScalar, quantities::Quantity<RDimension>>,
            Frame, Rank>
operator/(Multivector<LScalar, Frame, Rank> const& left,
          quantities::Quantity<RDimension> const& right);

template<typename Scalar, typename Frame, unsigned int Rank>
void operator+=(Multivector<Scalar, Frame, Rank>& left,
                Multivector<Scalar, Frame, Rank> const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
void operator-=(Multivector<Scalar, Frame, Rank>& left,
                Multivector<Scalar, Frame, Rank> const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
void operator*=(Multivector<Scalar, Frame, Rank>& left,
                quantities::Dimensionless const& right);
template<typename Scalar, typename Frame, unsigned int Rank>
void operator/=(Multivector<Scalar, Frame, Rank>& left,
                quantities::Dimensionless const& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/Grassmann-body.hpp"
