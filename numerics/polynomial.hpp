
#pragma once

#include <tuple>
#include <utility>

#include "base/not_null.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using base::not_null;
using geometry::Point;
using quantities::Derivative;
using quantities::NthDerivative;

template<typename Value, typename Argument, typename>
struct NthDerivativesGenerator;
template<typename Value, typename Argument, int... orders>
struct NthDerivativesGenerator<Value,
                               Argument,
                               std::integer_sequence<int, orders...>> {
  using Type = std::tuple<NthDerivative<Value, Argument, orders>...>;
};

template<typename Value, typename Argument, typename Sequence>
using NthDerivatives =
    typename NthDerivativesGenerator<Value, Argument, Sequence>::Type;

// |Value| must belong to an affine space.  |Argument| must belong to a ring or
// to Point based on a ring.
// TODO(phl): We would like the base case to be the affine case (not limited to
// Point) and the specialized case to check for the existence of Sum and Product
// for Argument, and that works with Clang but not with VS2015.  Revisit once
// MSFT has fixed their bugs.
template<typename Value, typename Argument>
class Polynomial {
 public:
  virtual ~Polynomial() = default;

  virtual Value Evaluate(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

  // Only useful for benchmarking or analyzing performance.  Do not use in real
  // code.
  virtual int degree() const = 0;

  virtual void WriteToMessage(
      not_null<serialization::Polynomial*> message) const = 0;

  // The evaluator is not part of the serialization because it's fine to read
  // with a different evaluator than the one the polynomial was written with.
  template<template<typename, typename, int> class Evaluator>
  static not_null<std::unique_ptr<Polynomial>> ReadFromMessage(
      serialization::Polynomial const& message);
};

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
class PolynomialInMonomialBasis : public Polynomial<Value, Argument> {
 public:
  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients =
      NthDerivatives<Value,
                     Argument,
                     std::make_integer_sequence<int, degree_ + 1>>;

  // The coefficients are applied to powers of argument.
  explicit PolynomialInMonomialBasis(Coefficients const& coefficients);

  FORCE_INLINE(inline) Value
  Evaluate(Argument const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Argument const& argument) const override;

  constexpr int degree() const override;

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;
};

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
class PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>
    : public Polynomial<Value, Point<Argument>> {
 public:
  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients =
      NthDerivatives<Value,
                     Argument,
                     std::make_integer_sequence<int, degree_ + 1>>;

  // The coefficients are relative to origin; in other words they are applied to
  // powers of (argument - origin).
  PolynomialInMonomialBasis(Coefficients const& coefficients,
                            Point<Argument> const& origin);

  FORCE_INLINE(inline) Value
  Evaluate(Point<Argument> const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Point<Argument> const& argument) const override;

  constexpr int degree() const override;

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;
  Point<Argument> origin_;
};

}  // namespace internal_polynomial

using internal_polynomial::NthDerivative;
using internal_polynomial::Polynomial;
using internal_polynomial::PolynomialInMonomialBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
