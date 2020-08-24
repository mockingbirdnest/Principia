
#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <utility>

#include "base/not_null.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
FORWARD_DECLARE_FROM(polynomial,
                     TEMPLATE(typename Value, typename Argument, int degree_,
                              template<typename, typename, int> class Evaluator)
                              class,
                     PolynomialInMonomialBasis);
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value, typename Argument, int degree_,
             template<typename, typename, int> class Evaluator,
             typename OptionalExpressIn = std::nullopt_t) std::string,
    ToMathematica,
    (numerics::
         PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
             polynomial,
     OptionalExpressIn express_in = std::nullopt));
}  // namespace mathematica

namespace numerics {
namespace internal_polynomial {

using base::not_constructible;
using base::not_null;
using geometry::Point;
using quantities::Derivative;
using quantities::Derivatives;
using quantities::Difference;
using quantities::Primitive;
using quantities::Product;
using quantities::Quotient;

// |Value| must belong to an affine space.  |Argument| must belong to a ring or
// to Point based on a ring.
// TODO(phl): We would like the base case to be the affine case (not limited to
// Point) and the specialized case to check for the existence of Sum and Product
// for Argument, and that works with Clang but not with VS2015.  Revisit once
// MSFT has fixed their bugs.
template<typename Value, typename Argument>
class Polynomial {
 public:
  // This virtual destructor makes this class and its subclasses non-literal, so
  // constexpr-ness is a bit of a lie for polynomials.
  // TODO(phl): Consider providing an explicit deleter function that would allow
  // making the destructor protected and nonvirtual.
  virtual ~Polynomial() = default;

  virtual Value Evaluate(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

  // Only useful for benchmarking or analyzing performance.  Do not use in real
  // code.
  //TODO(phl):Fix comment.
  virtual int degree() const = 0;

  // Only useful for logging.  Do not use in real code.
  virtual bool is_zero() const = 0;

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
  using Coefficients = Derivatives<Value, Argument, degree_ + 1>;

  // The coefficients are applied to powers of argument.
  explicit constexpr PolynomialInMonomialBasis(
      Coefficients coefficients);

  // A polynomial may be explicitly converted to a higher degree (possibly with
  // a different evaluator).
  template<int higher_degree_,
           template<typename, typename, int> class HigherEvaluator>
  explicit operator PolynomialInMonomialBasis<
      Value, Argument, higher_degree_, HigherEvaluator>() const;

  FORCE_INLINE(inline) Value
  Evaluate(Argument const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Argument const& argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  template<int order = 1>
  PolynomialInMonomialBasis<
      Derivative<Value, Argument, order>, Argument, degree_ - order, Evaluator>
  Derivative() const;

  // The constant term of the result is zero.
  PolynomialInMonomialBasis<
      Primitive<Difference<Value>, Argument>, Argument, degree_ + 1, Evaluator>
  Primitive() const;

  PolynomialInMonomialBasis& operator+=(PolynomialInMonomialBasis const& right);
  PolynomialInMonomialBasis& operator-=(PolynomialInMonomialBasis const& right);

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;

  template<typename V, typename A, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator-(PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E>
  friend operator+(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<S, V>, A, d, E>
  friend operator*(S const& left,
                   PolynomialInMonomialBasis<V, A, d, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<V, S>, A, d, E>
  friend operator*(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Quotient<V, S>, A, d, E>
  friend operator/(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<L, R>, A, l + r, E>
  friend operator*(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename V, typename A, int d,
           template<typename, typename, int> class E>
  friend std::ostream& operator<<(
      std::ostream& out,
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial);
  template<typename V, typename A, int d,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematica(
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
      O express_in);
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
  using Coefficients = Derivatives<Value, Argument, degree_ + 1>;

  // The coefficients are relative to origin; in other words they are applied to
  // powers of (argument - origin).
  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      Point<Argument> const& origin);

  // A polynomial may be explicitly converted to a higher degree (possibly with
  // a different evaluator).
  template<int higher_degree_,
           template<typename, typename, int> class HigherEvaluator>
  explicit operator PolynomialInMonomialBasis<
      Value, Point<Argument>, higher_degree_, HigherEvaluator>() const;

  FORCE_INLINE(inline) Value
  Evaluate(Point<Argument> const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Point<Argument> const& argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  Point<Argument> const& origin() const;

  // Returns a copy of this polynomial adjusted to the given origin.
  PolynomialInMonomialBasis AtOrigin(Point<Argument> const& origin) const;

  template<int order = 1>
  PolynomialInMonomialBasis<
      Derivative<Value, Argument, order>, Point<Argument>, degree_ - order,
      Evaluator>
  Derivative() const;

  // The constant term of the result is zero.
  PolynomialInMonomialBasis<
      Primitive<Difference<Value>, Argument>, Point<Argument>,
      degree_ + 1, Evaluator>
  Primitive() const;

  PolynomialInMonomialBasis& operator+=(const PolynomialInMonomialBasis& right);
  PolynomialInMonomialBasis& operator-=(const PolynomialInMonomialBasis& right);

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;
  Point<Argument> origin_;

  template<typename V, typename A, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator-(PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E>
  friend operator+(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<S, V>, A, d, E>
  friend operator*(S const& left,
                   PolynomialInMonomialBasis<V, A, d, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<V, S>, A, d, E>
  friend operator*(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Quotient<V, S>, A, d, E>
  friend operator/(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> class E>
  constexpr PolynomialInMonomialBasis<Product<L, R>, A, l + r, E>
  friend operator*(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename V, typename A, int d,
           template<typename, typename, int> class E>
  friend std::ostream& operator<<(
      std::ostream& out,
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial);
  template<typename V, typename A, int d,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::internal_mathematica::ToMathematica(
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
      O express_in);
};

// Vector space of polynomials.

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator+(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right);

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right);

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator>
operator+(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right);

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator>
operator-(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Product<Scalar, Value>, Argument,
                                    degree_, Evaluator>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Product<Value, Scalar>, Argument,
                                    degree_, Evaluator>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument,
                                    degree_, Evaluator>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right);

// Algebra of polynomials.

template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                    ldegree_ + rdegree_, Evaluator>
operator*(
    PolynomialInMonomialBasis<LValue, Argument, ldegree_, Evaluator> const&
        left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_, Evaluator> const&
        right);

// Output.

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
        polynomial);

}  // namespace internal_polynomial

using internal_polynomial::Polynomial;
using internal_polynomial::PolynomialInMonomialBasis;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
