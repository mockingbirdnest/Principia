// The files containing the tree of child classes of `Polynomial` must be
// included in the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#include "numerics/polynomial.hpp"
#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_IN_MONOMIAL_BASIS_HPP_
#define PRINCIPIA_NUMERICS_POLYNOMIAL_IN_MONOMIAL_BASIS_HPP_

#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "base/algebra.hpp"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/point.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

// The presence of an operator+ below causes a bizarre compilation error in
// seemingly unrelated code in PolynomialTestInMonomialBasis.VectorSpace.
// In recent versions of Visual Studio, this problem has been resolved; however,
// adding constraints to Polynomial brings it back unless the relevant operators
// are constrained.
#define PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS \
  !PRINCIPIA_COMPILER_MSVC || !(_MSC_FULL_VER == 192'930'036 || \
                                _MSC_FULL_VER == 192'930'037 || \
                                _MSC_FULL_VER == 192'930'038 || \
                                _MSC_FULL_VER == 192'930'133 || \
                                _MSC_FULL_VER == 192'930'139 || \
                                _MSC_FULL_VER == 192'930'143 || \
                                _MSC_FULL_VER == 192'930'147 || \
                                _MSC_FULL_VER == 193'431'937 || \
                                _MSC_FULL_VER == 193'431'942 || \
                                _MSC_FULL_VER == 193'431'944 || \
                                _MSC_FULL_VER == 193'532'216 || \
                                _MSC_FULL_VER == 193'532'217 || \
                                _MSC_FULL_VER == 193'632'532 || \
                                _MSC_FULL_VER == 193'632'535 || \
                                _MSC_FULL_VER == 194'234'435 || \
                                _MSC_FULL_VER == 194'435'213)

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value, typename Argument, int degree,
             TEMPLATE(typename, typename, int) typename Evaluator) class,
    PolynomialInMonomialBasis,
    FROM(polynomial_in_monomial_basis));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value, typename Argument, int degree,
             TEMPLATE(typename, typename, int) typename Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_polynomial_in_monomial_basis::
        PolynomialInMonomialBasis<Value, Argument, degree, Evaluator> const&
            polynomial,
    OptionalExpressIn express_in),
    FROM(mathematica));
}  // namespace mathematica

namespace numerics {
namespace _polynomial_in_monomial_basis {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_point;
using namespace principia::numerics::_polynomial;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::quantities::_tuples;

// Used to decide which evaluator to use for a particular polynomial.
class Policy {
 public:
  template<typename Value, typename Argument, int degree,
           template<typename, typename, int> typename Evaluator_>
  not_null<std::unique_ptr<Polynomial<Value, Argument>>>
  WithEvaluator(PolynomialInMonomialBasis<Value, Argument, degree, Evaluator_>&&
                    polynomial) const;

  static constexpr Policy AlwaysEstrin();
  static constexpr Policy AlwaysEstrinWithoutFMA();

  void WriteToMessage(
      not_null<serialization::PolynomialInMonomialBasis::Policy*> message)
      const;
  static Policy ReadFromMessage(
      serialization::PolynomialInMonomialBasis::Policy const& message);

 private:
  explicit constexpr Policy(
      serialization::PolynomialInMonomialBasis::Policy::Kind kind);

  serialization::PolynomialInMonomialBasis::Policy::Kind kind_;
};

template<typename Value, typename Argument, int degree>
using DefaultEvaluator = std::conditional_t<degree <= 3,  // NOLINT
                                            Horner<Value, Argument, degree>,
                                            Estrin<Value, Argument, degree>>;


template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator_ =
             DefaultEvaluator>
class PolynomialInMonomialBasis : public Polynomial<Value_, Argument_> {
 public:
  using Argument = Argument_;
  using Value = Value_;
  template<typename V, typename A, int d>
  using Evaluator = Evaluator_<V, A, d>;

  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients = Derivatives<Value, Argument, degree_ + 1>;

  // The coefficients are relative to origin; in other words they are applied to
  // powers of (argument - origin).

  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      Argument const& origin);
  constexpr explicit PolynomialInMonomialBasis(Coefficients coefficients)
    requires additive_group<Argument>;
  constexpr PolynomialInMonomialBasis() requires additive_group<Argument>;

  friend constexpr bool operator==(PolynomialInMonomialBasis const& left,
                                   PolynomialInMonomialBasis const& right) =
      default;
  friend constexpr bool operator!=(PolynomialInMonomialBasis const& left,
                                   PolynomialInMonomialBasis const& right) =
      default;

  // A polynomial may be explicitly converted to a higher degree.
  template<int higher_degree_, typename Self>
  explicit operator PolynomialInMonomialBasis<Value, Argument, higher_degree_,
                                              Evaluator_>(this Self&& self);

  // A polynomial may be explicitly converted to use a different evaluator.
  template<template<typename, typename, int> typename OtherEvaluator,
           typename Self>
  explicit operator PolynomialInMonomialBasis<Value, Argument, degree_,
                                              OtherEvaluator>(this Self&& self);

  PolynomialInMonomialBasis& operator+=(const PolynomialInMonomialBasis& right);
  PolynomialInMonomialBasis& operator-=(const PolynomialInMonomialBasis& right);

  template<typename S>
  constexpr PolynomialInMonomialBasis& operator*=(S const& right)
    requires module<Value, S>;

  template<typename S>
  constexpr PolynomialInMonomialBasis& operator/=(S const& right)
    requires vector_space<Value, S>;

  Value PRINCIPIA_VECTORCALL operator()(Argument argument) const;
  Derivative<Value, Argument> PRINCIPIA_VECTORCALL EvaluateDerivative(
      Argument argument) const;

  void PRINCIPIA_VECTORCALL EvaluateWithDerivative(
      Argument argument,
      Value& value,
      Derivative<Value, Argument>& derivative) const;

  constexpr int degree() const override;
  bool is_zero() const override;

  Coefficients const& coefficients() const;
  Argument const& origin() const;

  // Returns a copy of this polynomial adjusted to the given origin.
  PolynomialInMonomialBasis AtOrigin(Argument const& origin) const;

  template<int order = 1>
  PolynomialInMonomialBasis<
      Derivative<Value, Argument, order>, Argument, degree_ - order, Evaluator_>
  Derivative() const;

  // The constant term of the result is zero.
  PolynomialInMonomialBasis<Primitive<Value, Argument>, Argument, degree_ + 1,
                            Evaluator_>
  Primitive() const
    requires additive_group<Value>;

  base::_algebra::Primitive<Value, Argument> Integrate(
      Argument const& argument1,
      Argument const& argument2) const
    requires additive_group<Value>;

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;

  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 protected:
  // We redefined the good name `Derivative` above, so here we are...
  template<typename V, typename A>
  using D = base::_algebra::Derivative<V, A>;

  Value PRINCIPIA_VECTORCALL VirtualEvaluate(Argument argument) const override;
  D<Value, Argument> PRINCIPIA_VECTORCALL
  VirtualEvaluateDerivative(Argument argument) const override;
  void PRINCIPIA_VECTORCALL
  VirtualEvaluateWithDerivative(Argument argument,
                                Value& value,
                                D<Value, Argument>& derivative) const override;

 private:
  Coefficients coefficients_;
  Argument origin_;

  template<additive_group V, affine A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator-(PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<additive_group V, affine A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E>
  friend operator+(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<affine V, affine A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E> friend operator+(
      PolynomialInMonomialBasis<Difference<V>, A, l, E> const& left,
      PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<affine V, affine A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r), E> friend operator+(
      PolynomialInMonomialBasis<V, A, l, E> const& left,
      PolynomialInMonomialBasis<Difference<V>, A, r, E> const& right);
  template<affine V, affine A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, std::max(l, r), E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename S,
           typename V, affine A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<S, V>, A, d, E>
  friend operator*(S const& left,
                   PolynomialInMonomialBasis<V, A, d, E> const& right);
  template<typename S,
           typename V, affine A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<V, S>, A, d, E>
  friend operator*(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename S,
           typename V, affine A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Quotient<V, S>, A, d, E>
  friend operator/(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename L, typename R, affine A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<L, R>, A, l + r, E>
  friend operator*(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
  template<typename V, affine A, int l,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, l, E> friend operator+(
      PolynomialInMonomialBasis<Difference<V>, A, l, E> const& left,
      V const& right);
#else
  template<typename V,
           std::same_as<Difference<V>> VDifference,
           affine A,
           int l,
           template<typename, typename, int>
           typename E>
  constexpr PolynomialInMonomialBasis<V, A, l, E> friend operator+(
      PolynomialInMonomialBasis<VDifference, A, l, E> const& left,
      V const& right);
#endif
  template<typename V, affine A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator+(
      V const& left,
      PolynomialInMonomialBasis<Difference<V>, A, r, E> const& right);
  template<typename V, affine A, int l,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, l, E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   V const& right);
  template<typename V, affine A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, r, E>
  friend operator-(V const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<L, A, l * r, E>
  friend Compose(PolynomialInMonomialBasis<L, R, l, E> const& left,
                 PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename L, typename R, affine A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<
      typename Hilbert<L, R>::InnerProductType, A, l + r, E>
  friend PointwiseInnerProduct(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename V, affine A, int d,
           template<typename, typename, int> typename E>
  friend std::ostream& operator<<(
      std::ostream& out,
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial);
  template<typename V, typename A, int d,
           template<typename, typename, int> typename E,
           typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
      O express_in);
};

// Vector space of polynomials.

template<additive_group Value, affine Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator_>
operator+(PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

template<additive_group Value, affine Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator_>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

template<additive_group Value, affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator_>
operator+(PolynomialInMonomialBasis<Value, Argument, ldegree_,
                                    Evaluator_> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

template<affine Value, affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator_>
operator+(PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree_,
                                    Evaluator_> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

template<affine Value, affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator_>
operator+(PolynomialInMonomialBasis<Value, Argument, ldegree_,
                                    Evaluator_> const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree_,
                                    Evaluator_> const& right);

template<affine Value, affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    std::max(ldegree_, rdegree_), Evaluator_>
operator-(PolynomialInMonomialBasis<Value, Argument, ldegree_,
                                    Evaluator_> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

// TODO(phl): The operators below work fine when `Value` is a vector space, but
// they are bogus when it is merely an affine space.  For instance, multiplying
// by a scalar a polynomial whose `Value` is a `Position` makes no sense.  We
// need to fix this.

template<typename Scalar,
         typename Value, affine Argument, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_,
                                    Evaluator_>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_,
                                    Evaluator_> const& right);

template<typename Scalar,
         typename Value, affine Argument, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_,
                                    Evaluator_>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_,
                                    Evaluator_> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, affine Argument, int degree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree_,
                                    Evaluator_>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_,
                                    Evaluator_> const& left,
          Scalar const& right);

// Algebra of polynomials.

template<typename LValue, typename RValue,
         affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                    ldegree_ + rdegree_, Evaluator_>
operator*(PolynomialInMonomialBasis<LValue, Argument, ldegree_,
                                    Evaluator_> const& left,
          PolynomialInMonomialBasis<RValue, Argument, rdegree_,
                                    Evaluator_> const& right);

// Additive operators polynomial ± constant.

#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
template<typename Value, affine Argument, int ldegree_,
         template<typename, typename, int>
         typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator_>
operator+(PolynomialInMonomialBasis<Difference<Value>,
                                    Argument,
                                    ldegree_,
                                    Evaluator_> const& left,
          Value const& right);
#else
template<typename Value,
         std::same_as<Difference<Value>> ValueDifference,
         affine Argument,
         int ldegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator_>
operator+(PolynomialInMonomialBasis<ValueDifference, Argument,
                                    ldegree_, Evaluator_> const& left,
          Value const& right);
#endif

template<typename Value, affine Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator_>
operator+(Value const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    rdegree_, Evaluator_> const& right);

template<typename Value, affine Argument, int ldegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree_,
                                    Evaluator_>
operator-(PolynomialInMonomialBasis<Value, Argument, ldegree_,
                                    Evaluator_> const& left,
          Value const& right);

template<typename Value, affine Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree_,
                                    Evaluator_>
operator-(Value const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_,
                                    Evaluator_> const& right);

// Application monoid.

template<typename LValue, typename RValue,
         typename RArgument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<LValue, RArgument, ldegree_ * rdegree_,
                                    Evaluator_>
Compose(PolynomialInMonomialBasis<LValue, RValue, ldegree_,
                                    Evaluator_> const& left,
        PolynomialInMonomialBasis<RValue, RArgument, rdegree_,
                                    Evaluator_> const& right);

// Returns a scalar polynomial obtained by pointwise inner product of two
// vector-valued polynomials.
template<typename LValue, typename RValue,
         affine Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator_>
constexpr PolynomialInMonomialBasis<
    typename Hilbert<LValue, RValue>::InnerProductType, Argument,
    ldegree_ + rdegree_, Evaluator_>
PointwiseInnerProduct(PolynomialInMonomialBasis<LValue, Argument, ldegree_,
                                                Evaluator_> const& left,
                      PolynomialInMonomialBasis<RValue, Argument, rdegree_,
                                                Evaluator_> const& right);

// Output.

template<typename Value, affine Argument, int degree_,
         template<typename, typename, int> typename Evaluator_>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree_,
                              Evaluator_> const& polynomial);

}  // namespace internal

using internal::DefaultEvaluator;
using internal::Policy;
using internal::PolynomialInMonomialBasis;

}  // namespace _polynomial_in_monomial_basis
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_in_monomial_basis_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_IN_MONOMIAL_BASIS_HPP_
