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

#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/point.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

// The presence of an operator+ below causes a bizarre compilation error in
// seemingly unrelated code in PolynomialTestInMonomialBasis.VectorSpace.
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
                                _MSC_FULL_VER == 193'632'535)

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value, typename Argument, int degree_) class,
    PolynomialInMonomialBasis,
    FROM(polynomial_in_monomial_basis));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value, typename Argument, int degree_,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_polynomial_in_monomial_basis::
        PolynomialInMonomialBasis<Value, Argument, degree_> const&
            polynomial,
    OptionalExpressIn express_in),
    FROM(mathematica));
}  // namespace mathematica

namespace numerics {
namespace _polynomial_in_monomial_basis {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_point;
using namespace principia::numerics::_polynomial;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_tuples;

// Used to decide which evaluator to use for a particular polynomial.
class Policy {
 public:
  template<typename Value, typename Argument, int degree>
  PolynomialInMonomialBasis<Value, Argument, degree>&& WithEvaluator(
      PolynomialInMonomialBasis<Value, Argument, degree>&& polynomial) const;

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


template<typename Value_, typename Argument_, int degree_>
class PolynomialInMonomialBasis : public Polynomial<Value_, Argument_> {
 public:
  using Argument = Argument_;
  using Value = Value_;

  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients = Derivatives<Value, Argument, degree_ + 1>;

  // The coefficients are relative to origin; in other words they are applied to
  // powers of (argument - origin).

  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      Argument const& origin);
  template<template<typename, typename, int> typename Evaluator>
  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      Argument const& origin,
                                      with_evaluator_t<Evaluator>);

  constexpr explicit PolynomialInMonomialBasis(Coefficients coefficients)
    requires additive_group<Argument>;
  template<template<typename, typename, int> typename Evaluator>
  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      with_evaluator_t<Evaluator>)
    requires additive_group<Argument>;

  friend constexpr bool operator==(PolynomialInMonomialBasis const& left,
                                   PolynomialInMonomialBasis const& right) =
      default;
  friend constexpr bool operator!=(PolynomialInMonomialBasis const& left,
                                   PolynomialInMonomialBasis const& right) =
      default;

  // A polynomial may be explicitly converted to a higher degree.
  template<int higher_degree_>
  explicit operator PolynomialInMonomialBasis<Value, Argument, higher_degree_>()
      const;

  PolynomialInMonomialBasis& operator+=(const PolynomialInMonomialBasis& right);
  PolynomialInMonomialBasis& operator-=(const PolynomialInMonomialBasis& right);

  Value PRINCIPIA_VECTORCALL operator()(Argument argument) const final;
  Derivative<Value, Argument> PRINCIPIA_VECTORCALL EvaluateDerivative(
      Argument argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  Coefficients const& coefficients() const;
  Argument const& origin() const;

  // Returns a copy of this polynomial adjusted to the given origin.
  PolynomialInMonomialBasis AtOrigin(Argument const& origin) const;

  template<int order = 1>
  PolynomialInMonomialBasis<
      Derivative<Value, Argument, order>, Argument, degree_ - order>
  Derivative() const;

  // The constant term of the result is zero.
  PolynomialInMonomialBasis<Primitive<Value, Argument>, Argument, degree_ + 1>
  Primitive() const
    requires additive_group<Value>;

  quantities::_named_quantities::Primitive<Value, Argument> Integrate(
      Argument const& argument1,
      Argument const& argument2) const
    requires additive_group<Value>;

  // Changes the evaluator of this object.  Useful on the result of an operator
  // or of `ReadFromMessage`, as these functions use the default evaluator.
  template<template<typename, typename, int> typename Evaluator>
  PolynomialInMonomialBasis&& WithEvaluator() &&;

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;

  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);
  // Compatibility deserialization, when the evaluator is not present in
  // `message`.
  template<template<typename, typename, int> typename Evaluator>
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  static constexpr not_null<
      Evaluator<Value_, Difference<Argument_>, degree_> const*>
  DefaultEvaluator();

  // The evaluator is only nonnull on the compatibility path.
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message,
      Evaluator<Value, Difference<Argument>, degree_> const* evaluator);

  Coefficients coefficients_;
  Argument origin_;
  // TODO(phl): The `Evaluator` class should be able to take a `Point`.
  not_null<Evaluator<Value_, Difference<Argument_>, degree_> const*> evaluator_;

  template<typename V, typename A, int r>
  constexpr PolynomialInMonomialBasis<V, A, r>
  friend operator-(PolynomialInMonomialBasis<V, A, r> const& right);
  template<typename V, typename A, int l, int r>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r)>
  friend operator+(PolynomialInMonomialBasis<V, A, l> const& left,
                   PolynomialInMonomialBasis<V, A, r> const& right);
  template<typename V, typename A, int l, int r>
  constexpr PolynomialInMonomialBasis<V, A, std::max(l, r)>
  friend operator-(PolynomialInMonomialBasis<V, A, l> const& left,
                   PolynomialInMonomialBasis<V, A, r> const& right);
  template<typename S,
           typename V, typename A, int d>
  constexpr PolynomialInMonomialBasis<Product<S, V>, A, d>
  friend operator*(S const& left,
                   PolynomialInMonomialBasis<V, A, d> const& right);
  template<typename S,
           typename V, typename A, int d>
  constexpr PolynomialInMonomialBasis<Product<V, S>, A, d>
  friend operator*(PolynomialInMonomialBasis<V, A, d> const& left,
                   S const& right);
  template<typename S,
           typename V, typename A, int d>
  constexpr PolynomialInMonomialBasis<Quotient<V, S>, A, d>
  friend operator/(PolynomialInMonomialBasis<V, A, d> const& left,
                   S const& right);
  template<typename L, typename R, typename A,
           int l, int r>
  constexpr PolynomialInMonomialBasis<Product<L, R>, A, l + r>
  friend operator*(
      PolynomialInMonomialBasis<L, A, l> const& left,
      PolynomialInMonomialBasis<R, A, r> const& right);
#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
  template<typename V, typename A, int l>
  constexpr PolynomialInMonomialBasis<V, A, l>
  friend operator+(
      PolynomialInMonomialBasis<Difference<V>, A, l> const& left,
      V const& right);
#endif
  template<typename V, typename A, int r>
  constexpr PolynomialInMonomialBasis<V, A, r>
  friend operator+(
      V const& left,
      PolynomialInMonomialBasis<Difference<V>, A, r> const& right);
  template<typename V, typename A, int l>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, l>
  friend operator-(PolynomialInMonomialBasis<V, A, l> const& left,
                   V const& right);
  template<typename V, typename A, int r>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, r>
  friend operator-(V const& left,
                   PolynomialInMonomialBasis<V, A, r> const& right);
  template<typename L, typename R, typename A,
           int l, int r>
  constexpr PolynomialInMonomialBasis<L, A, l * r>
  friend Compose(PolynomialInMonomialBasis<L, R, l> const& left,
                 PolynomialInMonomialBasis<R, A, r> const& right);
  template<typename L, typename R, typename A,
           int l, int r>
  constexpr PolynomialInMonomialBasis<
      typename Hilbert<L, R>::InnerProductType, A, l + r>
  friend PointwiseInnerProduct(
      PolynomialInMonomialBasis<L, A, l> const& left,
      PolynomialInMonomialBasis<R, A, r> const& right);
  template<typename V, typename A, int d>
  friend std::ostream& operator<<(
      std::ostream& out,
      PolynomialInMonomialBasis<V, A, d> const& polynomial);
  template<typename V, typename A, int d,
           typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PolynomialInMonomialBasis<V, A, d> const& polynomial,
      O express_in);
};

// Vector space of polynomials.

template<typename Value, typename Argument, int rdegree_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_>
operator+(PolynomialInMonomialBasis<Value, Argument, rdegree_> const& right);

template<typename Value, typename Argument, int rdegree_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree_> const& right);

template<typename Value, typename Argument, int ldegree_, int rdegree_>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_)>
operator+(PolynomialInMonomialBasis<Value, Argument, ldegree_> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_> const& right);

template<typename Value, typename Argument, int ldegree_, int rdegree_>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    std::max(ldegree_, rdegree_)>
operator-(PolynomialInMonomialBasis<Value, Argument, ldegree_> const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_> const& right);

template<typename Scalar,
         typename Value, typename Argument, int degree_>
constexpr PolynomialInMonomialBasis<Product<Scalar, Value>, Argument, degree_>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_> const& right);

template<typename Scalar,
         typename Value, typename Argument, int degree_>
constexpr PolynomialInMonomialBasis<Product<Value, Scalar>, Argument, degree_>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_> const& left,
          Scalar const& right);

template<typename Scalar,
         typename Value, typename Argument, int degree_>
constexpr PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument, degree_>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_> const& left,
          Scalar const& right);

// Algebra of polynomials.

template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_>
constexpr PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                    ldegree_ + rdegree_>
operator*(PolynomialInMonomialBasis<LValue, Argument, ldegree_> const& left,
          PolynomialInMonomialBasis<RValue, Argument, rdegree_> const& right);

// Additive operators polynomial ± constant.

#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
template<typename Value, typename Argument, int ldegree_>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree_>
operator+(PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    ldegree_> const& left,
          Value const& right);
#endif

template<typename Value, typename Argument, int rdegree_>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_>
operator+(Value const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    rdegree_> const& right);

template<typename Value, typename Argument, int ldegree_>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, ldegree_>
operator-(PolynomialInMonomialBasis<Value, Argument, ldegree_> const& left,
          Value const& right);

template<typename Value, typename Argument, int rdegree_>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument, rdegree_>
operator-(Value const& left,
          PolynomialInMonomialBasis<Value, Argument, rdegree_> const& right);

// Application monoid.

template<typename LValue, typename RValue,
         typename RArgument, int ldegree_, int rdegree_>
constexpr PolynomialInMonomialBasis<LValue, RArgument, ldegree_ * rdegree_>
Compose(PolynomialInMonomialBasis<LValue, RValue, ldegree_> const& left,
        PolynomialInMonomialBasis<RValue, RArgument, rdegree_> const& right);

// Returns a scalar polynomial obtained by pointwise inner product of two
// vector-valued polynomials.
template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_>
constexpr PolynomialInMonomialBasis<
    typename Hilbert<LValue, RValue>::InnerProductType, Argument,
    ldegree_ + rdegree_>
PointwiseInnerProduct(
    PolynomialInMonomialBasis<LValue, Argument, ldegree_> const& left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_> const& right);

// Output.

template<typename Value, typename Argument, int degree_>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree_> const& polynomial);

}  // namespace internal

using internal::Policy;
using internal::PolynomialInMonomialBasis;
using internal::with_evaluator;

}  // namespace _polynomial_in_monomial_basis
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_in_monomial_basis_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_IN_MONOMIAL_BASIS_HPP_
