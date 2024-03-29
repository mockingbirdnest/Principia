// The files containing the tree of child classes of |Polynomial| must be
// included in the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#include "numerics/polynomial.hpp"
#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_
#define PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_

#include <memory>
#include <optional>
#include <string>

#include "absl/container/btree_set.h"
#include "base/not_null.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value, typename Argument, auto degree) class,
    PolynomialInЧебышёвBasis,
    FROM(polynomial_in_чебышёв_basis));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value,
             typename Argument,
             int degree,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_polynomial_in_чебышёв_basis::
         PolynomialInЧебышёвBasis<Value, Argument, degree> const& series,
     OptionalExpressIn express_in),
    FROM(mathematica));
}  // namespace mathematica

namespace serialization {
using PolynomialInЧебышёвBasis = PolynomialInChebyshevBasis;
using ЧебышёвSeries = ChebyshevSeries;
}  // namespace serialization

namespace numerics {
namespace _polynomial_in_чебышёв_basis {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_polynomial;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;

template<typename Value_, typename Argument_, auto _ = std::nullopt>
class PolynomialInЧебышёвBasis;

// Degree-agnostic base class defining the contract of polynomials in the
// Чебышёв basis and used for polymorphic storage.
template<typename Value_, typename Argument_>
class PolynomialInЧебышёвBasis<Value_, Argument_, std::nullopt>
    : public Polynomial<Value_, Argument_> {
 public:
  using Argument = Argument_;
  using Value = Value_;

  friend constexpr bool operator==(PolynomialInЧебышёвBasis const& left,
                                   PolynomialInЧебышёвBasis const& right) =
      default;
  friend constexpr bool operator!=(PolynomialInЧебышёвBasis const& left,
                                   PolynomialInЧебышёвBasis const& right) =
      default;

  // Bounds passed at construction.
  virtual Argument const& lower_bound() const = 0;
  virtual Argument const& upper_bound() const = 0;

  // Returns true if this polynomial may (but doesn't necessarily) have real
  // roots.  Returns false it is guaranteed not to have real roots.  This is
  // significantly faster than calling |RealRoots|.  If |error_estimate| is
  // given, false is only returned if the envelope of the series at a distance
  // of |error_estimate| has no real roots.  This is useful if the series is an
  // approximation of some function with an L∞ error less than |error_estimate|.
  bool MayHaveRealRoots(Value error_estimate = Value{}) const
    requires convertible_to_quantity<Value_>;

  // Returns the real roots of the polynomial, computed as the eigenvalues of
  // the Frobenius companion matrix.
  absl::btree_set<Argument> RealRoots(double ε) const
    requires convertible_to_quantity<Value_>;

  // Compatibility deserialization: this class is the equivalent of the old
  // ЧебышёвSeries.
  static std::unique_ptr<PolynomialInЧебышёвBasis> ReadFromMessage(
      serialization::ЧебышёвSeries const& pre_канторович_message);

 protected:
  // Beware! The concept of real roots is only meaningful for scalar-valued
  // polynomials.  These functions will fail if called on a polynomial that is
  // not scalar-valued.
  virtual bool MayHaveRealRootsOrDie(Value error_estimate) const = 0;
  virtual absl::btree_set<Argument> RealRootsOrDie(double ε) const = 0;
};

template<typename Value_, typename Argument_, int degree_>
class PolynomialInЧебышёвBasis<Value_, Argument_, degree_>
    : public PolynomialInЧебышёвBasis<Value_, Argument_> {
 public:
  static_assert(degree_ >= 0);
  using Argument = Argument_;
  using Value = Value_;

  // The elements of the basis are dimensionless, unlike the monomial basis.
  using Coefficients = std::array<Value, degree_ + 1>;

  // The polynomials are only defined over [lower_bound, upper_bound].
  constexpr PolynomialInЧебышёвBasis(Coefficients coefficients,
                                     Argument const& lower_bound,
                                     Argument const& upper_bound);

  friend constexpr bool operator==(PolynomialInЧебышёвBasis const& left,
                                   PolynomialInЧебышёвBasis const& right) =
      default;
  friend constexpr bool operator!=(PolynomialInЧебышёвBasis const& left,
                                   PolynomialInЧебышёвBasis const& right) =
      default;

  Value __vectorcall operator()(Argument argument) const override;
  Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  Argument const& lower_bound() const override;
  Argument const& upper_bound() const override;

  // Returns the Frobenius companion matrix suitable for the Чебышёв basis.
  FixedMatrix<double, degree_, degree_> FrobeniusCompanionMatrix() const
    requires convertible_to_quantity<Value_>;

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInЧебышёвBasis ReadFromMessage(
      serialization::Polynomial const& message);

 protected:
  bool MayHaveRealRootsOrDie(Value error_estimate) const override;
  absl::btree_set<Argument> RealRootsOrDie(double ε) const override;

 private:
  Coefficients coefficients_;
  Argument lower_bound_;
  Argument upper_bound_;
  Difference<Argument> width_;
  // Precomputed to save operations at the expense of some accuracy loss.
  Inverse<Difference<Argument>> one_over_width_;

  template<typename V, typename A, int d, typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PolynomialInЧебышёвBasis<V, A, d> const& series,
      O express_in);
};

}  // namespace internal

using internal::PolynomialInЧебышёвBasis;

}  // namespace _polynomial_in_чебышёв_basis
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_in_чебышёв_basis_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_
