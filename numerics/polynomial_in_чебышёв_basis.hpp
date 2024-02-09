// The files containing the tree of child classes of |Polynomial| must be
// included in the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#include "numerics/polynomial.hpp"
#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_HPP_
#ifndef PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_
#define PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_

#include "absl/container/btree_set.h"
#include "base/not_null.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/numerics.pb.h"

// TODO(phl): Mathematica support.
namespace principia {

namespace serialization {
using PolynomialInЧебышёвBasis = PolynomialInChebyshevBasis;
}  // namespace serialization

namespace numerics {
namespace _polynomial_in_чебышёв_basis {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_polynomial;
using namespace principia::quantities::_named_quantities;

template<typename Value_, typename Argument_, int degree_>
class PolynomialInЧебышёвBasis : public Polynomial<Value_, Argument_> {
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

  Value operator()(Argument const& argument) const override;
  Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  Argument const& lower_bound() const;
  Argument const& upper_bound() const;

  // Returns the Frobenius companion matrix suitable for the Чебышёв basis.
  FixedMatrix<double, degree_, degree_>
  FrobeniusCompanionMatrix() const;

  // Returns true if this polynomial may (but doesn't necessarily) have real
  // roots.  Returns false it is guaranteed not to have real roots.  This is
  // significantly faster than calling |RealRoots|.  If |error_estimate| is
  // given, false is only returned if the envelope of the series at a distance
  // of |error_estimate| has no real roots.  This is useful if the series is an
  // approximation of some function with an L∞ error less than |error_estimate|.
  bool MayHaveRealRoots(Value error_estimate = Value{}) const;

  // Returns the real roots of the polynomial, computed as the eigenvalues of
  // the Frobenius companion matrix.
  absl::btree_set<Argument> RealRoots(double ε) const;

  void WriteToMessage(not_null<serialization::Polynomial*> message) const;
  static PolynomialInЧебышёвBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;
  Argument lower_bound_;
  Argument upper_bound_;
  Difference<Argument> width_;
  // Precomputed to save operations at the expense of some accuracy loss.
  Inverse<Difference<Argument>> one_over_width_;

  template<typename V, typename A, int d>
  friend constexpr bool operator==(
      PolynomialInЧебышёвBasis<V, A, d> const& left,
      PolynomialInЧебышёвBasis<V, A, d> const& right);
};

template<typename Value, typename Argument, int degree>
constexpr bool operator==(
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& left,
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& right);
template<typename Value, typename Argument, int degree>
constexpr bool operator!=(
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& left,
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& right);

}  // namespace internal

using internal::PolynomialInЧебышёвBasis;

}  // namespace _polynomial_in_чебышёв_basis
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_in_чебышёв_basis_body.hpp"

#endif  // PRINCIPIA_NUMERICS_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_HPP_
