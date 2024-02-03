#pragma once

#include "numerics/polynomial_in_чебышёв_basis.hpp"

namespace principia {
namespace numerics {
namespace polynomial_in_чебышёв_basis {
namespace internal {

template<typename Value_, typename Argument_, int degree_>
constexpr PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
PolynomialInЧебышёвBasis(Coefficients coefficients,
                         Argument const& lower_bound,
                         Argument const& upper_bound) {}

template<typename Value_, typename Argument_, int degree_>
Value_ PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::operator()(
    Argument const& argument) const {}

template<typename Value_, typename Argument_, int degree_>
Derivative<Value_, Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::EvaluateDerivative(
    Argument const& argument) const {}

template<typename Value_, typename Argument_, int degree_>
constexpr
int PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::degree() const {
  return degree_;
}

template<typename Value_, typename Argument_, int degree_>
bool PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::is_zero() const {
  return coefficients_ == Coefficients{};
}

template<typename Value_, typename Argument_, int degree_>
Argument_ const&
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::lower_bound() const {
  return lower_bound_;
}

template<typename Value_, typename Argument_, int degree_>
Argument_ const&
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::upper_bound() const {
  return upper_bound_;
}

template<typename Value_, typename Argument_, int degree_>
FixedMatrix<double, degree_ + 1, degree_ + 1>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::FrobeniusCompanionMatrix()
    const {
}

template<typename Value_, typename Argument_, int degree_>
bool PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::MayHaveRealRoots(
    Value const error_estimate) const {
}

template<typename Value_, typename Argument_, int degree_>
absl::btree_set<Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
RealRoots(double const ε) const {
}

template<typename Value_, typename Argument_, int degree_>
void PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::WriteToMessage(
    not_null<serialization::PolynomialInЧебышёвBasis*> message) const {}

template<typename Value_, typename Argument_, int degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::ReadFromMessage(
    serialization::PolynomialInЧебышёвBasis const& message) {
}

}  // namespace internal
}  // namespace polynomial_in_чебышёв_basis
}  // namespace numerics
}  // namespace principia
