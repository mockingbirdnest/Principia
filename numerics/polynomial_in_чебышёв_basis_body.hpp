#pragma once

#include "numerics/polynomial_in_чебышёв_basis.hpp"

#include <memory>
#include <utility>

#include "base/not_constructible.hpp"
#include "base/tags.hpp"
#include "base/traits.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "numerics/combinatorics.hpp"
#include "numerics/matrix_computations.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_in_чебышёв_basis {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_tags;
using namespace principia::base::_traits;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_matrix_computations;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

template<typename Value_, typename Argument_>
bool PolynomialInЧебышёвBasis<Value_, Argument_, std::nullopt>::
    MayHaveRealRoots(Value const error_estimate) const
  requires convertible_to_quantity<Value_> {
  return MayHaveRealRootsOrDie(error_estimate);
}

template<typename Value_, typename Argument_>
absl::btree_set<Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, std::nullopt>::RealRoots(
    double const ε) const
  requires convertible_to_quantity<Value_> {
  return RealRootsOrDie(ε);
}

#define PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(degree) \
  case (degree):                                                             \
    return std::make_unique<                                                 \
        PolynomialInЧебышёвBasis<Value, Argument, degree>>(                  \
        PolynomialInЧебышёвBasis<Value, Argument, degree>::ReadFromMessage(  \
            message))

template<typename Value_, typename Argument_>
std::unique_ptr<PolynomialInЧебышёвBasis<Value_, Argument_, std::nullopt>>
PolynomialInЧебышёвBasis<Value_, Argument_, std::nullopt>::ReadFromMessage(
    serialization::ЧебышёвSeries const& pre_канторович_message) {
  LOG(WARNING) << "Reading pre-Канторович PolynomialInЧебышёвBasis";
  serialization::Polynomial message;
  message.set_degree(pre_канторович_message.coefficient_size() - 1);
  auto* const extension = message.MutableExtension(
      serialization::PolynomialInЧебышёвBasis::extension);
  for (auto const& coefficient : pre_канторович_message.coefficient()) {
    switch (coefficient.coefficient_case()) {
      case serialization::ЧебышёвSeries::Coefficient::kDouble:
        extension->add_coefficient()->set_double_(coefficient.double_());
        break;
      case serialization::ЧебышёвSeries::Coefficient::kQuantity:
        *extension->add_coefficient()->mutable_quantity() =
            coefficient.quantity();
        break;
      case serialization::ЧебышёвSeries::Coefficient::kMultivector:
        *extension->add_coefficient()->mutable_multivector() =
            coefficient.multivector();
        break;
      case serialization::ЧебышёвSeries::Coefficient::COEFFICIENT_NOT_SET:
        LOG(FATAL) << pre_канторович_message.DebugString();
    };
  }
  *extension->mutable_lower_bound()->mutable_point() =
      pre_канторович_message.lower_bound();
  *extension->mutable_upper_bound()->mutable_point() =
      pre_канторович_message.upper_bound();
  switch (pre_канторович_message.coefficient_size() - 1) {
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(0);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(1);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(2);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(3);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(4);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(5);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(6);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(7);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(8);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(9);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(10);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(11);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(12);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(13);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(14);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(15);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(16);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(17);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(18);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(19);
    PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE(20);
    default:
      LOG(FATAL) << "Unexpected degree: "
                 << pre_канторович_message.DebugString();
#if PRINCIPIA_COMPILER_MSVC && \
    (_MSC_FULL_VER == 193'933'523 || \
     _MSC_FULL_VER == 194'033'813 || \
     _MSC_FULL_VER == 194'134'120 || \
     _MSC_FULL_VER == 194'134'123)
      std::abort();
#endif
  }
}

#undef PRINCIPIA_POLYNOMIAL_IN_ЧЕБЫШЁВ_BASIS_DESERIALIZATION_DEGREE

template<typename Value_, typename Argument_, int degree_>
constexpr PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
PolynomialInЧебышёвBasis(Coefficients coefficients,
                         Argument const& lower_bound,
                         Argument const& upper_bound)
    : coefficients_(std::move(coefficients)),
      lower_bound_(lower_bound),
      upper_bound_(upper_bound),
      width_(upper_bound - lower_bound),
      one_over_width_(1 / width_) {}

template<typename Value_, typename Argument_, int degree_>
Value_ PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::operator()(
    Argument const argument) const {
  // This formula ensures continuity at the edges by producing -1 or +1 within
  // 2 ulps for `lower_bound_` and `upper_bound_`.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) * one_over_width_;
  // We have to allow `scaled_argument` to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  DCHECK_LE(scaled_argument, 1.1);
  DCHECK_GE(scaled_argument, -1.1);

  double const two_scaled_argument = scaled_argument + scaled_argument;
  Value const& c₀ = coefficients_[0];
  switch (degree_) {
    case 0:
      return c₀;
    case 1:
      return c₀ + scaled_argument * coefficients_[1];
    default:
      Value bₖ₊₂ = coefficients_[degree_];
      Value bₖ₊₁ = coefficients_[degree_ - 1] + two_scaled_argument * bₖ₊₂;
      for (int k = degree_ - 2; k >= 1; --k) {
        Value const bₖ = coefficients_[k] + two_scaled_argument * bₖ₊₁ - bₖ₊₂;
        bₖ₊₂ = bₖ₊₁;
        bₖ₊₁ = bₖ;
      }
      return c₀ + scaled_argument * bₖ₊₁ - bₖ₊₂;
  }
}

template<typename Value_, typename Argument_, int degree_>
Derivative<Value_, Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::EvaluateDerivative(
    Argument const argument) const {
  // See comments above.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) * one_over_width_;
  double const two_scaled_argument = scaled_argument + scaled_argument;
  DCHECK_LE(scaled_argument, 1.1);
  DCHECK_GE(scaled_argument, -1.1);

  Value bₖ₊₂{};
  Value bₖ₊₁{};
  for (int k = degree_ - 1; k >= 1; --k) {
    Value const bₖ =
        coefficients_[k + 1] * (k + 1) + two_scaled_argument * bₖ₊₁ - bₖ₊₂;
    bₖ₊₂ = bₖ₊₁;
    bₖ₊₁ = bₖ;
  }
  return (coefficients_[1] + two_scaled_argument * bₖ₊₁ - bₖ₊₂) *
         (one_over_width_ + one_over_width_);
}

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
FixedMatrix<double, degree_, degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::FrobeniusCompanionMatrix()
    const requires convertible_to_quantity<Value_> {
  int const N = degree();
  FixedMatrix<double, degree_, degree_> A(uninitialized);

  // j = 1.
  for (int k = 1; k <= N; ++k) {
    auto& Aⱼₖ = A(0, k - 1);
    Aⱼₖ = KroneckerDelta(2, k);
  }

  for (int j = 2; j < N; ++j) {
    for (int k = 1; k <= N; ++k) {
      auto& Aⱼₖ = A(j - 1, k - 1);
      Aⱼₖ = 0.5 * (KroneckerDelta(j, k + 1) + KroneckerDelta(j, k - 1));
    }
  }

  // j = N.
  auto const two_aN = 2 * coefficients_[N];
  for (int k = 1; k <= N; ++k) {
    auto& Aⱼₖ = A(N - 1, k - 1);
    // Note that [Boy13] formula (B.2) has aⱼ₋₁ instead of aₖ₋₁ below, but
    // that's probably a typo because it's immediately corrected in formula
    // (B.3).
    Aⱼₖ = -coefficients_[k - 1] / two_aN + 0.5 * KroneckerDelta(k, N - 1);
  }

  return A;
}

template<typename Value_, typename Argument_, int degree_>
void PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::WriteToMessage(
    not_null<serialization::Polynomial*> message) const {
  if constexpr (is_instance_of_v<R3Element, Value>) {
    // `R3Element` is a low-level structure, so we don't want polynomials to use
    // it directly: they should go through `Multivector`.  However, it's useful
    // to be able to run benchmarks using them.
    LOG(FATAL) << "R3Element only supported for tests";
  } else {
    message->set_degree(degree_);
    auto* const extension = message->MutableExtension(
        serialization::PolynomialInЧебышёвBasis::extension);

    using CoefficientSerializer = DoubleOrQuantityOrMultivectorSerializer<
        Value,
        serialization::PolynomialInЧебышёвBasis::Coefficient>;
    using ArgumentSerializer = DoubleOrQuantityOrPointOrMultivectorSerializer<
        Argument,
        serialization::PolynomialInЧебышёвBasis::Argument>;

    for (auto const& coefficient : coefficients_) {
      CoefficientSerializer::WriteToMessage(coefficient,
                                            extension->add_coefficient());
    }
    ArgumentSerializer::WriteToMessage(lower_bound_,
                                       extension->mutable_lower_bound());
    ArgumentSerializer::WriteToMessage(upper_bound_,
                                       extension->mutable_upper_bound());}
  }

template<typename Value_, typename Argument_, int degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::ReadFromMessage(
    serialization::Polynomial const& message) {
  CHECK_EQ(degree_, message.degree()) << message.DebugString();
  CHECK(
      message.HasExtension(serialization::PolynomialInЧебышёвBasis::extension))
      << message.DebugString();
  auto const& extension =
      message.GetExtension(serialization::PolynomialInЧебышёвBasis::extension);

  using CoefficientSerializer = DoubleOrQuantityOrMultivectorSerializer<
      Value,
      serialization::PolynomialInЧебышёвBasis::Coefficient>;
  using ArgumentSerializer = DoubleOrQuantityOrPointOrMultivectorSerializer<
      Argument,
      serialization::PolynomialInЧебышёвBasis::Argument>;

  Coefficients coefficients;
  for (int i = 0; i < coefficients.size(); ++i) {
    coefficients[i] =
        CoefficientSerializer::ReadFromMessage(extension.coefficient(i));
  }
  return PolynomialInЧебышёвBasis(
      coefficients,
      ArgumentSerializer::ReadFromMessage(extension.lower_bound()),
      ArgumentSerializer::ReadFromMessage(extension.upper_bound()));
}

template<typename Value_, typename Argument_, int degree_>
bool PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
MayHaveRealRootsOrDie(Value const error_estimate) const {
  if constexpr (convertible_to_quantity<Value>) {
    CHECK_LE(Value{}, error_estimate);
    // This code follow [Boy06], theorem 2.  Note that [Boy06] has another
    // criterion, B₁ and concludes: “There was no detectable difference between
    // the two criteria, so the first criterion, which requires only summing the
    // absolute values of all coefficients but the first, should be used to the
    // exclusion of the other zero-free test”.  My own experiments failed to
    // locate, for N = 5, any series for which B₁ would be better than B₀, after
    // trying millions of random series.
    int const N = degree();
    Value B₀{};
    for (int j = 1; j <= N; ++j) {
      auto const abs_aⱼ = Abs(coefficients_[j]);
      B₀ += abs_aⱼ;
    }
    auto const abs_a₀ = Abs(coefficients_[0]);
    // The error may shift the curve vertically.  Note that the following
    // comparison is valid if the right-hand side is negative.
    return B₀ >= abs_a₀ - error_estimate;
  } else {
    LOG(FATAL) << "Real roots only meaningful for scalar-valued polynomials";
  }
}

template<typename Value_, typename Argument_, int degree_>
absl::btree_set<Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
RealRootsOrDie(double const ε) const {
  if constexpr (convertible_to_quantity<Value>) {
    auto const companion_matrix = FrobeniusCompanionMatrix();
    auto const real_schur_decomposition =
        RealSchurDecomposition(companion_matrix, ε);
    absl::btree_set<double> const& scaled_real_roots =
        real_schur_decomposition.real_eigenvalues;

    // Rescale from [-1, 1] to [lower_bound_, upper_bound_].
    absl::btree_set<Argument> real_roots;
    auto const midpoint = Barycentre({lower_bound_, upper_bound_});
    auto const half_width = 0.5 * (upper_bound_ - lower_bound_);
    for (auto const& scaled_real_root : scaled_real_roots) {
      // Чебышёв polynomials don't make sense outside of [-1, 1] but they may
      // stil have roots there.
      if (-1 <= scaled_real_root && scaled_real_root <= 1) {
        real_roots.insert(scaled_real_root * half_width + midpoint);
      }
    }
    return real_roots;
  } else {
    LOG(FATAL) << "Real roots only meaningful for scalar-valued polynomials";
  }
}

}  // namespace internal
}  // namespace _polynomial_in_чебышёв_basis
}  // namespace numerics
}  // namespace principia
