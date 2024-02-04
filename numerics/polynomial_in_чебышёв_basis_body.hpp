#pragma once

#include "numerics/polynomial_in_чебышёв_basis.hpp"

#include "base/tags.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/serialization.hpp"
#include "numerics/combinatorics.hpp"
#include "numerics/matrix_computations.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_in_чебышёв_basis {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_tags;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_matrix_computations;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// A helper class for implementing |Evaluate| that can be specialized for speed.
template<typename Value_, int degree_>
struct EvaluationHelper : not_constructible {
  using Value = Value_;
  using Coefficients = std::array<Value, degree_ + 1>;

  static Value Evaluate(Coefficients const& coefficients,
                        double scaled_argument);
};

// The compiler does a much better job on an |R3Element<double>| than on a
// |Vector<Quantity>| so we specialize this case.
template<typename Scalar_, typename Frame_, int rank_, int degree_>
struct EvaluationHelper<Multivector<Scalar_, Frame_, rank_>, degree_>
    : not_constructible {
  using Value = Multivector<Scalar_, Frame_, rank_>;
  using Coefficients = std::array<Value, degree_ + 1>;

  static Value Evaluate(Coefficients const& coefficients,
                        double scaled_argument);
};

template<typename Value_, int degree_>
auto EvaluationHelper<Value_, degree_>::Evaluate(
    Coefficients const& coefficients,
    double const scaled_argument) -> Value {
  double const two_scaled_argument = scaled_argument + scaled_argument;
  Value const c_0 = coefficients[0];
  switch (degree_) {
    case 0:
      return c_0;
    case 1:
      return c_0 + scaled_argument * coefficients[1];
    default:
      // b_degree   = c_degree.
      Value b_i = coefficients[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      Value b_j = coefficients[degree_ - 1] + two_scaled_argument * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        b_i = coefficients[k + 1] + two_scaled_argument * b_j - b_i;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        b_j = coefficients[k] + two_scaled_argument * b_i - b_j;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients[1] + two_scaled_argument * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_argument * b_i - b_j;
      } else {
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_argument * b_j - b_i;
      }
  }
}

template<typename Scalar_, typename Frame_, int rank_, int degree_>
auto EvaluationHelper<Multivector<Scalar_, Frame_, rank_>, degree_>::Evaluate(
    Coefficients const& coefficients,
    double const scaled_argument) -> Value {
  double const two_scaled_argument = scaled_argument + scaled_argument;
  R3Element<double> const c_0 = coefficients[0];
  switch (degree_) {
    case 0:
      return Multivector<double, Frame_, rank_>(c_0) * si::Unit<Scalar_>;
    case 1:
      return Multivector<double, Frame_, rank_>(
                 c_0 + scaled_argument * coefficients[1]) * si::Unit<Scalar_>;
    default:
      // b_degree   = c_degree.
      R3Element<double> b_i = coefficients[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      R3Element<double> b_j =
          coefficients[degree_ - 1] + two_scaled_argument * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        R3Element<double> const c_kplus1 = coefficients[k + 1];
        b_i.x = c_kplus1.x + two_scaled_argument * b_j.x - b_i.x;
        b_i.y = c_kplus1.y + two_scaled_argument * b_j.y - b_i.y;
        b_i.z = c_kplus1.z + two_scaled_argument * b_j.z - b_i.z;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        R3Element<double> const c_k = coefficients[k];
        b_j.x = c_k.x + two_scaled_argument * b_i.x - b_j.x;
        b_j.y = c_k.y + two_scaled_argument * b_i.y - b_j.y;
        b_j.z = c_k.z + two_scaled_argument * b_i.z - b_j.z;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients[1] + two_scaled_argument * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame_, rank_>(
                   c_0 + scaled_argument * b_i - b_j) * si::Unit<Scalar_>;
      } else {
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame_, rank_>(
                   c_0 + scaled_argument * b_j - b_i) * si::Unit<Scalar_>;
      }
    }
}

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
    Argument const& argument) const {
  // This formula ensures continuity at the edges by producing -1 or +1 within
  // 2 ulps for |lower_bound_| and |upper_bound_|.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) *
      one_over_width_;
  // We have to allow |scaled_argument| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  DCHECK_LE(scaled_argument, 1.1);
  DCHECK_GE(scaled_argument, -1.1);

  return EvaluationHelper<Value, degree_>::Evaluate(coefficients_,
                                                    scaled_argument);
}

template<typename Value_, typename Argument_, int degree_>
Derivative<Value_, Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::EvaluateDerivative(
    Argument const& argument) const {
  // See comments above.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) *
      one_over_width_;
  double const two_scaled_argument = scaled_argument + scaled_argument;
  DCHECK_LE(scaled_argument, 1.1);
  DCHECK_GE(scaled_argument, -1.1);

  Value b_kplus2_vector{};
  Value b_kplus1_vector{};
  Value* b_kplus2 = &b_kplus2_vector;
  Value* b_kplus1 = &b_kplus1_vector;
  Value* const& b_k = b_kplus2;  // An overlay.
  for (int k = degree_ - 1; k >= 1; --k) {
    *b_k = coefficients_[k + 1] * (k + 1) +
           two_scaled_argument * *b_kplus1 - *b_kplus2;
    Value* const last_b_k = b_k;
    b_kplus2 = b_kplus1;
    b_kplus1 = last_b_k;
  }
  return (coefficients_[1] + two_scaled_argument * *b_kplus1 - *b_kplus2) *
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
FixedMatrix<double, degree_ + 1, degree_ + 1>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::FrobeniusCompanionMatrix()
    const {
  int const N = degree();
  FixedMatrix<double, degree_ + 1, degree_ + 1> A(uninitialized);

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
bool PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::MayHaveRealRoots(
    Value const error_estimate) const {
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
}

template<typename Value_, typename Argument_, int degree_>
absl::btree_set<Argument_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::
RealRoots(double const ε) const {
  auto const companion_matrix = FrobeniusCompanionMatrix();
  auto const real_schur_decomposition =
      RealSchurDecomposition(companion_matrix, ε);
  absl::btree_set<double> const& scaled_real_roots =
      real_schur_decomposition.real_eigenvalues;

  // Rescale from [-1, 1] to [lower_bound_, upper_bound_].
  absl::btree_set<Argument> real_roots;
  auto const midpoint = Barycentre(std::pair{lower_bound_, upper_bound_},
                                   std::pair{1, 1});
  auto const half_width = 0.5 * (upper_bound_ - lower_bound_);
  for (auto const& scaled_real_root : scaled_real_roots) {
    // Чебышёв polynomials don't make sense outside of [-1, 1] but they may
    // stil have roots there.
    if (-1 <= scaled_real_root && scaled_real_root <= 1) {
      real_roots.insert(scaled_real_root * half_width + midpoint);
    }
  }
  return real_roots;
}

template<typename Value_, typename Argument_, int degree_>
void PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::WriteToMessage(
    not_null<serialization::Polynomial*> message) const {
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
                                     extension->mutable_upper_bound());
}

template<typename Value_, typename Argument_, int degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>
PolynomialInЧебышёвBasis<Value_, Argument_, degree_>::ReadFromMessage(
    serialization::Polynomial const& message) {
  // TODO(phl): Add compatibility code with |ЧебышёвSeries|.
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

template<typename Value, typename Argument, int degree>
constexpr bool operator==(
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& left,
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& right) {
  return left.coefficients_ == right.coefficients_ &&
         left.lower_bound_ == right.lower_bound_ &&
         left.upper_bound_ == right.upper_bound_;
}

template<typename Value, typename Argument, int degree>
constexpr bool operator!=(
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& left,
    PolynomialInЧебышёвBasis<Value, Argument, degree> const& right) {
  return !operator==(left, right);
}

}  // namespace internal
}  // namespace _polynomial_in_чебышёв_basis
}  // namespace numerics
}  // namespace principia
