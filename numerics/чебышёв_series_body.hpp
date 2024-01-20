#include "numerics/чебышёв_series.hpp"

#include <vector>

#include "base/tags.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "glog/logging.h"
#include "numerics/combinatorics.hpp"
#include "numerics/matrix_computations.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _чебышёв_series {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_serialization;
using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_matrix_computations;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// The compiler does a much better job on an |R3Element<double>| than on a
// |Vector<Quantity>| so we specialize this case.
template<typename Scalar, typename Frame, int rank>
class EvaluationHelper<Multivector<Scalar, Frame, rank>> final {
 public:
  EvaluationHelper(
      std::vector<Multivector<Scalar, Frame, rank>> const& coefficients,
      int degree);
  EvaluationHelper(EvaluationHelper&& other) = default;
  EvaluationHelper& operator=(EvaluationHelper&& other) = default;

  Multivector<Scalar, Frame, rank> EvaluateImplementation(
      double scaled_argument) const;

  Multivector<Scalar, Frame, rank> coefficients(int index) const;
  int degree() const;

 private:
  std::vector<R3Element<double>> coefficients_;
  int degree_;
};

template<typename Value>
EvaluationHelper<Value>::EvaluationHelper(
    std::vector<Value> const& coefficients,
    int const degree) : coefficients_(coefficients), degree_(degree) {}

template<typename Value>
Value EvaluationHelper<Value>::EvaluateImplementation(
    double const scaled_argument) const {
  double const two_scaled_argument = scaled_argument + scaled_argument;
  Value const c_0 = coefficients_[0];
  switch (degree_) {
    case 0:
      return c_0;
    case 1:
      return c_0 + scaled_argument * coefficients_[1];
    default:
      // b_degree   = c_degree.
      Value b_i = coefficients_[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      Value b_j = coefficients_[degree_ - 1] + two_scaled_argument * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        b_i = coefficients_[k + 1] + two_scaled_argument * b_j - b_i;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        b_j = coefficients_[k] + two_scaled_argument * b_i - b_j;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients_[1] + two_scaled_argument * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_argument * b_i - b_j;
      } else {
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_argument * b_j - b_i;
      }
  }
}

template<typename Value>
Value EvaluationHelper<Value>::coefficients(int const index) const {
  return coefficients_[index];
}

template<typename Value>
int EvaluationHelper<Value>::degree() const {
  return degree_;
}

template<typename Scalar, typename Frame, int rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::EvaluationHelper(
    std::vector<Multivector<Scalar, Frame, rank>> const& coefficients,
    int const degree) : degree_(degree) {
  for (auto const& coefficient : coefficients) {
    coefficients_.push_back(coefficient.coordinates() / si::Unit<Scalar>);
  }
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::EvaluateImplementation(
    double const scaled_argument) const {
  double const two_scaled_argument = scaled_argument + scaled_argument;
  R3Element<double> const c_0 = coefficients_[0];
  switch (degree_) {
    case 0:
      return Multivector<double, Frame, rank>(c_0) * si::Unit<Scalar>;
    case 1:
      return Multivector<double, Frame, rank>(
                 c_0 + scaled_argument * coefficients_[1]) * si::Unit<Scalar>;
    default:
      // b_degree   = c_degree.
      R3Element<double> b_i = coefficients_[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      R3Element<double> b_j =
          coefficients_[degree_ - 1] + two_scaled_argument * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        R3Element<double> const c_kplus1 = coefficients_[k + 1];
        b_i.x = c_kplus1.x + two_scaled_argument * b_j.x - b_i.x;
        b_i.y = c_kplus1.y + two_scaled_argument * b_j.y - b_i.y;
        b_i.z = c_kplus1.z + two_scaled_argument * b_j.z - b_i.z;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        R3Element<double> const c_k = coefficients_[k];
        b_j.x = c_k.x + two_scaled_argument * b_i.x - b_j.x;
        b_j.y = c_k.y + two_scaled_argument * b_i.y - b_j.y;
        b_j.z = c_k.z + two_scaled_argument * b_i.z - b_j.z;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients_[1] + two_scaled_argument * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame, rank>(
                   c_0 + scaled_argument * b_i - b_j) * si::Unit<Scalar>;
      } else {
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame, rank>(
                   c_0 + scaled_argument * b_j - b_i) * si::Unit<Scalar>;
      }
    }
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::coefficients(
    int const index) const {
  return Multivector<double, Frame, rank>(
             coefficients_[index]) * si::Unit<Scalar>;
}

template<typename Scalar, typename Frame, int rank>
int EvaluationHelper<Multivector<Scalar, Frame, rank>>::degree() const {
  return degree_;
}

template<typename Value, typename Argument>
ЧебышёвSeries<Value, Argument>::ЧебышёвSeries(
    std::vector<Value> const& coefficients,
    Argument const& lower_bound,
    Argument const& upper_bound)
    : lower_bound_(lower_bound),
      upper_bound_(upper_bound),
      helper_(coefficients,
              /*degree=*/static_cast<int>(coefficients.size()) - 1) {
  CHECK_LE(0, helper_.degree()) << "Degree must be at least 0";
  CHECK_LT(lower_bound_, upper_bound_) << "Argument interval must not be empty";
  // Precomputed to save operations at the expense of some accuracy loss.
  auto const width = upper_bound_ - lower_bound_;
  one_over_width_ = 1 / width;
}


template<typename Value, typename Argument>
bool ЧебышёвSeries<Value, Argument>::operator==(
    ЧебышёвSeries const& right) const {
  if (helper_.degree() != right.helper_.degree()) {
    return false;
  }
  for (int k = 0; k < helper_.degree(); ++k) {
    if (helper_.coefficients(k) != right.helper_.coefficients(k)) {
      return false;
    }
  }
  return lower_bound_ == right.lower_bound_ &&
         upper_bound_ == right.upper_bound_;
}

template<typename Value, typename Argument>
bool ЧебышёвSeries<Value, Argument>::operator!=(
    ЧебышёвSeries const& right) const {
  return !ЧебышёвSeries::operator==(right);
}

template<typename Value, typename Argument>
Argument const& ЧебышёвSeries<Value, Argument>::lower_bound() const {
  return lower_bound_;
}

template<typename Value, typename Argument>
Argument const& ЧебышёвSeries<Value, Argument>::upper_bound() const {
  return upper_bound_;
}

template<typename Value, typename Argument>
int ЧебышёвSeries<Value, Argument>::degree() const {
  return helper_.degree();
}

template<typename Value, typename Argument>
Value ЧебышёвSeries<Value, Argument>::Evaluate(Argument const& argument) const {
  // This formula ensures continuity at the edges by producing -1 or +1 within
  // 2 ulps for |lower_bound_| and |upper_bound_|.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) *
      one_over_width_;
  // We have to allow |scaled_argument| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  DCHECK_LE(scaled_argument, 1.1);
  DCHECK_GE(scaled_argument, -1.1);

  return helper_.EvaluateImplementation(scaled_argument);
}

template<typename Value, typename Argument>
Derivative<Value, Argument> ЧебышёвSeries<Value, Argument>::EvaluateDerivative(
    Argument const& argument) const {
  // See comments above.
  double const scaled_argument =
      ((argument - upper_bound_) + (argument - lower_bound_)) *
      one_over_width_;
  double const two_scaled_argument = scaled_argument + scaled_argument;
#ifdef _DEBUG
  CHECK_LE(scaled_argument, 1.1);
  CHECK_GE(scaled_argument, -1.1);
#endif

  Value b_kplus2_vector{};
  Value b_kplus1_vector{};
  Value* b_kplus2 = &b_kplus2_vector;
  Value* b_kplus1 = &b_kplus1_vector;
  Value* const& b_k = b_kplus2;  // An overlay.
  for (int k = helper_.degree() - 1; k >= 1; --k) {
    *b_k = helper_.coefficients(k + 1) * (k + 1) +
           two_scaled_argument * *b_kplus1 - *b_kplus2;
    Value* const last_b_k = b_k;
    b_kplus2 = b_kplus1;
    b_kplus1 = last_b_k;
  }
  return (helper_.coefficients(1) + two_scaled_argument * *b_kplus1 -
          *b_kplus2) *
         (one_over_width_ + one_over_width_);
}

template<typename Value, typename Argument>
UnboundedMatrix<double>
ЧебышёвSeries<Value, Argument>::FrobeniusCompanionMatrix() const {
  int const N = degree();
  UnboundedMatrix<double> A(N, N, uninitialized);

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
  auto const two_aN = 2 * helper_.coefficients(N);
  for (int k = 1; k <= N; ++k) {
    auto& Aⱼₖ = A(N - 1, k - 1);
    // Note that [Boy13] formula (B.2) has aⱼ₋₁ instead of aₖ₋₁ below, but
    // that's probably a typo because it's immediately corrected in formula
    // (B.3).
    Aⱼₖ =
        -helper_.coefficients(k - 1) / two_aN + 0.5 * KroneckerDelta(k, N - 1);
  }

  return A;
}

template<typename Value, typename Argument>
bool ЧебышёвSeries<Value, Argument>::MayHaveRealRoots() const {
  // TODO(phl): This is a property of the series so we should cache it.
  // This code follow [Boy06], theorem 2.
  int const N = degree();
  Value B₀{};
  Value B₁{};
  for (int j = 1; j <= N; ++j) {
    auto const Abs_aⱼ = Abs(helper_.coefficients(j));
    B₀ += Abs_aⱼ;
    B₁ += j * Abs_aⱼ;
  }
  if (B₀ < Abs(helper_.coefficients(0))) {
    return false;
  }
  double const h = π / N;
  auto const B₁h_over_2 = 0.5 * B₁ * h;
  for (int j = 0; j <= N; ++j) {
    // Note that [Boy06] has a typo in the definition of tⱼ.
    auto const tⱼ = j * h;
    auto const fⱼ = Evaluate(Cos(tⱼ));
    if (fⱼ <= B₁h_over_2) {
      return true;
    }
  }
  return false;
}

template<typename Value, typename Argument>
absl::btree_set<Argument> ЧебышёвSeries<Value, Argument>::RealRoots(
    double const ε) const {
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

template<typename Value, typename Argument>
void ЧебышёвSeries<Value, Argument>::WriteToMessage(
    not_null<serialization::ЧебышёвSeries*> const message) const {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Value,
                          serialization::ЧебышёвSeries::Coefficient>;

  for (int k = 0; k <= helper_.degree(); ++k) {
    Serializer::WriteToMessage(helper_.coefficients(k),
                               message->add_coefficient());
  }
  lower_bound_.WriteToMessage(message->mutable_lower_bound());
  upper_bound_.WriteToMessage(message->mutable_upper_bound());
}

template<typename Value, typename Argument>
ЧебышёвSeries<Value, Argument> ЧебышёвSeries<Value, Argument>::ReadFromMessage(
    serialization::ЧебышёвSeries const& message) {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Value,
                          serialization::ЧебышёвSeries::Coefficient>;

  std::vector<Value> coefficients;
  coefficients.reserve(message.coefficient_size());
  for (auto const& coefficient : message.coefficient()) {
    coefficients.push_back(Serializer::ReadFromMessage(coefficient));
  }
  return ЧебышёвSeries(coefficients,
                       Argument::ReadFromMessage(message.lower_bound()),
                       Argument::ReadFromMessage(message.upper_bound()));
}

}  // namespace internal
}  // namespace _чебышёв_series
}  // namespace numerics
}  // namespace principia
