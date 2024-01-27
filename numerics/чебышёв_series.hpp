#pragma once

#include <string>
#include <vector>

#include "absl/container/btree_set.h"
#include "base/macros.hpp"  // 🧙 For forward declarations.
#include "base/not_null.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

// Spelling: Чебышёв ЧЕБЫШЁВ чебышёв

namespace principia {
namespace numerics {
FORWARD_DECLARE(
    TEMPLATE(typename Value, typename Argument) class,
    ЧебышёвSeries,
    FROM(чебышёв_series));
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION(
    TEMPLATE(typename Value,
             typename Argument,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_чебышёв_series::ЧебышёвSeries<Value, Argument> const& series,
     OptionalExpressIn express_in),
    FROM(mathematica));
}  // namespace mathematica

namespace serialization {
using ЧебышёвSeries = ChebyshevSeries;
}  // namespace serialization

namespace numerics {
namespace _чебышёв_series {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A helper class for implementing |Evaluate| that can be specialized for speed.
template<typename Value>
class EvaluationHelper final {
 public:
  EvaluationHelper(std::vector<Value> const& coefficients, int degree);
  EvaluationHelper(EvaluationHelper&& other) = default;
  EvaluationHelper& operator=(EvaluationHelper&& other) = default;

  Value EvaluateImplementation(double scaled_argument) const;

  Value coefficients(int index) const;
  int degree() const;

 private:
  std::vector<Value> coefficients_;
  int degree_;

  template<typename V, typename A, typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      ЧебышёвSeries<V, A> const& series,
      O express_in);
};

// A Чебышёв series with values in the vector space |Value|.  The argument is
// in the affine space |Argument|.
template<typename Value, typename Argument>
class ЧебышёвSeries final {
 public:
  // The element at position i in |coefficients| is the coefficient of Tᵢ.  The
  // polynomials are scaled to the interval [lower_bound, upper_bound], which
  // must be nonempty.
  ЧебышёвSeries(std::vector<Value> const& coefficients,
                Argument const& lower_bound,
                Argument const& upper_bound);
  ЧебышёвSeries(ЧебышёвSeries&& other) = default;
  ЧебышёвSeries& operator=(ЧебышёвSeries&& other) = default;

  bool operator==(ЧебышёвSeries const& right) const;
  bool operator!=(ЧебышёвSeries const& right) const;

  Argument const& lower_bound() const;
  Argument const& upper_bound() const;

  // Only useful for benchmarking or analyzing performance.  Do not use in real
  // code.
  int degree() const;

  // Uses the Clenshaw algorithm.  |argument| must be in the range
  // [lower_bound, upper_bound].
  Value Evaluate(Argument const& argument) const;
  Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const;

  // Returns the Frobenius companion matrix suitable for the Чебышёв basis.
  UnboundedMatrix<double> FrobeniusCompanionMatrix() const;

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

  void WriteToMessage(not_null<serialization::ЧебышёвSeries*> message) const;
  static ЧебышёвSeries ReadFromMessage(
      serialization::ЧебышёвSeries const& message);

 private:
  Argument lower_bound_;
  Argument upper_bound_;
  Inverse<Difference<Argument>> one_over_width_;
  EvaluationHelper<Value> helper_;

  template<typename V, typename A, typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      ЧебышёвSeries<V, A> const& series,
      O express_in);
};

}  // namespace internal

using internal::ЧебышёвSeries;

}  // namespace _чебышёв_series
}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_series_body.hpp"
