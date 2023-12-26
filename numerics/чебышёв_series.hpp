#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

// Spelling: Чебышёв ЧЕБЫШЁВ чебышёв
namespace principia {

namespace serialization {
using ЧебышёвSeries = ChebyshevSeries;
}  // namespace serialization

namespace numerics {
namespace _чебышёв_series {
namespace internal {

using namespace principia::base::_not_null;
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
  template<int degree>
  ЧебышёвSeries(FixedVector<Value, degree + 1> const& coefficients,
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

  void WriteToMessage(not_null<serialization::ЧебышёвSeries*> message) const;
  static ЧебышёвSeries ReadFromMessage(
      serialization::ЧебышёвSeries const& message);

 private:
  Argument lower_bound_;
  Argument upper_bound_;
  Inverse<Difference<Argument>> one_over_width_;
  EvaluationHelper<Value> helper_;
};

}  // namespace internal

using internal::ЧебышёвSeries;

}  // namespace _чебышёв_series
}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_series_body.hpp"
