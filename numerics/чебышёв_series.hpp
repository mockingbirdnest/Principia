
#pragma once

#include <vector>

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

// Spelling: Чебышёв ЧЕБЫШЁВ чебышёв
namespace principia {

namespace serialization {
using ЧебышёвSeries = ChebyshevSeries;
}  // namespace serialization

namespace numerics {
namespace internal_чебышёв_series {

using base::not_null;
using geometry::Instant;
using quantities::Inverse;
using quantities::Time;
using quantities::Variation;

// A helper class for implementing |Evaluate| that can be specialized for speed.
template<typename Vector>
class EvaluationHelper final {
 public:
  EvaluationHelper(std::vector<Vector> const& coefficients, int degree);
  EvaluationHelper(EvaluationHelper&& other) = default;
  EvaluationHelper& operator=(EvaluationHelper&& other) = default;

  Vector EvaluateImplementation(double scaled_t) const;

  Vector coefficients(int index) const;
  int degree() const;

 private:
  std::vector<Vector> coefficients_;
  int degree_;
};

// A Чебышёв series with values in the vector space |Vector|.  The argument is
// an |Instant|.
template<typename Vector>
class ЧебышёвSeries final {
 public:
  // The element at position i in |coefficients| is the coefficient of Tᵢ.  The
  // polynomials are scaled to the interval [t_min, t_max], which must be
  // nonempty.
  ЧебышёвSeries(std::vector<Vector> const& coefficients,
                Instant const& t_min,
                Instant const& t_max);
  ЧебышёвSeries(ЧебышёвSeries&& other) = default;
  ЧебышёвSeries& operator=(ЧебышёвSeries&& other) = default;

  bool operator==(ЧебышёвSeries const& right) const;
  bool operator!=(ЧебышёвSeries const& right) const;

  Instant const& t_min() const;
  Instant const& t_max() const;

  // Only useful for benchmarking or analyzing performance.  Do not use in real
  // code.
  int degree() const;

  // Uses the Clenshaw algorithm.  |t| must be in the range [t_min, t_max].
  Vector Evaluate(Instant const& t) const;
  Variation<Vector> EvaluateDerivative(Instant const& t) const;

  void WriteToMessage(not_null<serialization::ЧебышёвSeries*> message) const;
  static ЧебышёвSeries ReadFromMessage(
      serialization::ЧебышёвSeries const& message);

 private:
  Instant t_min_;
  Instant t_max_;
  Inverse<Time> one_over_duration_;
  EvaluationHelper<Vector> helper_;
};

}  // namespace internal_чебышёв_series

using internal_чебышёв_series::ЧебышёвSeries;

}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_series_body.hpp"
