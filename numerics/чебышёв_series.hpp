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

using geometry::Instant;
using quantities::Time;
using quantities::Variation;

namespace numerics {

// A Чебышёв series with values in the vector space |Vector|.  The argument is
// an |Instant|.
template<typename Vector>
class ЧебышёвSeries {
 public:
  // The element at position i in |coefficients| is the coefficient of Tᵢ.  The
  // polynomials are scaled to the interval [t_min, t_max], which must be
  // nonempty.
  explicit ЧебышёвSeries(std::vector<Vector> const& coefficients,
                         Instant const& t_min,
                         Instant const& t_max);

  bool operator==(ЧебышёвSeries const& right) const;
  bool operator!=(ЧебышёвSeries const& right) const;

  // Uses the Clenshaw algorithm.  |t| must be in the range [t_min, t_max].
  Vector Evaluate(Instant const& t) const;
  Variation<Vector> EvaluateDerivative(Instant const& t) const;

  void WriteToMessage(
      not_null<serialization::ЧебышёвSeries*> const message) const;
  static ЧебышёвSeries ReadFromMessage(
      serialization::ЧебышёвSeries const& message);

  // Computes a Newhall approximation of the given |degree|.  |q| and |v| are
  // the positions and velocities over a constant division of [t_min, t_max].
  static ЧебышёвSeries NewhallApproximation(
      int const degree,
      std::vector<Vector> const& q,
      std::vector<Variation<Vector>> const& v,
      Instant const& t_min,
      Instant const& t_max);

 private:
  std::vector<Vector> const coefficients_;
  int const degree_;
  Instant const t_min_;
  Instant const t_max_;
  Instant t_mean_;
  Time::Inverse two_over_duration_;
};

}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_series_body.hpp"
