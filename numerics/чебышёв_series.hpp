#pragma once

#include <vector>

#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"

// Spelling: Чебышёв ЧЕБЫШЁВ чебышёв
namespace principia {

using geometry::Instant;
using quantities::Time;

namespace numerics {

template<typename Scalar>
class ЧебышёвSeries {
 public:
  // The element at position i in |coefficients| is the coefficient of Tᵢ.  The
  // polynomials are scaled to the interval [t_min, t_max], which must be
  // nonempty.
  explicit ЧебышёвSeries(std::vector<Scalar> const& coefficients,
                         Instant const& t_min,
                         Instant const& t_max);

  // Uses the Clenshaw algorithm.  |t| must be in the range [t_min, t_max].
  Scalar Evaluate(Instant const& t) const;

 private:
  std::vector<Scalar> const coefficients_;
  int const degree_;
  Instant const& t_min_;
  Instant const& t_max_;
  Time::Inverse inverse_duration_;
};

}  // namespace numerics
}  // namespace principia

#include "numerics/чебышёв_series_body.hpp"
