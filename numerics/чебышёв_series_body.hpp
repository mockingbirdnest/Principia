
#include "numerics/чебышёв_series.hpp"

#include "glog/logging.h"

namespace principia {
namespace numerics {

template<typename Scalar>
ЧебышёвSeries<Scalar>::ЧебышёвSeries(std::vector<Scalar> const& coefficients,
                                     Instant const& t_min,
                                     Instant const& t_max)
    : coefficients_(coefficients),
      degree_(static_cast<int>(coefficients_.size()) - 1),
      t_min_(t_min),
      t_max_(t_max_) {
  CHECK_LE(0, degree_);
  CHECK_LT(t_min_, t_max_);
  // Precomputed to save the division at the expense of some accuracy loss.
  inverse_duration_ = 1 / (t_max_ - t_min_);
}

template<typename Scalar>
Scalar ЧебышёвSeries<Scalar>::Evaluate(Instant const& t) const {
  double const scaled_t = (t - t_min_) * inverse_duration_;

  double b_k = coefficients_[degree_];
  double b_kplus1 = 0.0;
  double b_kplus2 = 0.0;
  for (int i = degree_ - 1; i >= 1; --i) {
    b_kplus2 = b_kplus1;
    b_kplus1 = b_k;
    b_k = coefficients_[i] + 2 * scaled_t * b_kplus1 - b_kplus2;
  }
  return coefficients_[0] + scaled_t * b_kplus1 - b_kplus2;
}

}  // namespace numerics
}  // namespace principia
