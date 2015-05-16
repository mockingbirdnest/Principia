
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
      t_max_(t_max) {
  CHECK_LE(0, degree_) << "Degree must be at least 0";
  CHECK_LT(t_min_, t_max_) << "Time interval must not be empty";
  // Precomputed to save operations at the expense of some accuracy loss.
  Time const duration = t_max_ - t_min_;
  t_mean_ = t_min_ + 0.5 * duration;
  two_over_duration_ = 2 / duration;
}

template<typename Scalar>
Scalar ЧебышёвSeries<Scalar>::Evaluate(Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  double const two_scaled_t = scaled_t + scaled_t;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);

  double b_kplus2 = 0.0;
  double b_kplus1 = 0.0;
  double b_k = 0.0;
  for (int k = degree_; k >= 1; --k) {
    b_k = coefficients_[k] + two_scaled_t * b_kplus1 - b_kplus2;
    b_kplus2 = b_kplus1;
    b_kplus1 = b_k;
  }
  return coefficients_[0] + scaled_t * b_kplus1 - b_kplus2;
}

template<typename Scalar>
void ЧебышёвSeries<Scalar>::WriteToMessage(
    not_null<serialization::ЧебышёвSeries*> const message) const {
  for (auto const& coefficient : coefficients_) {
    coefficient.WriteToMessage(message->add_coefficient());
  }
  t_min_.WriteToMessage(message->mutable_t_min());
  t_max_.WriteToMessage(message->mutable_t_max());
}

template<typename Scalar>
static ЧебышёвSeries<Scalar> ЧебышёвSeries<Scalar>::ReadFromMessage(
    serialization::ЧебышёвSeries const& message) {
}

}  // namespace numerics
}  // namespace principia
