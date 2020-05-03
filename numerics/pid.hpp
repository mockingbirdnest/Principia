
#pragma once

#include <deque>
#include <optional>

#include "geometry/named_quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_pid {

using geometry::Instant;
using geometry::Normed;
using quantities::Difference;
using quantities::Product;
using quantities::Quotient;
using quantities::Time;
using quantities::Variation;

template<typename Input,
         typename Output,
         int horizon,
         int finite_difference_order>
class PID {
 public:
  static_assert(finite_difference_order <= horizon);

  using Kp = Quotient<typename Normed<Output>::NormType,
                      typename Normed<Difference<Input>>::NormType>;
  using Ki = Quotient<typename Normed<Output>::NormType,
                      typename Normed<Product<Difference<Input>, Time>>::
                          NormType>;
  using Kd = Quotient<typename Normed<Output>::NormType,
                      typename Normed<Difference<Variation<Input>>>::NormType>;

  PID(Kp const& kp, Ki const& ki, Kd const& kd);

  // Clears the state of the PID.
  void Clear();

  // Adds the error between the process variable and the set-point to the state
  // of the PID and returns the control variable.
  Output ComputeControlVariable(Input const& process_variable,
                                Input const& set_point,
                                Instant const& t);

 private:
  Kp const kp_;
  Ki const ki_;
  Kd const kd_;

  // The PID is not prepared to receive inputs that are not equally spaced.
  std::optional<Time> previous_Δt_;
  std::optional<Instant> previous_t_;

  // The front element is the oldest.
  std::deque<Difference<Input>> errors_;
};

}  // namespace internal_pid

using internal_pid::PID;

}  // namespace numerics
}  // namespace principia

#include "numerics/pid_body.hpp"
