
#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

namespace principia {
namespace integrators {

template<typename Position,
         typename Momentum,
         int order_,
         int evaluations_,
         bool time_reversible_,
         bool first_same_as_last_>
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          order_,
                                          evaluations_,
                                          time_reversible_,
                                          first_same_as_last_>::
    SymplecticPartitionedRungeKuttaIntegrator(
        FixedVector<double, stages_> const& a,
        FixedVector<double, stages_> const& b)
    : a_(a),
      b_(b) {
  if (first_same_as_last_) {
    CHECK_EQ(0.0, a_[stages_ - 1]);
  }
  if (time_reversible) {
    CHECK(first_same_as_last);
    for (int i = 0; i < stages_ - 1; ++i) {
      CHECK_EQ(a_[i], a_[stages_ - 2 - i]);
    }
    for (int i = 0; i < stages_; ++i) {
      CHECK_EQ(b_[i], b_[stages_ - 1 - i]);
    }
  }
}

template<typename Position,
         typename Momentum,
         int order_,
         int evaluations_,
         bool time_reversible_,
         bool first_same_as_last_>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      order_,
                                      time_reversible_,
                                      evaluations_,
                                      first_same_as_last_ ? BAB : BA> const&
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          order_,
                                          evaluations_,
                                          time_reversible_,
                                          first_same_as_last_>::BForceMethod()
    const {
  if (b_force_method_ == nullptr) {
    b_force_method_ = std::make_unique<
        SymplecticRungeKuttaNyströmIntegrator<Position,
                                              order_,
                                              time_reversible_,
                                              evaluations_,
                                              first_same_as_last_ ? BAB : BA>>(
        serialization::FixedStepSizeIntegrator::DUMMY,
        a_,
        b_);
  }
  return *b_force_method_;
}

template<typename Position,
         typename Momentum,
         int order_,
         int evaluations_,
         bool time_reversible_,
         bool first_same_as_last_>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      order_,
                                      time_reversible_,
                                      evaluations_,
                                      first_same_as_last_ ? ABA : BA> const&
SymplecticPartitionedRungeKuttaIntegrator<Position,
                                          Momentum,
                                          order_,
                                          evaluations_,
                                          time_reversible_,
                                          first_same_as_last_>::AForceMethod()
    const {
  if (a_force_method_ == nullptr) {
    auto shifted_a = a_;
    if (first_same_as_last) {
      // |*this| is a BAB method, with A and B interchangeable.  Exchanging A
      // and B does not only swap |a_| and |b_|, it shifts |a_| (because |ABA|
      // means b₀ vanishes, whereas |BAB| means aᵣ vanishes).
      for (int i = 0; i < stages_; ++i) {
        shifted_a[i] = a_[(i - 1) % stages_];
      }
    }
    a_force_method_ = std::make_unique<
        SymplecticRungeKuttaNyströmIntegrator<Position,
                                              order_,
                                              time_reversible_,
                                              evaluations_,
                                              first_same_as_last_ ? ABA : BA>>(
        serialization::FixedStepSizeIntegrator::DUMMY,
        b_,
        shifted_a);
  }
  return *a_force_method_;
}

}  // namespace integrators
}  // namespace principia
