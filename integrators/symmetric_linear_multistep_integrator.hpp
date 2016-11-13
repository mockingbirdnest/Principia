// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
#include "integrators/ordinary_differential_equations.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include "base/status.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {

using numerics::FixedVector;

namespace integrators {

template <typename Position, int order_>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>> {
  static int const half_order_ = order_ / 2 + 1;
public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;

  SymmetricLinearMultistepIntegrator(
      serialization::FixedStepSizeIntegrator::Kind const kind,
      FixedVector<double, half_order_> const& ɑ,
      FixedVector<double, half_order_> const& β);

  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override;

  static int const order = order_;

 private:
  FixedVector<double, half_order_> const ɑ_;
  FixedVector<double, half_order_> const β_;
};

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
QuinlanTremaine1990Order8();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8A();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8B();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/10> const&
QuinlanTremaine1990Order10();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/12> const&
QuinlanTremaine1990Order12();

template <typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/14> const&
QuinlanTremaine1990Order14();

}  // namespace integrators
}  // namespace principia

#include "symmetric_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_ORDINARY_DIFFERENTIAL_EQUATIONS_HPP_
