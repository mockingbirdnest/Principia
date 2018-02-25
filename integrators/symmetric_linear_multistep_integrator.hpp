
// The files containing the tree of of child classes of |Integrator| must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
//  parent.
#ifndef PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
#include "integrators/integrators.hpp"
#else
#ifndef PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#define PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_

#include <list>
#include <vector>

#include "base/status.hpp"
#include "integrators/backward_difference.hpp"
#include "integrators/cohen_hubbard_oesterwinter.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_symmetric_linear_multistep_integrator {

using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::DoublePrecision;
using numerics::FixedVector;
using quantities::Time;

template<typename Position, int order_>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>> {
  static constexpr int half_order_ = order_ / 2 + 1;
 public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;
  using AppendState = typename Integrator<ODE>::AppendState;

  static constexpr int order = order_;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    Status Solve(Instant const& t_final) override;
    SymmetricLinearMultistepIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;

   private:
    // The data for a previous step of the integration.  The |Displacement|s
    // here are really |Position|s, but we do complex computations on them and
    // it would be very inconvenient to cast these computations as barycentres.
    struct Step final {
      std::vector<DoublePrecision<typename ODE::Displacement>> displacements;
      std::vector<typename ODE::Acceleration> accelerations;
      DoublePrecision<Instant> time;

      void WriteToMessage(
          not_null<serialization::SymmetricLinearMultistepIntegratorInstance::
                       Step*> const message) const;
      static Step ReadFromMessage(
          serialization::SymmetricLinearMultistepIntegratorInstance::Step const&
              message);
    };

    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             SymmetricLinearMultistepIntegrator const& integrator);

    // For deserialization.
    Instance(IntegrationProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             int startup_step_index,
             std::list<Step> const& previous_steps,
             SymmetricLinearMultistepIntegrator const& integrator);

    // Performs the startup integration, i.e., computes enough states to either
    // reach |t_final| or to reach a point where |instance.previous_steps_| has
    // |order - 1| elements.  During startup |instance.current_state_| is
    // updated more frequently than once every |instance.step_|.
    void StartupSolve(Instant const& t_final);

    // Performs the velocity computation using a backward difference based on
    // the positions computed by the main integrator.
    void ComputeVelocityUsingBackwardDifferences();

    // Performs the velocity computation using the Cohen-Hubbard-Oesterwinter
    // method based on the accelerations computed by the main integrator.
    void ComputeVelocityUsingCohenHubbardOesterwinter();

    static void FillStepFromSystemState(ODE const& equation,
                                        typename ODE::SystemState const& state,
                                        Step& step);

    int startup_step_index_ = 0;
    std::list<Step> previous_steps_;  // At most |order_| elements.
    SymmetricLinearMultistepIntegrator const& integrator_;
    friend class SymmetricLinearMultistepIntegrator;
  };

  SymmetricLinearMultistepIntegrator(
      serialization::FixedStepSizeIntegrator::Kind kind,
      FixedStepSizeIntegrator<ODE> const& startup_integrator,
      FixedVector<double, half_order_> const& ɑ,
      FixedVector<double, half_order_> const& β_numerator,
      double β_denominator);

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

 private:
  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> ReadFromMessage(
      serialization::FixedStepSizeIntegratorInstance const& message,
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  BackwardDifference<order_ - 1> const& backward_difference_;
  CohenHubbardOesterwinter<order_> const& cohen_hubbard_oesterwinter_;
  FixedVector<double, half_order_> const ɑ_;
  FixedVector<double, half_order_> const β_numerator_;
  double const β_denominator_;
};

}  // namespace internal_symmetric_linear_multistep_integrator

using internal_symmetric_linear_multistep_integrator::
    SymmetricLinearMultistepIntegrator;

// This method and the next are from Quinlan (1999), Resonances and
// instabilities in symmetric multistep methods,
// https://arxiv.org/abs/astro-ph/9901136.
template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8A();

template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
Quinlan1999Order8B();

// These four methods are from Quinlan and Tremaine (1990), Symmetric multistep
// methods for the numerical integration of planetary orbits,
// http://adsabs.harvard.edu/full/1990AJ....100.1694Q.
template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/8> const&
QuinlanTremaine1990Order8();

template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/10> const&
QuinlanTremaine1990Order10();

template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/12> const&
QuinlanTremaine1990Order12();

template<typename Position>
SymmetricLinearMultistepIntegrator<Position,
                                   /*order=*/14> const&
QuinlanTremaine1990Order14();

}  // namespace integrators
}  // namespace principia

#include "symmetric_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
