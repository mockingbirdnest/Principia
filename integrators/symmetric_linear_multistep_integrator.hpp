
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

template<typename Method, typename Position>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>> {
 public:
  using ODE = SpecialSecondOrderDifferentialEquation<Position>;
  using AppendState = typename Integrator<ODE>::AppendState;

  static constexpr int order = Method::order;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    Status Solve(Instant const& t_final) override;
    SymmetricLinearMultistepIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    template<typename P = Position,
             typename = std::enable_if_t<base::is_serializable_v<P>>>
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::SymmetricLinearMultistepIntegratorInstance const&
            extension,
        IntegrationProblem<ODE> const& problem,
        AppendState const& append_state,
        Time const& step,
        SymmetricLinearMultistepIntegrator const& integrator);

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
                       Step*> message) const;
      template<typename P = Position,
               typename = std::enable_if_t<base::is_serializable_v<P>>>
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

  explicit SymmetricLinearMultistepIntegrator(
      FixedStepSizeIntegrator<ODE> const& startup_integrator);

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      IntegrationProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  static constexpr auto half_order_ = Method::Half(order);
  static constexpr auto ɑ_ = Method::ɑ;
  static constexpr auto β_numerator_ = Method::β_numerator;
  static constexpr auto β_denominator_ = Method::β_denominator;

  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  CohenHubbardOesterwinter<order> const& cohen_hubbard_oesterwinter_;
};

}  // namespace internal_symmetric_linear_multistep_integrator

template<typename Method, typename Position>
internal_symmetric_linear_multistep_integrator::
    SymmetricLinearMultistepIntegrator<Method, Position> const&
SymmetricLinearMultistepIntegrator();

}  // namespace integrators
}  // namespace principia

#include "symmetric_linear_multistep_integrator_body.hpp"

#endif  // PRINCIPIA_INTEGRATORS_SYMMETRIC_LINEAR_MULTISTEP_INTEGRATOR_HPP_
#endif  // PRINCIPIA_INTEGRATORS_INTEGRATORS_HPP_
