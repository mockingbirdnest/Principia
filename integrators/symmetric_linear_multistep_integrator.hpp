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

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "integrators/cohen_hubbard_oesterwinter.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "integrators/starter.hpp"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_symmetric_linear_multistep_integrator {

using numerics::DoublePrecision;
using numerics::FixedVector;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_named_quantities;
using namespace principia::quantities::_quantities;

// This implementation follows [QT90].
template<typename Method, typename ODE_>
class SymmetricLinearMultistepIntegrator
    : public FixedStepSizeIntegrator<ODE_> {
 public:
  using ODE = ODE_;
  static_assert(is_instance_of_v<SpecialSecondOrderDifferentialEquation, ODE>);
  using AppendState = typename Integrator<ODE>::AppendState;

  static constexpr int order = Method::order;

  class Instance : public FixedStepSizeIntegrator<ODE>::Instance {
   public:
    absl::Status Solve(Instant const& t_final) override;
    SymmetricLinearMultistepIntegrator const& integrator() const override;
    not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> Clone()
        const override;

    void WriteToMessage(
        not_null<serialization::IntegratorInstance*> message) const override;
    template<typename DV = typename ODE::DependentVariable,
             typename = std::enable_if_t<base::is_serializable_v<DV>>>
    static not_null<std::unique_ptr<Instance>> ReadFromMessage(
        serialization::SymmetricLinearMultistepIntegratorInstance const&
            extension,
        InitialValueProblem<ODE> const& problem,
        AppendState const& append_state,
        Time const& step,
        SymmetricLinearMultistepIntegrator const& integrator);

   private:
    // The data for a previous step of the integration.  The |Displacement|s
    // here are really |Position|s, but we do complex computations on them and
    // it would be very inconvenient to cast these computations as barycentres.
    struct Step final {
      std::vector<DoublePrecision<typename ODE::DependentVariableDifference>>
          displacements;
      typename ODE::DependentVariableDerivatives2 accelerations;
      DoublePrecision<Instant> time;

      void WriteToMessage(
          not_null<serialization::SymmetricLinearMultistepIntegratorInstance::
                       Step*> message) const;
      template<typename DV = typename ODE::DependentVariable,
               typename = std::enable_if_t<base::is_serializable_v<DV>>>
      static Step ReadFromMessage(
          serialization::SymmetricLinearMultistepIntegratorInstance::Step const&
              message);
    };

    class Starter : public integrators::Starter<ODE, Step, /*steps=*/order> {
     protected:
      using integrators::Starter<ODE, Step, order>::Starter;

      void FillStepFromState(ODE const& equation,
                             typename ODE::State const& state,
                             Step& step) const override;
      typename ODE::IndependentVariable independent_variable() const override;
    };

    Instance(InitialValueProblem<ODE> const& problem,
             AppendState const& append_state,
             Time const& step,
             SymmetricLinearMultistepIntegrator const& integrator);

    // Performs the velocity computation using the Cohen-Hubbard-Oesterwinter
    // method based on the accelerations computed by the main integrator.
    void ComputeVelocityUsingCohenHubbardOesterwinter();

    Starter starter_;
    SymmetricLinearMultistepIntegrator const& integrator_;
    friend class SymmetricLinearMultistepIntegrator;
  };

  explicit SymmetricLinearMultistepIntegrator(
      FixedStepSizeIntegrator<ODE> const& startup_integrator);

  not_null<std::unique_ptr<typename Integrator<ODE>::Instance>> NewInstance(
      InitialValueProblem<ODE> const& problem,
      AppendState const& append_state,
      Time const& step) const override;

  void WriteToMessage(
      not_null<serialization::FixedStepSizeIntegrator*> message) const override;

 private:
  static constexpr auto half_order_ = Method::Half(order);
  static constexpr auto α_ = Method::α;
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
