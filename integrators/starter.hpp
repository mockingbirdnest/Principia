#pragma once

#include <list>

#include "base/not_null.hpp"
#include "integrators/integrators.hpp"

namespace principia {
namespace integrators {
namespace internal_starter {

using namespace principia::base::_not_null;

// A helper object for starting a linear multistep integrator.  |ODE| is the
// equation being integrated.  |Step| is an object holding the data produced by
// executing a step of the integrator and used for subsequent steps.  The
// function |FillStepFromState| must fill this object from an |ODE::State|.
// |steps| is the number of steps of the integrator being started.
template<typename ODE, typename Step, int steps>
class Starter {
 public:
  // |startup_step_divisor| specifies how many steps of the |startup_integrator|
  // should be produced for a step of the main integrator.  |instance| is the
  // instance of the main integrator.
  Starter(FixedStepSizeIntegrator<ODE> const& startup_integrator,
          std::int64_t startup_step_divisor,
          not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const
              instance);

  // Performs the startup integration, i.e., computes enough states to either
  // reach |s_final| or to reach a point where |instance.previous_steps_| has
  // |order - 1| elements.  During startup |instance.current_state_| is
  // updated more frequently than once every |instance.step_|.
  void Solve(typename ODE::IndependentVariable const& s_final);

  // Appends a new step to |previous_steps_| and drops the oldest one if needed.
  // This object must be |started()|.
  void Push(Step step);

  // Returns the startup steps.  This object must be |started()|.
  std::list<Step> const& previous_steps() const;

  // Returns true iff the startup steps have all been computed.
  bool started() const;

  // Serialization helpers to write/read the starter data to/from a message.
  template<typename Message>
  void WriteToMessage(not_null<Message*> message) const;
  template<typename Message>
  void FillFromMessage(Message const& message);

 protected:
  // Must fill |step| from |state| and the right-hand side of the ODE.  |Step|
  // must contain all the information needed to compute subsequent steps of the
  // integrator.
  virtual void FillStepFromState(ODE const& equation,
                                 typename ODE::State const& state,
                                 Step& step) const = 0;

  // Returns the last value of the independent variable known to the instance.
  virtual typename ODE::IndependentVariable independent_variable() const = 0;

  // A helper to implement the previous function in subclasses.
  typename FixedStepSizeIntegrator<ODE>::Instance const& instance() const;

 private:
  FixedStepSizeIntegrator<ODE> const& startup_integrator_;
  std::int64_t const startup_step_divisor_;
  not_null<typename FixedStepSizeIntegrator<ODE>::Instance*> const
      instance_;

  int startup_step_index_ = 0;
  std::list<Step> previous_steps_;  // At most |order| elements.
};

}  // namespace internal_starter

using internal_starter::Starter;

}  // namespace integrators
}  // namespace principia

#include "integrators/starter_body.hpp"
