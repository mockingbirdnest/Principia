#pragma once

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/synchronization/mutex.h"
#include "base/recurring_thread.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/repeated_field.h"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/checkpointer.hpp"
#include "physics/clientele.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/geopotential.hpp"
#include "physics/integration_parameters.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/protector.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "serialization/numerics.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_ephemeris {

using base::not_null;
using base::RecurringThread;
using geometry::InfinitePast;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::ExplicitSecondOrderOrdinaryDifferentialEquation;
using integrators::FixedStepSizeIntegrator;
using integrators::InitialValueProblem;
using integrators::Integrator;
using integrators::SpecialSecondOrderDifferentialEquation;
using quantities::Acceleration;
using quantities::Length;
using quantities::SpecificEnergy;
using quantities::Speed;
using quantities::Time;

// Note on thread-safety: the integration functions (Prolong, FlowWithFixedStep,
// FlowWithAdaptiveStep) can be called concurrently as long as their parameters
// designated distinct objects.  No guarantee is offered for the other
// functions.
template<typename Frame>
class Ephemeris {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  using IntrinsicAcceleration =
      std::function<Vector<Acceleration, Frame>(Instant const& time)>;
  static std::nullptr_t constexpr NoIntrinsicAcceleration = nullptr;
  using GeneralizedIntrinsicAcceleration =
      std::function<Vector<Acceleration, Frame>(
          Instant const& time,
          DegreesOfFreedom<Frame> const& degrees_of_freedom)>;
  using IntrinsicAccelerations = std::vector<IntrinsicAcceleration>;
  static IntrinsicAccelerations const NoIntrinsicAccelerations;
  static std::int64_t constexpr unlimited_max_ephemeris_steps =
      std::numeric_limits<std::int64_t>::max();

  // The equations describing the motion of the |bodies_|.
  using NewtonianMotionEquation =
      SpecialSecondOrderDifferentialEquation<Position<Frame>>;
  using GeneralizedNewtonianMotionEquation =
      ExplicitSecondOrderOrdinaryDifferentialEquation<Position<Frame>>;

  using AdaptiveStepParameters =
      physics::AdaptiveStepParameters<NewtonianMotionEquation>;
  using FixedStepParameters =
      physics::FixedStepParameters<NewtonianMotionEquation>;
  using GeneralizedAdaptiveStepParameters =
      physics::AdaptiveStepParameters<GeneralizedNewtonianMotionEquation>;

  class AccuracyParameters final {
   public:
    AccuracyParameters(Length const& fitting_tolerance,
                       double geopotential_tolerance);

    void WriteToMessage(
        not_null<serialization::Ephemeris::AccuracyParameters*> message) const;
    static AccuracyParameters ReadFromMessage(
        serialization::Ephemeris::AccuracyParameters const& message);

   private:
    Length fitting_tolerance_;
    double geopotential_tolerance_ = 0;
    friend class Ephemeris<Frame>;
  };

  // Constructs an Ephemeris that owns the |bodies|.  The elements of vectors
  // |bodies| and |initial_state| correspond to one another.
  Ephemeris(std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
            std::vector<DegreesOfFreedom<Frame>> const& initial_state,
            Instant const& initial_time,
            AccuracyParameters const& accuracy_parameters,
            FixedStepParameters fixed_step_parameters);

  virtual ~Ephemeris();

  // Returns the bodies in the order in which they were given at construction.
  virtual std::vector<not_null<MassiveBody const*>> const& bodies() const;

  // Returns the trajectory for the given |body|.
  virtual not_null<ContinuousTrajectory<Frame> const*> trajectory(
      not_null<MassiveBody const*> body) const;

  // Returns true if at least one of the trajectories is empty.
  virtual bool empty() const;

  // The maximum of the |t_min|s of the trajectories.
  virtual Instant t_min() const EXCLUDES(lock_);
  // The mimimum of the |t_max|s of the trajectories.
  virtual Instant t_max() const;

  virtual FixedStepSizeIntegrator<NewtonianMotionEquation> const&
  planetary_integrator() const;

  virtual absl::Status last_severe_integration_status() const;

  // Prolongs the ephemeris up to at least |t|.  Returns an error iff the thread
  // is stopped.  After a successful call, |t_max() >= t|.
  virtual absl::Status Prolong(Instant const& t) EXCLUDES(lock_);

  // Asks the reanimator thread to asynchronously reconstruct the past so that
  // the |t_min()| of the ephemeris ultimately ends up at or before
  // |desired_t_min|.
  void RequestReanimation(Instant const& desired_t_min);

  // Same as |RequestReanimation|, but synchronous.  This function blocks until
  // the |t_min()| of the ephemeris is at or before |desired_t_min|.
  void AwaitReanimation(Instant const& desired_t_min);

  // Creates an instance suitable for integrating the given |trajectories| with
  // their |intrinsic_accelerations| using a fixed-step integrator parameterized
  // by |parameters|.
  virtual not_null<
      std::unique_ptr<typename Integrator<NewtonianMotionEquation>::Instance>>
  NewInstance(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      IntrinsicAccelerations const& intrinsic_accelerations,
      FixedStepParameters const& parameters);

  // Same as above, but returns an error status if the thread is stopped.
  virtual absl::StatusOr<not_null<
      std::unique_ptr<typename Integrator<NewtonianMotionEquation>::Instance>>>
  StoppableNewInstance(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      IntrinsicAccelerations const& intrinsic_accelerations,
      FixedStepParameters const& parameters);

  // Integrates, until exactly |t| (except for timeouts or singularities), the
  // |trajectory| followed by a massless body in the gravitational potential
  // described by |*this|.  If |t > t_max()|, calls |Prolong(t)| beforehand.
  // Prolongs the ephemeris by at most |max_ephemeris_steps|.  Returns OK if and
  // only if |*trajectory| was integrated until |t|.
  virtual absl::Status FlowWithAdaptiveStep(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      IntrinsicAcceleration intrinsic_acceleration,
      Instant const& t,
      AdaptiveStepParameters const& parameters,
      std::int64_t max_ephemeris_steps) EXCLUDES(lock_);

  // Same as above, but uses a generalized integrator.
  virtual absl::Status FlowWithAdaptiveStep(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      GeneralizedIntrinsicAcceleration intrinsic_acceleration,
      Instant const& t,
      GeneralizedAdaptiveStepParameters const& parameters,
      std::int64_t max_ephemeris_steps) EXCLUDES(lock_);

  // Integrates, until at most |t|, the trajectories followed by massless
  // bodies in the gravitational potential described by |*this|.  If
  // |t > t_max()|, calls |Prolong(t)| beforehand.  The trajectories and
  // integration parameters are given by the |instance|.
  virtual absl::Status FlowWithFixedStep(
      Instant const& t,
      typename Integrator<NewtonianMotionEquation>::Instance& instance)
      EXCLUDES(lock_);

  // Returns the gravitational acceleration on a massless body located at the
  // given |position| at time |t|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      Position<Frame> const& position,
      Instant const& t) const EXCLUDES(lock_);

  // Returns the gravitational acceleration on the massless body having the
  // given |trajectory| at time |t|.  |t| must be one of the times of the
  // |trajectory|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      Instant const& t) const EXCLUDES(lock_);

  // Returns the gravitational acceleration on the massive |body| at time |t|.
  // |body| must be one of the bodies of this object.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMassiveBody(
      not_null<MassiveBody const*> body,
      Instant const& t) const EXCLUDES(lock_);

  // Returns the potential at the given |position| at time |t|.
  SpecificEnergy ComputeGravitationalPotential(
      Position<Frame> const& position,
      Instant const& t) const EXCLUDES(lock_);

  // Computes the apsides of the relative trajectory of |body1| and |body2}.
  // Appends to the given trajectories two points for each apsis, one for
  // |body1| and one for |body2|.  The times of |apoapsides1| and |apoapsideds2|
  // are identical (are similarly for |periapsides1| and |periapsides2|).
  virtual void ComputeApsides(not_null<MassiveBody const*> body1,
                              not_null<MassiveBody const*> body2,
                              DiscreteTrajectory<Frame>& apoapsides1,
                              DiscreteTrajectory<Frame>& periapsides1,
                              DiscreteTrajectory<Frame>& apoapsides2,
                              DiscreteTrajectory<Frame>& periapsides2);

  // Returns the index of the given body in the serialization produced by
  // |WriteToMessage| and read by the |Read...| functions.  This index is not
  // suitable for other uses.
  virtual int serialization_index_for_body(
      not_null<MassiveBody const*> body) const;

  virtual not_null<MassiveBody const*> body_for_serialization_index(
      int serialization_index) const;

  virtual void WriteToMessage(
      not_null<serialization::Ephemeris*> message) const EXCLUDES(lock_);
  template<typename F = Frame,
           typename = std::enable_if_t<base::is_serializable_v<F>>>
  // The parameter |desired_t_min| indicates that the ephemeris must be restored
  // at a checkpoint such that, once the ephemeris is prolonged, its |t_min()|
  // is at or before |desired_t_min|.
  static not_null<std::unique_ptr<Ephemeris>> ReadFromMessage(
      Instant const& desired_t_min,
      serialization::Ephemeris const& message) EXCLUDES(lock_);

 protected:
  // For mocking purposes, leaves everything uninitialized and uses the given
  // |integrator|.
  explicit Ephemeris(FixedStepSizeIntegrator<typename Ephemeris<
                         Frame>::NewtonianMotionEquation> const& integrator);

 private:
  // Checkpointing support.
  void WriteToCheckpointIfNeeded(Instant const& time) const
      SHARED_LOCKS_REQUIRED(lock_);
  Checkpointer<serialization::Ephemeris>::Writer MakeCheckpointerWriter();
  Checkpointer<serialization::Ephemeris>::Reader MakeCheckpointerReader();

  // Called on a stoppable thread to reconstruct the past state of the ephemeris
  // and its trajectories starting in such a way that |t_min()| is at or before
  // |desired_t_min|.  The member variable |oldest_reanimated_checkpoint_| tells
  // the reanimator where to stop.
  absl::Status Reanimate(Instant const desired_t_min) EXCLUDES(lock_);

  // Reconstructs the past state of the ephemeris between |t_initial| and
  // |t_final| using the given checkpoint |message|.
  absl::Status ReanimateOneCheckpoint(
      serialization::Ephemeris::Checkpoint const& message,
      Instant const& t_initial,
      Instant const& t_final) EXCLUDES(lock_);

  // Callbacks for the integrators.
  void AppendMassiveBodiesState(
      typename NewtonianMotionEquation::State const& state)
      REQUIRES(lock_);
  template<typename ContinuousTrajectoryPtr>
  static std::vector<absl::Status> AppendMassiveBodiesStateToTrajectories(
      typename NewtonianMotionEquation::State const& state,
      std::vector<not_null<ContinuousTrajectoryPtr>> const& trajectories);
  static void AppendMasslessBodiesStateToTrajectories(
      typename NewtonianMotionEquation::State const& state,
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories);

  // Returns an equation suitable for the massive bodies contained in this
  // ephemeris.
  NewtonianMotionEquation MakeMassiveBodiesNewtonianMotionEquation();

  // Note the return by copy: the returned value is usable even if the
  // |instance_| is being integrated.
  Instant instance_time() const EXCLUDES(lock_);

  virtual Instant t_min_locked() const REQUIRES_SHARED(lock_);
  virtual Instant t_max_locked() const REQUIRES_SHARED(lock_);

  // Computes the accelerations between one body, |body1| (with index |b1| in
  // the |positions| and |accelerations| arrays) and the bodies |bodies2| (with
  // indices [b2_begin, b2_end[ in the |bodies2|, |positions| and
  // |accelerations| arrays).  The template parameters specify what we know
  // about the bodies, and therefore what forces apply.  Works for both owning
  // and non-owning pointers thanks to the |MassiveBodyConstPtr| template
  // parameter.
  template<bool body1_is_oblate,
           bool body2_is_oblate,
           typename MassiveBodyConstPtr>
  static void ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
      Instant const& t,
      MassiveBody const& body1,
      std::size_t b1,
      std::vector<not_null<MassiveBodyConstPtr>> const& bodies2,
      std::size_t b2_begin,
      std::size_t b2_end,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations,
      std::vector<Geopotential<Frame>> const& geopotentials);

  // Computes the accelerations due to one body, |body1| (with index |b1| in the
  // |bodies_| and |trajectories_| arrays) on massless bodies at the given
  // |positions|.  The template parameter specifies what we know about the
  // massive body, and therefore what forces apply.  Returns an integer for
  // efficiency.
  template<bool body1_is_oblate>
  std::underlying_type_t<absl::StatusCode>
  ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
      Instant const& t,
      MassiveBody const& body1,
      std::size_t b1,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const
      REQUIRES_SHARED(lock_);

  // Computes the potential resulting from one body, |body1| (with index |b1| in
  // the |bodies_| and |trajectories_| arrays) at the given |positions|.  The
  // template parameter specifies what we know about the massive body, and
  // therefore what potential applies.
  template<bool body1_is_oblate>
  void ComputeGravitationalPotentialsOfMassiveBody(
      Instant const& t,
      MassiveBody const& body1,
      std::size_t b1,
      std::vector<Position<Frame>> const& positions,
      std::vector<SpecificEnergy>& potentials) const
      REQUIRES_SHARED(lock_);

  // Computes the accelerations between all the massive bodies in |bodies_|.
  absl::Status ComputeGravitationalAccelerationBetweenAllMassiveBodies(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const;

  // Computes the acceleration exerted by the massive bodies in |bodies_| on
  // massless bodies.  The massless bodies are at the given |positions|.
  // Returns an error iff a collision occurred, i.e., the massless body is
  // inside one of the |bodies_|.
  absl::StatusCode
  ComputeGravitationalAccelerationByAllMassiveBodiesOnMasslessBodies(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const
      EXCLUDES(lock_);

  // Computes the potential resulting from the massive bodies in |bodies_|.  The
  // potentials are computed at the given |positions|.
  void ComputeGravitationalPotentialsOfAllMassiveBodies(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<SpecificEnergy>& potentials) const
      EXCLUDES(lock_);

  // Flows the given ODE with an adaptive step integrator.
  template<typename ODE>
  absl::Status FlowODEWithAdaptiveStep(
      typename ODE::RightHandSideComputation compute_acceleration,
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      Instant const& t,
      physics::AdaptiveStepParameters<ODE> const& parameters,
      std::int64_t max_ephemeris_steps) EXCLUDES(lock_);

  // Computes an estimate of the ratio |tolerance / error|.
  static double ToleranceToErrorRatio(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance,
      Time const& current_step_size,
      typename NewtonianMotionEquation::State const& /*state*/,
      typename NewtonianMotionEquation::State::Error const& error);

  // The bodies in the order in which they were given at construction.
  std::vector<not_null<MassiveBody const*>> unowned_bodies_;

  // The indices of bodies in |unowned_bodies_|.
  std::map<not_null<MassiveBody const*>, int> unowned_bodies_indices_;

  // The oblate bodies precede the spherical bodies in this vector.  The system
  // state is indexed in the same order.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies_;

  // Only has entries for the oblate bodies, at the same indices as |bodies_|.
  std::vector<Geopotential<Frame>> geopotentials_;

  // The indices in |bodies_| correspond to those in |trajectories_|.
  std::vector<not_null<ContinuousTrajectory<Frame>*>> trajectories_;

  std::map<not_null<MassiveBody const*>,
           not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>>
      bodies_to_trajectories_;

  AccuracyParameters const accuracy_parameters_;
  FixedStepParameters const fixed_step_parameters_;

  int number_of_oblate_bodies_ = 0;
  int number_of_spherical_bodies_ = 0;

  not_null<
      std::unique_ptr<Checkpointer<serialization::Ephemeris>>> checkpointer_;

  // This member must only be accessed by the |reanimator_| thread, or before
  // the |reanimator_| thread is started.  An ephemeris that is constructed de
  // novo won't ever need reanimation, so all the checkpoints are animate at
  // birth.
  Instant oldest_reanimated_checkpoint_ = InfinitePast;

  // The techniques and terminology follow [Lov22].
  RecurringThread<Instant> reanimator_;
  Clientele<Instant> reanimator_clientele_;

  // The fields above this line are fixed at construction and therefore not
  // protected.  Note that |ContinuousTrajectory| is thread-safe.  |lock_| is
  // also used to protect sections where the trajectories are not mutually
  // consistent (e.g., during Prolong).
  mutable absl::Mutex lock_;

  // Parameter passed to the last call to |RequestReanimation|, if any.
  std::optional<Instant> last_desired_t_min_ GUARDED_BY(lock_);

  std::unique_ptr<typename Integrator<NewtonianMotionEquation>::Instance>
      instance_ GUARDED_BY(lock_);

  absl::Status last_severe_integration_status_ GUARDED_BY(lock_);
};

}  // namespace internal_ephemeris

using internal_ephemeris::Ephemeris;

}  // namespace physics
}  // namespace principia

#include "physics/ephemeris_body.hpp"
