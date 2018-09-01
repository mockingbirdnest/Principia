
#pragma once

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "base/shared_lock_guard.hpp"
#include "base/status.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/repeated_field.h"
#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/geopotential.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_ephemeris {

using base::not_null;
using base::Status;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::FixedStepSizeIntegrator;
using integrators::Integrator;
using integrators::IntegrationProblem;
using integrators::SpecialSecondOrderDifferentialEquation;
using quantities::Acceleration;
using quantities::Length;
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
  using IntrinsicAccelerations = std::vector<IntrinsicAcceleration>;
  static IntrinsicAccelerations const NoIntrinsicAccelerations;
  static std::int64_t constexpr unlimited_max_ephemeris_steps =
      std::numeric_limits<std::int64_t>::max();

  // The equation describing the motion of the |bodies_|.
  using NewtonianMotionEquation =
      SpecialSecondOrderDifferentialEquation<Position<Frame>>;

  class PHYSICS_DLL AdaptiveStepParameters final {
   public:
    // The |length_| and |speed_integration_tolerance|s are used to compute the
    // |tolerance_to_error_ratio| for step size control.  The number of steps is
    // limited to |max_steps|.
    AdaptiveStepParameters(
        AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
        std::int64_t max_steps,
        Length const& length_integration_tolerance,
        Speed const& speed_integration_tolerance);

    AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const&
    integrator() const;
    std::int64_t max_steps() const;
    Length length_integration_tolerance() const;
    Speed speed_integration_tolerance() const;

    void set_max_steps(std::int64_t max_steps);
    void set_length_integration_tolerance(
        Length const& length_integration_tolerance);
    void set_speed_integration_tolerance(
        Speed const& speed_integration_tolerance);

    void WriteToMessage(
        not_null<serialization::Ephemeris::AdaptiveStepParameters*> const
            message) const;
    static AdaptiveStepParameters ReadFromMessage(
        serialization::Ephemeris::AdaptiveStepParameters const& message);

   private:
    // This will refer to a static object returned by a factory.
    not_null<AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const*>
        integrator_;
    std::int64_t max_steps_;
    Length length_integration_tolerance_;
    Speed speed_integration_tolerance_;
    friend class Ephemeris<Frame>;
  };

  class PHYSICS_DLL FixedStepParameters final {
   public:
    FixedStepParameters(
        FixedStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
        Time const& step);

    Time const& step() const;

    void WriteToMessage(
        not_null<serialization::Ephemeris::FixedStepParameters*> message) const;
    static FixedStepParameters ReadFromMessage(
        serialization::Ephemeris::FixedStepParameters const& message);

   private:
    // This will refer to a static object returned by a factory.
    not_null<FixedStepSizeIntegrator<NewtonianMotionEquation> const*>
        integrator_;
    Time step_;
    friend class Ephemeris<Frame>;
  };

  // Constructs an Ephemeris that owns the |bodies|.  The elements of vectors
  // |bodies| and |initial_state| correspond to one another.
  Ephemeris(std::vector<not_null<std::unique_ptr<MassiveBody const>>>&& bodies,
            std::vector<DegreesOfFreedom<Frame>> const& initial_state,
            Instant const& initial_time,
            Length const& fitting_tolerance,
            FixedStepParameters const& parameters);

  virtual ~Ephemeris() = default;

  // Returns the bodies in the order in which they were given at construction.
  virtual std::vector<not_null<MassiveBody const*>> const& bodies() const;

  // Returns the trajectory for the given |body|.
  virtual not_null<ContinuousTrajectory<Frame> const*> trajectory(
      not_null<MassiveBody const*> body) const;

  // Returns true if at least one of the trajectories is empty.
  virtual bool empty() const EXCLUDES(lock_);

  // The maximum of the |t_min|s of the trajectories.
  virtual Instant t_min() const EXCLUDES(lock_);
  // The mimimum of the |t_max|s of the trajectories.
  virtual Instant t_max() const EXCLUDES(lock_);

  virtual FixedStepSizeIntegrator<NewtonianMotionEquation> const&
  planetary_integrator() const;

  virtual Status last_severe_integration_status() const;

  // Calls |ForgetBefore| on all trajectories.  On return |t_min() == t|.
  virtual void ForgetBefore(Instant const& t);

  // Prolongs the ephemeris up to at least |t|.  After the call, |t_max() >= t|.
  virtual void Prolong(Instant const& t) EXCLUDES(lock_);

  // Creates an instance suitable for integrating the given |trajectories| with
  // their |intrinsic_accelerations| using a fixed-step integrator parameterized
  // by |parameters|.
  virtual not_null<
      std::unique_ptr<typename Integrator<NewtonianMotionEquation>::Instance>>
  NewInstance(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      IntrinsicAccelerations const& intrinsic_accelerations,
      FixedStepParameters const& parameters);

  // Integrates, until exactly |t| (except for timeouts or singularities), the
  // |trajectory| followed by a massless body in the gravitational potential
  // described by |*this|.  If |t > t_max()|, calls |Prolong(t)| beforehand.
  // Prolongs the ephemeris by at most |max_ephemeris_steps|.  If
  // |last_point_only| is true, only the last point is appended to the
  // trajectory.  Returns OK if and only if |*trajectory| was integrated until
  // |t|.
  virtual Status FlowWithAdaptiveStep(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      IntrinsicAcceleration intrinsic_acceleration,
      Instant const& t,
      AdaptiveStepParameters const& parameters,
      std::int64_t max_ephemeris_steps,
      bool last_point_only);

  // Integrates, until at most |t|, the trajectories followed by massless
  // bodies in the gravitational potential described by |*this|.  If
  // |t > t_max()|, calls |Prolong(t)| beforehand.  The trajectories and
  // integration parameters are given by the |instance|.
  virtual Status FlowWithFixedStep(
      Instant const& t,
      typename Integrator<NewtonianMotionEquation>::Instance& instance);

  // Returns the gravitational acceleration on a massless body located at the
  // given |position| at time |t|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      Position<Frame> const& position,
      Instant const& t) const;

  // Returns the gravitational acceleration on the massless body having the
  // given |trajectory| at time |t|.  |t| must be one of the times of the
  // |trajectory|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      not_null<DiscreteTrajectory<Frame>*> trajectory,
      Instant const& t) const;

  // Returns the gravitational acceleration on the massive |body| at time |t|.
  // |body| must be one of the bodies of this object.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMassiveBody(
      not_null<MassiveBody const*> body,
      Instant const& t) const REQUIRES_SHARED(lock_);

  // Computes the apsides of the relative trajectory of |body1| and |body2}.
  // Appends to the given trajectories two point for each apsis, one for |body1|
  // and one for |body2|.  The times of |apoapsides1| and |apoapsideds2| are
  // identical (are similarly for |periapsides1| and |periapsides2|).
  virtual void ComputeApsides(not_null<MassiveBody const*> const body1,
                              not_null<MassiveBody const*> const body2,
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
      not_null<serialization::Ephemeris*> message) const;
  static not_null<std::unique_ptr<Ephemeris>> ReadFromMessage(
      serialization::Ephemeris const& message);

 protected:
  // For mocking purposes, leaves everything uninitialized and uses the given
  // |integrator|.
  explicit Ephemeris(FixedStepSizeIntegrator<typename Ephemeris<
                         Frame>::NewtonianMotionEquation> const& integrator);

 private:
  // The state of the integration and of the continuous trajectory at a
  // particular time that we might want to use for compact serialization.
  struct Checkpoint final {
    std::unique_ptr<
        typename Integrator<NewtonianMotionEquation>::Instance> instance;
    std::vector<typename ContinuousTrajectory<Frame>::Checkpoint> checkpoints;
  };

  void AppendMassiveBodiesState(
      typename NewtonianMotionEquation::SystemState const& state)
      REQUIRES(lock_);
  static void AppendMasslessBodiesState(
      typename NewtonianMotionEquation::SystemState const& state,
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories);

  Checkpoint GetCheckpoint() REQUIRES_SHARED(lock_);

  // Same as t_max, but |lock_| must be held.
  Instant t_max_locked() const REQUIRES_SHARED(lock_);

  // Note the return by copy: the returned value is usable even if the
  // |instance_| is being integrated.
  Instant instance_time() const EXCLUDES(lock_);

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
      std::size_t const b1,
      std::vector<not_null<MassiveBodyConstPtr>> const& bodies2,
      std::size_t const b2_begin,
      std::size_t const b2_end,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations,
      std::vector<Geopotential<Frame>> const& geopotentials);

  // Computes the accelerations due to one body, |body1| (with index |b1| in the
  // |bodies_| and |trajectories_| arrays) on massless bodies at the given
  // |positions|.  The template parameter specifies what we know about the
  // massive body, and therefore what forces apply.  Returns false iff a
  // collision occurred, i.e., the massless body is inside |body1|.
  template<bool body1_is_oblate>
  bool ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
      Instant const& t,
      MassiveBody const& body1,
      std::size_t const b1,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const
      REQUIRES_SHARED(lock_);

  // Computes the accelerations between all the massive bodies in |bodies_|.
  void ComputeMassiveBodiesGravitationalAccelerations(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const;

  // Computes the acceleration exerted by the massive bodies in |bodies_| on
  // massless bodies.  The massless bodies are at the given |positions|.
  // Returns false iff a collision occurred, i.e., the massless body is inside
  // one of the |bodies_|.
  bool ComputeMasslessBodiesGravitationalAccelerations(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const
      EXCLUDES(lock_);

  // Same as above, but the massless bodies have intrinsic accelerations.
  // |intrinsic_accelerations| may be empty.  Returns false iff a collision
  // occurred, i.e., the massless body is inside one of the |bodies_|.
  bool ComputeMasslessBodiesTotalAccelerations(
      std::vector<IntrinsicAcceleration> const& intrinsic_accelerations,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      std::vector<Vector<Acceleration, Frame>>& accelerations) const;

  // Computes an estimate of the ratio |tolerance / error|.
  static double ToleranceToErrorRatio(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance,
      Time const& current_step_size,
      typename NewtonianMotionEquation::SystemStateError const& error);

  // Guards |instance_|, |trajectories_|, and |bodies_to_trajectories_| during
  // integration.  Note that the thread-safety annotations are incomplete
  // because we do not attempt to protect all the operations, only integration.
  mutable base::shared_mutex lock_;

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

  FixedStepParameters const parameters_;
  Length const fitting_tolerance_;
  std::unique_ptr<
      typename Integrator<NewtonianMotionEquation>::Instance> instance_;

  // These are the states other that the last which we preserve in order to
  // implement compact serialization.  The vector is time-ordered.
  std::vector<Checkpoint> checkpoints_;

  int number_of_oblate_bodies_ = 0;
  int number_of_spherical_bodies_ = 0;

  Status last_severe_integration_status_;
};

}  // namespace internal_ephemeris

using internal_ephemeris::Ephemeris;

}  // namespace physics
}  // namespace principia

#if !PHYSICS_DLL_IMPORT
#include "physics/ephemeris_body.hpp"
#endif
