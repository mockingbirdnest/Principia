
#pragma once

#include <list>
#include <map>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "integrators/integrators.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE_FROM(part, class, Part);

namespace internal_pile_up {

using base::not_null;
using geometry::Frame;
using geometry::Instant;
using geometry::Vector;
using integrators::Integrator;
using physics::DiscreteTrajectory;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::MasslessBody;
using physics::RelativeDegreesOfFreedom;
using quantities::Force;
using quantities::Mass;

// A |PileUp| handles a connected component of the graph of |Parts| under
// physical contact.  It advances the history and prolongation of its component
// |Parts|, modeling them as a massless body at their centre of mass.
class PileUp {
 public:
  PileUp(
      std::list<not_null<Part*>>&& parts,
      Instant const& t,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
      not_null<Ephemeris<Barycentric>*> ephemeris);

  virtual ~PileUp() = default;

  // This class is moveable.
  PileUp(PileUp&& pile_up) = default;
  PileUp& operator=(PileUp&& pile_up) = default;

  void set_mass(Mass const& mass);
  void set_intrinsic_force(Vector<Force, Barycentric> const& intrinsic_force);

  std::list<not_null<Part*>> const& parts() const;

  // Set the |degrees_of_freedom| for the given |part|.  These degrees of
  // freedom are *apparent* in the sense that they were reported by the game but
  // we know better since we are doing science.
  void SetPartApparentDegreesOfFreedom(
      not_null<Part*> part,
      DegreesOfFreedom<ApparentBubble> const& degrees_of_freedom);

  // Update the degrees of freedom (in |RigidPileUp|) of all the parts by
  // translating the *apparent* degrees of freedom so that their centre of mass
  // matches that computed by integration.
  // |SetPartApparentDegreesOfFreedom| must have been called for each part in
  // the pile-up, or for none.
  // The degrees of freedom set by this method are used by |NudgeParts|.
  // NOTE(egg): Eventually, this will also change their velocities and angular
  // velocities so that the angular momentum matches that which has been
  // computed for |this| |PileUp|.
  void DeformPileUpIfNeeded();

  // Flows the history authoritatively as far as possible up to |t|, advances
  // the histories of the parts and updates the degrees of freedom of the parts
  // if the pile-up is in the bubble.  After this call, the tail (of |*this|)
  // and of its parts have a (possibly ahistorical) final point exactly at |t|.
  void AdvanceTime(Instant const& t);

  // Adjusts the degrees of freedom of all parts in this pile up based on the
  // degrees of freedom of the pile-up computed by |AdvanceTime| and on the
  // |RigidPileUp| degrees of freedom of the parts, as set by
  // |DeformPileUpIfNeeded|.
  void NudgeParts() const;

  // Returns the last |Instant| passed to |AdvanceTime|.  Note that this may
  // correspond to a non-authoritative point.
  Instant const& time();

  void WriteToMessage(not_null<serialization::PileUp*> message) const;
  static PileUp ReadFromMessage(
      serialization::PileUp const& message,
      std::function<not_null<Part*>(PartId)> const& part_id_to_part,
      not_null<Ephemeris<Barycentric>*> ephemeris);

 private:
  // For deserialization.
  PileUp(
      std::list<not_null<Part*>>&& parts,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
      not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> psychohistory,
      not_null<Ephemeris<Barycentric>*> ephemeris);

  void AppendToPartTails(DiscreteTrajectory<Barycentric>::Iterator it,
                         bool authoritative) const;

  std::list<not_null<Part*>> parts_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::FixedStepParameters fixed_step_parameters_;

  // An optimization: the sum of the masses and intrinsic forces of the
  // |parts_|, computed by the union-find.
  Mass mass_;
  Vector<Force, Barycentric> intrinsic_force_;

  // TODO(phl): replace by an instance.  Specifically, this should contain
  // either an adaptive step instance (if the last call to AdvanceTime occurred
  // with nonzero intrinsic force), or a fixed step instance otherwise (with the
  // prolongation being computed by an instance local to the body of
  // |AdvanceTime|).
  // |psychohistory_.Size()| is either 1 or 2.  The first point is
  // authoritative, and the second point, if any, is not.
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> psychohistory_;

  // When present, this instance is used to integrate the trajectory of this
  // pile-up using a fixed-step integrator.  This instance is destroyed
  // if a variable-step integrator needs to be used because of an intrinsic
  // acceleration.
  std::unique_ptr<typename Integrator<
      Ephemeris<Barycentric>::NewtonianMotionEquation>::Instance>
      fixed_instance_;

  // The |PileUp| is seen as a (currently non-rotating) rigid body; the degrees
  // of freedom of the parts in the frame of that body can be set, however their
  // motion is not integrated; this is simply applied as an offset from the
  // rigid body motion of the |PileUp|.
  // The origin of |RigidPileUp| is the centre of mass of the pile up.
  // Its axes are those of Barycentric for now; eventually we will probably want
  // to use the inertia ellipsoid.
  using RigidPileUp = Frame<serialization::Frame::PluginTag,
                            serialization::Frame::RIGID_PILE_UP,
                            /*frame_is_inertial=*/false>;

  PartTo<DegreesOfFreedom<RigidPileUp>> actual_part_degrees_of_freedom_;
  PartTo<DegreesOfFreedom<ApparentBubble>> apparent_part_degrees_of_freedom_;

  friend class TestablePileUp;
};

}  // namespace internal_pile_up

using internal_pile_up::PileUp;

}  // namespace ksp_plugin
}  // namespace principia
