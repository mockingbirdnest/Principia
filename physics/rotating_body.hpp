
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.  This class will end
// up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#include "physics/massive_body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#define PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_

#include "base/macros.hpp"
#include OPTIONAL_HEADER

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::AngularVelocity;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Angle;
using quantities::AngularFrequency;

namespace physics {

template<typename Frame>
class RotatingBody : public MassiveBody {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  class Parameters {
   public:
    // |reference_angle| is the angle of the prime meridian at
    // |reference_instant|.  |angular_velocity| gives the direction and speed of
    // the rotation of the body.
    Parameters(Angle const& reference_angle,
               Instant const& reference_instant,
               AngularVelocity<Frame> const& angular_velocity);

   private:
    Angle const reference_angle_;
    Instant const reference_instant_;
    AngularVelocity<Frame> const angular_velocity_;
    template<typename F>
    friend class RotatingBody;
  };

  RotatingBody(MassiveBody::Parameters const& massive_body_parameters,
               Parameters const& parameters);
  ~RotatingBody() = default;

  // Returns the angular velocity passed at construction.
  AngularVelocity<Frame> const& angular_velocity() const;

  // Returns the position at time |t|.
  Angle AngleAt(Instant const& t) const;

  // Returns the rotation at time |t|.
  Rotation<Frame, Frame> RotationAt(Instant const& t) const;

  // Returns false.
  bool is_massless() const override;

  // Returns true.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const override;

  // Fails if the |RotatingBody| extension is absent from the message.
  static not_null<std::unique_ptr<RotatingBody<Frame>>> ReadFromMessage(
      serialization::RotatingBody const& message,
      MassiveBody::Parameters const& massive_body_parameters);

 private:
  Parameters const parameters_;
};

}  // namespace physics
}  // namespace principia

#include "physics/rotating_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
