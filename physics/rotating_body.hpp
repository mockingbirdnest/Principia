
// The files containing the tree of of child classes of |Body| must be included
// in the order of inheritance to avoid circular dependencies.  This class will
// end up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#include "physics/massive_body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#define PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_

#include <optional.hpp>

#include <vector>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::Instant;
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
    Parameters(Angle const& reference_angle,
               Instant const& reference_time,
               AngularFrequency const& ω,
               Vector<double, Frame> const& axis);

   private:
    Angle const reference_angle_;
    Instant const reference_time_;
    AngularFrequency const ω_;
    Vector<double, Frame> const axis_;
    friend class RotatingBody;
  };

  RotatingBody(MassiveBody::Parameters const& massive_body_parameters,
               Parameters const& parameters);
  ~RotatingBody() = default;

  // Returns the axis passed at construction.
  Vector<double, Frame> const& axis() const;

  // Returns false.
  bool is_massless() const override;

  // Returns true.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const override;

  // Fails unless |message.has_massive_body()|.
  static not_null<std::unique_ptr<RotatingBody<Frame>>> ReadFromMessage(
      serialization::Body const& message);

  // Fails if the |RotatingBody| extension is absent from the message.
  static not_null<std::unique_ptr<RotatingBody<Frame>>> ReadFromMessage(
      serialization::MassiveBody const& message);

 private:
  Parameters const parameters_;
};

}  // namespace physics
}  // namespace principia

#include "physics/rotating_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
