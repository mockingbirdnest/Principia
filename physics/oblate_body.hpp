
// The files containing the tree of of child classes of |Body| must be included
// in the order of inheritance to avoid circular dependencies.  This class will
// end up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#include "physics/massive_body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
#define PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_

#include <vector>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using geometry::Vector;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;
using quantities::Order2ZonalCoefficient;

namespace physics {

template<typename Frame>
class OblateBody : public MassiveBody {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  OblateBody(GravitationalParameter const& gravitational_parameter,
             double const j2,
             Length const& radius,
             Vector<double, Frame> const& axis);
  OblateBody(Mass const& mass,
             double const j2,
             Length const& radius,
             Vector<double, Frame> const& axis);
  OblateBody(GravitationalParameter const& gravitational_parameter,
             Order2ZonalCoefficient const& j2,
             Vector<double, Frame> const& axis);
  OblateBody(Mass const& mass,
             Order2ZonalCoefficient const& j2,
             Vector<double, Frame> const& axis);
  ~OblateBody() = default;

  // Returns the j2 coefficient.
  Order2ZonalCoefficient const& j2() const;

  // Returns the axis passed at construction.
  Vector<double, Frame> const& axis() const;

  // Returns false.
  bool is_massless() const;

  // Returns true.
  bool is_oblate() const;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const override;

  // Fails unless |message.has_massive_body()|.
  static not_null<std::unique_ptr<OblateBody<Frame>>> ReadFromMessage(
      serialization::Body const& message);

  // Fails if the |OblateBody| extension is absent from the message.
  static not_null<std::unique_ptr<OblateBody<Frame>>> ReadFromMessage(
      serialization::MassiveBody const& message);

 private:
  Order2ZonalCoefficient const j2_;
  Vector<double, Frame> const axis_;
};

}  // namespace physics
}  // namespace principia

#include "physics/oblate_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
