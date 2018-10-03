
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#include "physics/massive_body.hpp"
#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#ifndef PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#define PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/rotation.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_rotating_body {

using base::not_null;
using geometry::AngularVelocity;
using geometry::DefinesFrame;
using geometry::EulerAngles;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::si::Radian;

template<typename Frame>
class RotatingBody : public MassiveBody {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  class PHYSICS_DLL Parameters final {
   public:
    // |reference_angle| is the angle of the prime meridian at
    // |reference_instant|.  |angular_frequency| gives the rate of rotation of
    // the body around the pole (it may be negative, as is the convention for
    // planets and satellites whose rotation is retrograde).  The direction of
    // the pole is specified in |Frame| using |right_ascension_of_pole| and
    // |declination_of_pole|.
    Parameters(Length const& mean_radius,
               Angle const& reference_angle,
               Instant const& reference_instant,
               AngularFrequency const& angular_frequency,
               Angle const& right_ascension_of_pole,
               Angle const& declination_of_pole);

   private:
    Length const mean_radius_;
    Angle const reference_angle_;
    Instant const reference_instant_;
    AngularFrequency const angular_frequency_;
    Angle const right_ascension_of_pole_;
    Angle const declination_of_pole_;
    template<typename F>
    friend class RotatingBody;
  };

  RotatingBody(MassiveBody::Parameters const& massive_body_parameters,
               Parameters const& parameters);

  // Returns the radius passed at construction.
  Length mean_radius() const override;

  // Returns the direction defined by the right ascension and declination passed
  // at construction.
  Vector<double, Frame> const& polar_axis() const;

  // Two unit vectors in the equatorial plane of the body.  |biequatorial| is
  // also in the equatorial plane of Frame.  The basis |biequatorial|,
  // |equatorial|, |polar_axis| is direct.
  Vector<double, Frame> const& biequatorial() const;
  Vector<double, Frame> const& equatorial() const;

  // Returns the right ascension passed at construction.
  Angle const& right_ascension_of_pole() const;

  // Returns the declination at construction.
  Angle const& declination_of_pole() const;

  // Returns the angular frequency passed at construction.
  AngularFrequency const& angular_frequency() const;

  // Returns the angular velocity defined by the right ascension, declination,
  // and angular frequency passed at construction.
  AngularVelocity<Frame> const& angular_velocity() const;

  // Returns the position at time |t|.
  Angle AngleAt(Instant const& t) const;

  // Returns the rotation relating the reference frame of the surface of this
  // body to |Frame|.  The reference frame of the surface is defined as follows:
  //   - the z axis is the |polar_axis|;
  //   - the x axis points from the centre of the body to the reference
  //     meridian;
  //   - the reference frame has the same handedness as |Frame|.
  // Following this definition,
  //   Displacement<SurfaceFrame>(RadiusLatitudeLongitude(r, φ, λ))
  // converts planetocentric coordinates to a displacement from the body centre
  // in the surface frame.  Note that RadiusLatitudeLongitude should *not* be
  // used for planetographic coordinates, which use a latitude defined from
  // a reference ellipsoid, and a longitude around the negative pole---except
  // for the Earth, the Moon, and the Sun.
  // In the case of the Earth, see geodetic vs. geocentric latitudes.
  template<typename SurfaceFrame>
  Rotation<SurfaceFrame, Frame> FromSurfaceFrame(Instant const& t) const;
  template<typename SurfaceFrame>
  Rotation<Frame, SurfaceFrame> ToSurfaceFrame(Instant const& t) const;

  // Returns the rotation at time |t|.
  Rotation<Frame, Frame> RotationAt(Instant const& t) const;

  // Returns false.
  bool is_massless() const override;

  // Returns false.
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
  Vector<double, Frame> const polar_axis_;
  Vector<double, Frame> biequatorial_;
  Vector<double, Frame> equatorial_;
  AngularVelocity<Frame> const angular_velocity_;
};

// Define template member functions even when importing: these are not
// dllexport.

template<typename Frame>
template<typename SurfaceFrame>
Rotation<SurfaceFrame, Frame> RotatingBody<Frame>::FromSurfaceFrame(
    Instant const& t) const {
  return Rotation<SurfaceFrame, Frame>(
      π / 2 * Radian + right_ascension_of_pole(),
      π / 2 * Radian - declination_of_pole(),
      AngleAt(t),
      EulerAngles::ZXZ,
      DefinesFrame<SurfaceFrame>{});
}

template<typename Frame>
template<typename SurfaceFrame>
Rotation<Frame, SurfaceFrame> RotatingBody<Frame>::ToSurfaceFrame(
    Instant const& t) const {
  return FromSurfaceFrame<SurfaceFrame>(t).Inverse();
}

}  // namespace internal_rotating_body

using internal_rotating_body::RotatingBody;

}  // namespace physics
}  // namespace principia

#if !PHYSICS_DLL_IMPORT
#include "physics/rotating_body_body.hpp"
#endif

#endif  // PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
