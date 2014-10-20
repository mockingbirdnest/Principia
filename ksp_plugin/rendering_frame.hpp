#pragma once

#include "geometry/named_quantities.hpp"
#include "ksp_plugin/celestial.hpp"

using principia::geometry::AngularVelocity;

namespace principia {
namespace ksp_plugin {

struct Barycentre;

class RenderingFrame {
 public:
  virtual Trajectory<Barycentre> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const = 0;
};

class BodyCentredNonRotating : RenderingFrame {
 public:
  explicit BodyCentredNonRotating(Celestial<Barycentre> const& body);

  Trajectory<Barycentre> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& body_;
};

class BodyCentredRotatingWithSurface : RenderingFrame {
 public:
  BodyCentredRotatingWithSurface(Celestial<Barycentre> const& body,
                                 AngularVelocity<Barycentre> const& rotation_);

  Trajectory<Barycentre> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& body_;
  AngularVelocity<Barycentre> const rotation_;
};

class BarycentricRotating : RenderingFrame {
 public:
  BarycentricRotating(Celestial<Barycentre> const& primary,
                      Celestial<Barycentre> const& secondary_);

  Trajectory<Barycentre> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& primary_;
  Celestial<Barycentre> const& secondary_;
};

}  // namespace ksp_plugin
}  // namespace principia
