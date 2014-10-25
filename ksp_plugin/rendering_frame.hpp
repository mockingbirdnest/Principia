#pragma once

#include "geometry/named_quantities.hpp"
#include "ksp_plugin/celestial.hpp"

using principia::geometry::AngularVelocity;

// TODO(egg): the contents of this file make little sense, and the current
// abstractions get in the way of optimizatitons.  |RenderingFrame| should
// return a |Trajectory| in the rendering frame, rather than in |Barycentre|,
// so that this trajectory may be prolonged as the history is computed, rather
// than computed every time. Since the actual rendering frame is not known at
// compile time, a wrapper (|ApparentTrajectory|?) for the trajectory *and*
// the |RenderingFrame| is needed to operate on the trajectory correctly.
// While the realization of the apparent trajectory in WorldSpace will have to
// be done at every frame, this means we will be able to do the conversion
// Barycentre -> rendering frame incrementally (except when switching reference
// frames).

namespace principia {
namespace ksp_plugin {

struct Barycentre;

class RenderingFrame {
 public:
  virtual std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const = 0;
};

class BodyCentredNonRotatingFrame : RenderingFrame {
 public:
  explicit BodyCentredNonRotatingFrame(Celestial<Barycentre> const& body);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& body_;
};

class BodyCentredRotatingWithSurface : RenderingFrame {
 public:
  BodyCentredRotatingWithSurface(Celestial<Barycentre> const& body,
                                 AngularVelocity<Barycentre> const& rotation);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& body_;
  AngularVelocity<Barycentre> const rotation_;
};

class BarycentricRotating : RenderingFrame {
 public:
  BarycentricRotating(Celestial<Barycentre> const& primary,
                      Celestial<Barycentre> const& secondary);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& primary_;
  Celestial<Barycentre> const& secondary_;
};

}  // namespace ksp_plugin
}  // namespace principia
