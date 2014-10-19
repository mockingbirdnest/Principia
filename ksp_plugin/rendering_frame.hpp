#pragma once

#include "geometry/named_quantities.hpp"
#include "ksp_plugin/celestial.hpp"

using principia::geometry::AngularVelocity;

namespace principia {
namespace ksp_plugin {

struct Barycentre;

class RenderingFrame {
  virtual Position<Barycentre>* const RenderedPosition(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const = 0;
  virtual Velocity<Barycentre>* const RenderedVelocity(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const = 0;
};

class BodyCentredNonRotating : RenderingFrame {
 public:
  explicit BodyCentredNonRotating(Celestial<Barycentre> const& body);

  Position<Barycentre>* const RenderedPosition(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
  Velocity<Barycentre>* const RenderedVelocity(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
 private:
  Celestial<Barycentre> const& body_;
};

class BodyCentredRotatingWithSurface : RenderingFrame {
 public:
  BodyCentredRotatingWithSurface(Celestial<Barycentre> const& body,
                                 AngularVelocity<Barycentre> const& rotation_);

  Position<Barycentre>* const RenderedPosition(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
  Velocity<Barycentre>* const RenderedVelocity(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
 private:
  Celestial<Barycentre> const& body_;
  AngularVelocity<Barycentre> const rotation_;
};

class BarycentricRotating : RenderingFrame {
 public:
  BarycentricRotating(Celestial<Barycentre> const& primary,
                      Celestial<Barycentre> const& secondary_);

  Position<Barycentre>* const RenderedPosition(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
  Velocity<Barycentre>* const RenderedVelocity(
      Trajectory<Barycentre> const& trajectory,
      Instant const& t) const override;
 private:
  Celestial<Barycentre> const& primary_;
  Celestial<Barycentre> const& secondary_;
};

}  // namespace ksp_plugin
}  // namespace principia
