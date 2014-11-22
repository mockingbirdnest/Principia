#pragma once

#include "physics/transforms.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "glog/logging.h"
#include "physics/massive_body.hpp"

using principia::geometry::AffineMap;
using principia::geometry::Bivector;
using principia::geometry::Displacement;
using principia::geometry::Identity;
using principia::geometry::Permutation;
using principia::geometry::Position;

namespace principia {
namespace physics {

namespace {

// TODO(egg): Move this somewhere more appropriate, wrap it, etc.
struct Matrix {
  geometry::R3Element<double> row_x;
  geometry::R3Element<double> row_y;
  geometry::R3Element<double> row_z;

  template <typename Scalar>
  geometry::R3Element<Scalar> operator()(
      geometry::R3Element<Scalar> const& right) const {
    return {geometry::Dot(row_x, right),
            geometry::Dot(row_y, right),
            geometry::Dot(row_z, right)};
  }
};

Matrix FromColumns(geometry::R3Element<double> const& column_x,
                   geometry::R3Element<double> const& column_y,
                   geometry::R3Element<double> const& column_z) {
  return {{column_x.x, column_y.x, column_z.x},
          {column_x.y, column_y.y, column_z.y},
          {column_x.z, column_y.z, column_z.z}};
}

Matrix Transpose(Matrix const& m) {
  return FromColumns(m.row_x, m.row_y, m.row_z);
}

// Returns the rotation matrix that maps the standard basis to the basis of the
// barycentric frame.  |barycentre_degrees_of_freedom| must be a convex
// combination of the two other parameters.
template<typename Frame>
Matrix FromStandardBasisToBasisOfBarycentricFrame(
    DegreesOfFreedom<Frame> const& barycentre_degrees_of_freedom,
    DegreesOfFreedom<Frame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<Frame> const& secondary_degrees_of_freedom) {
  Displacement<Frame> const reference_direction =
      primary_degrees_of_freedom.position -
      barycentre_degrees_of_freedom.position;
  Vector<double, Frame> const normalized_reference_direction =
      Normalize(reference_direction);
  Velocity<Frame> const reference_coplanar =
      primary_degrees_of_freedom.velocity -
      barycentre_degrees_of_freedom.velocity;
  Vector<double, Frame> const normalized_reference_coplanar =
      Normalize(reference_coplanar);
  // Modified Gram-Schmidt.
  Vector<double, Frame> const reference_normal =
      normalized_reference_coplanar -
          InnerProduct(normalized_reference_coplanar,
                       normalized_reference_direction) *
              normalized_reference_direction;
  // TODO(egg): should we normalize this?
  Bivector<double, Frame> const reference_binormal =
      Wedge(normalized_reference_direction, reference_normal);
  return FromColumns(normalized_reference_direction.coordinates(),
                     reference_normal.coordinates(),
                     reference_binormal.coordinates());
}

}  // namespace

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
std::unique_ptr<Transforms<FromFrame, ThroughFrame, ToFrame>>
Transforms<FromFrame, ThroughFrame, ToFrame>::BodyCentredNonRotating(
    Trajectory<FromFrame> const& from_centre_trajectory,
    Trajectory<ToFrame> const& to_centre_trajectory) {
  std::unique_ptr<Transforms> transforms(new Transforms);

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  Transforms* that = transforms.get();
  transforms->first_ =
      [&from_centre_trajectory, that](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
      DegreesOfFreedom<ThroughFrame> {
    // First check if the result is cached.
    auto cache_it = that->first_cache_.find(t);
    if (cache_it != that->first_cache_.end()) {
      return cache_it->second;
    }

    // on_or_after() is Ln(N), but it doesn't matter unless the map gets very
    // big, in which case we'll have cache misses anyway.
    Trajectory<FromFrame>::NativeIterator const centre_it =
        from_centre_trajectory.on_or_after(t);
    CHECK_EQ(centre_it.time(), t)
        << "Time " << t << " not in centre trajectory";
    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        centre_it.degrees_of_freedom();

    AffineMap<FromFrame, ThroughFrame, Length, Identity> position_map(
        centre_degrees_of_freedom.position,
        ThroughFrame::origin,
        Identity<FromFrame, ThroughFrame>());
    Identity<FromFrame, ThroughFrame> velocity_map;
    DegreesOfFreedom<ThroughFrame> through_degrees_of_freedom =
        {position_map(from_degrees_of_freedom.position),
         velocity_map(from_degrees_of_freedom.velocity -
                      centre_degrees_of_freedom.velocity)};

    // Cache the result before returning it.
    that->first_cache_.emplace(t, through_degrees_of_freedom);
    return std::move(through_degrees_of_freedom);
  };

  transforms->second_ =
      [&to_centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom) ->
      DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<ToFrame> const& last_centre_degrees_of_freedom =
        to_centre_trajectory.last().degrees_of_freedom();

    AffineMap<ThroughFrame, ToFrame, Length, Identity> position_map(
        ThroughFrame::origin,
        last_centre_degrees_of_freedom.position,
        Identity<ThroughFrame, ToFrame>());
    Identity<ThroughFrame, ToFrame> velocity_map;
    return {position_map(through_degrees_of_freedom.position),
            velocity_map(through_degrees_of_freedom.velocity)};
  };

  return transforms;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
std::unique_ptr<Transforms<FromFrame, ThroughFrame, ToFrame>>
Transforms<FromFrame, ThroughFrame, ToFrame>::BarycentricRotating(
      Trajectory<FromFrame> const& from_primary_trajectory,
      Trajectory<ToFrame> const& to_primary_trajectory,
      Trajectory<FromFrame> const& from_secondary_trajectory,
      Trajectory<ToFrame> const& to_secondary_trajectory) {
  std::unique_ptr<Transforms> transforms(new Transforms);

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  Transforms* that = transforms.get();
  transforms->first_ =
      [&from_primary_trajectory, &from_secondary_trajectory, that](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
      DegreesOfFreedom<ThroughFrame> {
    // First check if the result is cached.
    auto cache_it = that->first_cache_.find(t);
    if (cache_it != that->first_cache_.end()) {
      return cache_it->second;
    }

    // on_or_after() is Ln(N).
    Trajectory<FromFrame>::NativeIterator const primary_it =
        from_primary_trajectory.on_or_after(t);
    CHECK_EQ(primary_it.time(), t)
        << "Time " << t << " not in primary trajectory";
    Trajectory<FromFrame>::NativeIterator secondary_it =
        from_secondary_trajectory.on_or_after(t);
    CHECK_EQ(secondary_it.time(), t)
        << "Time " << t << " not in secondary trajectory";

    DegreesOfFreedom<FromFrame> const& primary_degrees_of_freedom =
        primary_it.degrees_of_freedom();
    DegreesOfFreedom<FromFrame> const& secondary_degrees_of_freedom =
        secondary_it.degrees_of_freedom();
    DegreesOfFreedom<FromFrame> const barycentre_degrees_of_freedom =
        Barycentre<FromFrame, GravitationalParameter>(
            {primary_degrees_of_freedom,
             secondary_degrees_of_freedom},
            {from_primary_trajectory.body<MassiveBody>().
                 gravitational_parameter(),
             from_secondary_trajectory.body<MassiveBody>().
                 gravitational_parameter()});
    Matrix const from_basis_of_barycentric_frame_to_standard_basis =
        Transpose(FromStandardBasisToBasisOfBarycentricFrame(
                      barycentre_degrees_of_freedom,
                      primary_degrees_of_freedom,
                      secondary_degrees_of_freedom));
    // TODO(phl): There should be an affine map here too, once we have properly
    // 'framed' the matrix.
    DegreesOfFreedom<ThroughFrame> through_degrees_of_freedom =
        {Displacement<ThroughFrame>(
             from_basis_of_barycentric_frame_to_standard_basis(
                 (from_degrees_of_freedom.position -
                      barycentre_degrees_of_freedom.position).
                  coordinates())) + ThroughFrame::origin,
         Velocity<ThroughFrame>(
             from_basis_of_barycentric_frame_to_standard_basis(
                 (from_degrees_of_freedom.velocity -
                     barycentre_degrees_of_freedom.velocity).
                  coordinates()))};

    // Cache the result before returning it.
    that->first_cache_.emplace(t, through_degrees_of_freedom);
    return std::move(through_degrees_of_freedom);
  };

  transforms->second_ =
      [&to_primary_trajectory, &to_secondary_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom) ->
      DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
        to_primary_trajectory.last().degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
        to_secondary_trajectory.last().degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const last_barycentre =
        Barycentre<ToFrame, GravitationalParameter>(
            {last_primary_degrees_of_freedom,
             last_secondary_degrees_of_freedom},
            {to_primary_trajectory.body<MassiveBody>().
                 gravitational_parameter(),
             to_secondary_trajectory.body<MassiveBody>().
                 gravitational_parameter()});
    Matrix const from_standard_basis_to_basis_of_last_barycentric_frame =
        FromStandardBasisToBasisOfBarycentricFrame(
            last_barycentre,
            last_primary_degrees_of_freedom,
            last_secondary_degrees_of_freedom);
    // TODO(phl): There should be an affine map here too, once we have properly
    // 'framed' the matrix.
    return {Displacement<ToFrame>(
                from_standard_basis_to_basis_of_last_barycentric_frame(
                    (through_degrees_of_freedom.position -
                     ThroughFrame::origin).coordinates())) +
                last_barycentre.position,
            Velocity<ToFrame>(
                from_standard_basis_to_basis_of_last_barycentric_frame(
                    through_degrees_of_freedom.velocity.coordinates()))};
  };

  return transforms;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::first(
    Trajectory<FromFrame> const* from_trajectory) {
  return CHECK_NOTNULL(from_trajectory)->first_with_transform(first_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<ThroughFrame>::template TransformingIterator<ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::second(
    Trajectory<ThroughFrame> const* through_trajectory) {
  return CHECK_NOTNULL(through_trajectory)->
             first_with_transform(second_);
}

}  // namespace physics
}  // namespace principia
