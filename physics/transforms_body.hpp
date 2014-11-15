#pragma once

#include "physics/transforms.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "glog/logging.h"

using principia::geometry::AffineMap;
using principia::geometry::Bivector;
using principia::geometry::Displacement;
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
Transforms<FromFrame, ThroughFrame, ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::BodyCentredNonRotating(
    Trajectory<FromFrame> const& centre_trajectory) {
  Transforms transforms;

  transforms.first_transform_ =
      [&centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
      DegreesOfFreedom<ThroughFrame> {
    // on_or_after() is Ln(N), but it doesn't matter unless the map gets very
    // big, in which case we'll have cache misses anyway.
    Trajectory<FromFrame>::NativeIterator const centre_it =
        centre_trajectory.on_or_after(t);
    CHECK_EQ(centre_it.time(), t)
        << "Time " << t << " not in centre trajectory";
    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        centre_it.degrees_of_freedom();

    AffineMap<FromFrame, ThroughFrame, Length, Permutation> position_map(
        centre_degrees_of_freedom.position,
        Position<ThroughFrame>(),
        Permutation<FromFrame, ThroughFrame>::Identity());
    return {position_map(from_degrees_of_freedom.position),
            from_degrees_of_freedom.velocity -
                centre_degrees_of_freedom.velocity};
  };

  transforms.second_transform_ =
      [&centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom) ->
      DegreesOfFreedom<ToFrame> {
    // NOTE(phl): This seems to assume that ToFrame is identical to FromFrame.
    DegreesOfFreedom<ToFrame> const& last_centre_degrees_of_freedom =
        centre_trajectory.last().degrees_of_freedom();

    AffineMap<ThroughFrame, ToFrame, Length, Permutation> position_map(
        Position<ThroughFrame>(),
        last_centre_degrees_of_freedom.position,
        Permutation<ThroughFrame, ToFrame>::Identity());
    return {position_map(through_degrees_of_freedom.position),
            through_degrees_of_freedom.velocity};
  };

  return transforms;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::BarycentricRotating(
      Trajectory<FromFrame> const& primary_trajectory,
      Trajectory<FromFrame> const& secondary_trajectory) {
  Transforms transforms;

  transforms.first_transform_ =
      [&primary_trajectory,
       &secondary_trajectory](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
      DegreesOfFreedom<ThroughFrame> {
    // on_or_after() is Ln(N).
    Trajectory<FromFrame>::NativeIterator const primary_it =
        primary_trajectory.on_or_after(t);
    CHECK_EQ(primary_it.time(), t)
        << "Time " << t << " not in primary trajectory";
    Trajectory<FromFrame>::NativeIterator secondary_it =
        secondary_trajectory.on_or_after(t);
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
            {primary_trajectory.body().gravitational_parameter(),
             secondary_trajectory.body().gravitational_parameter()});
    Matrix const from_basis_of_barycentric_frame_to_standard_basis =
        Transpose(FromStandardBasisToBasisOfBarycentricFrame(
                      barycentre_degrees_of_freedom,
                      primary_degrees_of_freedom,
                      secondary_degrees_of_freedom));
    // TODO(phl): There should be an affine map here too, once we have properly
    // 'framed' the matrix.
    return {Displacement<ThroughFrame>(
                from_basis_of_barycentric_frame_to_standard_basis(
                    (from_degrees_of_freedom.position -
                        barycentre_degrees_of_freedom.position).
                            coordinates())) + Position<ThroughFrame>(),
            Velocity<ThroughFrame>(
                from_basis_of_barycentric_frame_to_standard_basis(
                    (from_degrees_of_freedom.velocity -
                        barycentre_degrees_of_freedom.velocity).
                            coordinates()))};
  };

  transforms.second_transform_ =
      [&primary_trajectory,
       &secondary_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom) ->
      DegreesOfFreedom<ToFrame> {
    // NOTE(phl): This seems to assume that ToFrame is identical to FromFrame.
    DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
        primary_trajectory.last().degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
        secondary_trajectory.last().degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const last_barycentre =
        Barycentre<ToFrame, GravitationalParameter>(
            {last_primary_degrees_of_freedom,
             last_secondary_degrees_of_freedom},
            {primary_trajectory.body().gravitational_parameter(),
             secondary_trajectory.body().gravitational_parameter()});
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
                     Position<ThroughFrame>()).coordinates())) +
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
  return CHECK_NOTNULL(from_trajectory)->first_with_transform(first_transform_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<ThroughFrame>::template TransformingIterator<ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::second(
    Trajectory<ThroughFrame> const* through_trajectory) {
  return CHECK_NOTNULL(through_trajectory)->
             first_with_transform(second_transform_);
}

}  // namespace physics
}  // namespace principia
