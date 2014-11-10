#pragma once

#include "physics/transforms.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"

using principia::geometry::Bivector;
using principia::geometry::Displacement;

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

// Returns the rotation matrix that maps the standard basis to the barycentric
// frame, as well as the state of the barycentre itself.
template<typename Frame>
std::pair<Matrix, DegreesOfFreedom<Frame>>
FromStandardBasisToBarycentricFrame(
    DegreesOfFreedom<Frame> const& primary_degrees_of_freedom,
    GravitationalParameter const& primary_gravitational_parameter,
    DegreesOfFreedom<Frame> const& secondary_degrees_of_freedom,
    GravitationalParameter const& secondary_gravitational_parameter) {
  Position<Frame> const barycentre_position =
      geometry::Barycentre<Displacement<Frame>, GravitationalParameter>(
          {primary_degrees_of_freedom.position,
           secondary_degrees_of_freedom.position},
          {primary_gravitational_parameter,
           secondary_gravitational_parameter});
  // TODO(phl): a barycentre on vectors (and on degrees of freedom).
  Velocity<Frame> const barycentre_velocity =
      (primary_gravitational_parameter *
           primary_degrees_of_freedom.velocity +
           secondary_gravitational_parameter *
           secondary_degrees_of_freedom.velocity) /
           (primary_gravitational_parameter +
            secondary_gravitational_parameter);
  Displacement<Frame> const reference_direction =
      primary_degrees_of_freedom.position - barycentre_position;
  Vector<double, Frame> const normalized_reference_direction =
      reference_direction / reference_direction.Norm();
  Velocity<Frame> const reference_coplanar =
      primary_degrees_of_freedom.velocity - barycentre_velocity;
  Vector<double, Frame> const normalized_reference_coplanar =
      reference_coplanar / reference_coplanar.Norm();
  // Modified Gram-Schmidt.
  Vector<double, Frame> const reference_normal =
      normalized_reference_coplanar -
          InnerProduct(normalized_reference_coplanar,
                       normalized_reference_direction) *
              normalized_reference_direction;
  // TODO(egg): should we normalize this?
  Bivector<double, Frame> const reference_binormal =
      Wedge(normalized_reference_direction, reference_normal);
  return {FromColumns(normalized_reference_direction.coordinates(),
                      reference_normal.coordinates(),
                      reference_binormal.coordinates()),
          {barycentre_position, barycentre_velocity}};
}

}  // namespace

template<typename FromFrame, typename ToFrame>
typename Trajectory<FromFrame>::TransformingIterator<ToFrame>
BodyCentredNonRotatingTransformingIterator(
    Trajectory<FromFrame> const& centre_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory) {
  CHECK_NOTNULL(transformed_trajectory);
  Trajectory<FromFrame>::Transform<ToFrame> transform =
      [&centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
    DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<FromFrame> const& last_centre_degrees_of_freedom =
        centre_trajectory.last().degrees_of_freedom();
    // on_or_after() is Ln(N), but it doesn't matter unless the map gets very
    // big, in which case we'll have cache misses anyway.
    Trajectory<FromFrame>::NativeIterator const centre_it =
        centre_trajectory.on_or_after(t);
    CHECK_EQ(centre_it.time(), t)
        << "Time " << t << " not in centre trajectory";
    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        centre_it.degrees_of_freedom();
    return {from_degrees_of_freedom.position -
                centre_degrees_of_freedom.position +
                last_centre_degrees_of_freedom.position,
            from_degrees_of_freedom.velocity -
                centre_degrees_of_freedom.velocity};
  };
  return transformed_trajectory->first_with_transform(transform);
}

template<typename FromFrame, typename ToFrame>
typename Trajectory<FromFrame>::TransformingIterator<ToFrame>
BarycentricRotatingTransformingIterator(
    Trajectory<FromFrame> const& primary_trajectory,
    Trajectory<FromFrame> const& secondary_trajectory,
    Trajectory<FromFrame> const* transformed_trajectory) {
  CHECK_NOTNULL(transformed_trajectory);

  // Start by computing the matrix that transforms from the standard basis to
  // the last barycentric frame.  We pass it by copy to the lambda so that it
  // doesn't recompute it each time.
  DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
      primary_trajectory.last().degrees_of_freedom();
  DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
      secondary_trajectory.last().degrees_of_freedom();

  auto const last_matrix_and_barycentre =
      FromStandardBasisToBarycentricFrame(
          last_primary_degrees_of_freedom,
          primary_trajectory.body().gravitational_parameter(),
          last_secondary_degrees_of_freedom,
          secondary_trajectory.body().gravitational_parameter());
  Matrix const& from_standard_basis_to_last_barycentric_frame =
      last_matrix_and_barycentre.first;
  DegreesOfFreedom<ToFrame> const& last_barycentre =
      last_matrix_and_barycentre.second;

  Trajectory<FromFrame>::Transform<ToFrame> transform =
      [&primary_trajectory,
       &secondary_trajectory,
       from_standard_basis_to_last_barycentric_frame,
       last_barycentre](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
    DegreesOfFreedom<ToFrame> {
    // on_or_after() is Ln(N).
    Trajectory<FromFrame>::NativeIterator const primary_it =
        primary_trajectory.on_or_after(t);
    CHECK_EQ(primary_it.time(), t)
        << "Time " << t << " not in primary trajectory";
    Trajectory<FromFrame>::NativeIterator secondary_it =
        secondary_trajectory.on_or_after(t);
    CHECK_EQ(secondary_it.time(), t)
        << "Time " << t << " not in secondary trajectory";

    // OUARCH!
    DegreesOfFreedom<ToFrame> const& primary_degrees_of_freedom =
        primary_it.degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const& secondary_degrees_of_freedom =
        secondary_it.degrees_of_freedom();
    auto const matrix_and_barycentre =
        FromStandardBasisToBarycentricFrame(
            primary_degrees_of_freedom,
            primary_trajectory.body().gravitational_parameter(),
            secondary_degrees_of_freedom,
            secondary_trajectory.body().gravitational_parameter());
    Matrix const& from_barycentric_frame_to_standard_basis =
        Transpose(matrix_and_barycentre.first);
    DegreesOfFreedom<ToFrame> const barycentre =
        matrix_and_barycentre.second;
    return {Displacement<ToFrame>(
                from_standard_basis_to_last_barycentric_frame(
                    from_barycentric_frame_to_standard_basis(
                        (from_degrees_of_freedom.position -
                            barycentre.position).coordinates()))) +
                last_barycentre.position,
            Velocity<ToFrame>(
                from_standard_basis_to_last_barycentric_frame(
                    from_barycentric_frame_to_standard_basis(
                        (from_degrees_of_freedom.velocity -
                            barycentre.velocity).coordinates())))};
  };
  return transformed_trajectory->first_with_transform(transform);
}

}  // namespace physics
}  // namespace principia
