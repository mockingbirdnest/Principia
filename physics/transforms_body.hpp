#pragma once

#include "physics/transforms.hpp"

#include "glog/logging.h"

namespace principia {
namespace physics {

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
  Trajectory<FromFrame>::Transform<ToFrame> transform =
      [&centre_trajectory](
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom) ->
    DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
        primary_trajectory.last().degrees_of_freedom();
    DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
        secondary_trajectory.last().degrees_of_freedom();

    // on_or_after() is Ln(N).
    Trajectory<FromFrame>::NativeIterator const primary_it =
        primary_trajectory.on_or_after(t);
    CHECK_EQ(primary_it.time(), t)
        << "Time " << t << " not in primary trajectory";
    DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
    Trajectory<FromFrame>::NativeIterator secondary_it =
        secondary_trajectory.on_or_after(t);
    CHECK_EQ(secondary_it.time(), t)
        << "Time " << t << " not in secondary trajectory";

    // OUARCH!
    auto const last_matrix_and_barycentre =
        FromStandardBasisToBarycentricFrame(
            last_primary_degrees_of_freedom,
            primary_trajectory.body().gravitational_parameter(),
            last_secondary_degrees_of_freedom,
            secondary_trajectory.body().gravitational_parameter());
    Matrix const& from_standard_basis_to_last_barycentric_frame =
        last_matrix_and_barycentre.first;
    DegreesOfFreedom<ToFrame> const& current_barycentre =
        last_matrix_and_barycentre.second;

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
    return {Displacement<Barycentre>(
                from_standard_basis_to_last_barycentric_frame(
                    from_barycentric_frame_to_standard_basis(
                        (from_degrees_of_freedom.position -
                            barycentre.position).coordinates()))) +
                current_barycentre.position,
            Velocity<Barycentre>(
                from_standard_basis_to_last_barycentric_frame(
                    from_barycentric_frame_to_standard_basis(
                        (from_degrees_of_freedom.velocity -
                            barycentre.velocity).coordinates())))});
  };
  return transformed_trajectory->first_with_transform(transform);
}

}  // namespace physics
}  // namespace principia
