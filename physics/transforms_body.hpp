#pragma once

#include "physics/transforms.hpp"

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/permutation.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::make_not_null_unique;
using geometry::AffineMap;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Identity;
using geometry::Permutation;
using geometry::Position;
using geometry::R3x3Matrix;
using geometry::Rotation;
using geometry::Wedge;
using quantities::AngularFrequency;
using quantities::Pow;
using si::Radian;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace physics {

namespace {

// Fills |*rotation| with the rotation that maps the basis of the barycentric
// frame to the standard basis.  Fills |*angular_frequency| with the
// corresponding angular velocity.  These pointers must be nonnull, and there is
// no transfer of ownership.  |barycentre_degrees_of_freedom| must be a convex
// combination of the two other degrees of freedom.
template<typename FromFrame, typename ToFrame>
void FromBasisOfBarycentricFrameToStandardBasis(
    DegreesOfFreedom<FromFrame> const& barycentre_degrees_of_freedom,
    DegreesOfFreedom<FromFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<FromFrame> const& secondary_degrees_of_freedom,
    not_null<Rotation<FromFrame, ToFrame>*> const rotation,
    not_null<Bivector<AngularFrequency, FromFrame>*> const angular_frequency) {
  RelativeDegreesOfFreedom<FromFrame> const reference =
      primary_degrees_of_freedom - barycentre_degrees_of_freedom;
  Displacement<FromFrame> const& reference_direction =
      reference.displacement();
  Velocity<FromFrame> reference_normal = reference.velocity();
  reference_direction.template Orthogonalize<Speed, FromFrame>(
      &reference_normal);
  Bivector<Product<Length, Speed>, FromFrame> const reference_binormal =
      Wedge(reference_direction, reference_normal);
  *rotation = Rotation<FromFrame, ToFrame>(
                  R3x3Matrix(Normalize(reference_direction).coordinates(),
                             Normalize(reference_normal).coordinates(),
                             Normalize(reference_binormal).coordinates()));
  *angular_frequency =
      (Radian / Pow<2>(reference_direction.Norm())) * reference_binormal;
}

template<typename ThroughFrame, typename ToFrame, typename LazyTrajectory>
void FromStandardBasisToBasisOfLastBarycentricFrame(
    LazyTrajectory const& to_primary_trajectory,
    LazyTrajectory const& to_secondary_trajectory,
    not_null<Rotation<ThroughFrame, ToFrame>*> const rotation,
    not_null<DegreesOfFreedom<ToFrame>*> const
        last_barycentre_degrees_of_freedom) {
  DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
      to_primary_trajectory().last().degrees_of_freedom();
  DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
      to_secondary_trajectory().last().degrees_of_freedom();
  *last_barycentre_degrees_of_freedom =
      Barycentre<ToFrame, GravitationalParameter>(
          {last_primary_degrees_of_freedom,
           last_secondary_degrees_of_freedom},
          {to_primary_trajectory().template body<MassiveBody>()->
               gravitational_parameter(),
           to_secondary_trajectory().template body<MassiveBody>()->
               gravitational_parameter()});
  Rotation<ToFrame, ThroughFrame>
      from_basis_of_last_barycentric_frame_to_standard_basis =
          Rotation<ToFrame, ThroughFrame>::Identity();
  Bivector<AngularFrequency, ToFrame> angular_frequency;
  FromBasisOfBarycentricFrameToStandardBasis<ToFrame, ThroughFrame>(
      *last_barycentre_degrees_of_freedom,
      last_primary_degrees_of_freedom,
      last_secondary_degrees_of_freedom,
      &from_basis_of_last_barycentric_frame_to_standard_basis,
      &angular_frequency);
  *rotation = from_basis_of_last_barycentric_frame_to_standard_basis.Inverse();
}

}  // namespace

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>>>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::BodyCentredNonRotating(
    Mobile const& centre,
    LazyTrajectory<ToFrame> const& to_trajectory) {
  not_null<std::unique_ptr<Transforms>> transforms =
      make_not_null_unique<Transforms>();

  transforms->coordinate_frame_ = CoordinateFrame<ToFrame>();

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  not_null<Transforms*> that = transforms.get();
  transforms->first_ =
      [&centre, that](
          LazyTrajectory<FromFrame> const& from_trajectory,
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom,
          not_null<Trajectory<FromFrame> const*> const trajectory) ->
      DegreesOfFreedom<ThroughFrame> {
    // First check if the result is cached.
    DegreesOfFreedom<ThroughFrame>* cached_through_degrees_of_freedom = nullptr;
    if (that->first_cache_.Lookup(trajectory, t,
                                  &cached_through_degrees_of_freedom)) {
      return *cached_through_degrees_of_freedom;
    }

    // on_or_after() is Ln(N), but it doesn't matter unless the map gets very
    // big, in which case we'll have cache misses anyway.
    TYPENAME Trajectory<FromFrame>::NativeIterator const centre_it =
        (centre.*from_trajectory)().on_or_after(t);
    CHECK_EQ(centre_it.time(), t)
        << "Time " << t << " not in centre trajectory";
    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        centre_it.degrees_of_freedom();

    AffineMap<FromFrame, ThroughFrame, Length, Identity> const position_map(
        centre_degrees_of_freedom.position(),
        ThroughFrame::origin,
        Identity<FromFrame, ThroughFrame>());
    // TODO(phl): Should |velocity_map| be an affine map?
    Identity<FromFrame, ThroughFrame> const velocity_map;
    DegreesOfFreedom<ThroughFrame> through_degrees_of_freedom =
        {position_map(from_degrees_of_freedom.position()),
         velocity_map(from_degrees_of_freedom.velocity() -
                      centre_degrees_of_freedom.velocity())};

    // Cache the result before returning it.
    that->first_cache_.Insert(trajectory, t, through_degrees_of_freedom);
    return through_degrees_of_freedom;
  };

  transforms->second_ =
      [&centre, to_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom,
          Trajectory<ThroughFrame> const* trajectory) ->
      DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<ToFrame> const& last_centre_degrees_of_freedom =
        (centre.*to_trajectory)().last().degrees_of_freedom();

    AffineMap<ThroughFrame, ToFrame, Length, Identity> const position_map(
        ThroughFrame::origin,
        last_centre_degrees_of_freedom.position(),
        Identity<ThroughFrame, ToFrame>());
    Identity<ThroughFrame, ToFrame> const velocity_map;
    return {position_map(through_degrees_of_freedom.position()),
            velocity_map(through_degrees_of_freedom.velocity())};
  };

  return transforms;
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>>>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::BarycentricRotating(
    Mobile const& primary,
    Mobile const& secondary,
    LazyTrajectory<ToFrame> const& to_trajectory) {
  not_null<std::unique_ptr<Transforms>> transforms =
      make_not_null_unique<Transforms>();

  transforms->coordinate_frame_ =
      [to_primary_trajectory, to_secondary_trajectory](
          Position<ToFrame> const& q) {
    Rotation<ThroughFrame, ToFrame>
        from_standard_basis_to_basis_of_last_barycentric_frame =
            Rotation<ThroughFrame, ToFrame>::Identity();
    DegreesOfFreedom<ToFrame> dummy = {ToFrame::origin, Velocity<ToFrame>()};
    FromStandardBasisToBasisOfLastBarycentricFrame<ThroughFrame, ToFrame>(
        to_primary_trajectory,
        to_secondary_trajectory,
        &from_standard_basis_to_basis_of_last_barycentric_frame,
        &dummy);

    return from_standard_basis_to_basis_of_last_barycentric_frame *
               Rotation<ToFrame, ThroughFrame>::Identity();
  };

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  not_null<Transforms*> that = transforms.get();
  transforms->first_ =
      [&primary, &secondary, that](
          LazyTrajectory<FromFrame> const& from_trajectory,
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom,
          not_null<Trajectory<FromFrame> const*> const trajectory) ->
      DegreesOfFreedom<ThroughFrame> {
    // First check if the result is cached.
    DegreesOfFreedom<ThroughFrame>* cached_through_degrees_of_freedom = nullptr;
    if (that->first_cache_.Lookup(trajectory, t,
                                  &cached_through_degrees_of_freedom)) {
      return *cached_through_degrees_of_freedom;
    }

    // |on_or_after()| is Ln(N).
    TYPENAME Trajectory<FromFrame>::NativeIterator const primary_it =
        (primary.*from_trajectory)().on_or_after(t);
    CHECK_EQ(primary_it.time(), t)
        << "Time " << t << " not in primary trajectory";
    TYPENAME Trajectory<FromFrame>::NativeIterator secondary_it =
        (secondary.*from_trajectory)().on_or_after(t);
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
            {(primary.*from_trajectory)().template body<MassiveBody>()->
                 gravitational_parameter(),
             (secondary.*from_trajectory)().template body<MassiveBody>()->
                 gravitational_parameter()});
    Rotation<FromFrame, ThroughFrame>
        from_basis_of_barycentric_frame_to_standard_basis =
            Rotation<FromFrame, ThroughFrame>::Identity();
    Bivector<AngularFrequency, FromFrame> angular_frequency;
    FromBasisOfBarycentricFrameToStandardBasis<FromFrame, ThroughFrame>(
        barycentre_degrees_of_freedom,
        primary_degrees_of_freedom,
        secondary_degrees_of_freedom,
        &from_basis_of_barycentric_frame_to_standard_basis,
        &angular_frequency);

    AffineMap<FromFrame, ThroughFrame, Length, Rotation> const position_map(
        barycentre_degrees_of_freedom.position(),
        ThroughFrame::origin,
        from_basis_of_barycentric_frame_to_standard_basis);
    // TODO(phl): This is where we wonder if |velocity_map| should be an affine
    // map.  Also, the filioque.
    Rotation<FromFrame, ThroughFrame> const& velocity_map =
        from_basis_of_barycentric_frame_to_standard_basis;
    DegreesOfFreedom<ThroughFrame> through_degrees_of_freedom =
        {position_map(from_degrees_of_freedom.position()),
         velocity_map(from_degrees_of_freedom.velocity() -
                      barycentre_degrees_of_freedom.velocity() -
                        angular_frequency *
                          (from_degrees_of_freedom.position() -
                           barycentre_degrees_of_freedom.position()) / Radian)};

    // Cache the result before returning it.
    that->first_cache_.Insert(trajectory, t, through_degrees_of_freedom);
    return through_degrees_of_freedom;
  };

  transforms->second_ =
      [&primary, &secondary, to_trajectory](
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom,
          Trajectory<ThroughFrame> const* trajectory) ->
      DegreesOfFreedom<ToFrame> {
    Rotation<ThroughFrame, ToFrame>
        from_standard_basis_to_basis_of_last_barycentric_frame =
            Rotation<ThroughFrame, ToFrame>::Identity();
    DegreesOfFreedom<ToFrame> last_barycentre_degrees_of_freedom =
        {ToFrame::origin, Velocity<ToFrame>()};
    FromStandardBasisToBasisOfLastBarycentricFrame<ThroughFrame, ToFrame>(
        to_primary_trajectory,
        to_secondary_trajectory,
        &from_standard_basis_to_basis_of_last_barycentric_frame,
        &last_barycentre_degrees_of_freedom);

    AffineMap<ThroughFrame, ToFrame, Length, Rotation> const position_map(
        ThroughFrame::origin,
        last_barycentre_degrees_of_freedom.position(),
        from_standard_basis_to_basis_of_last_barycentric_frame);
    Rotation<ThroughFrame, ToFrame> const& velocity_map =
        from_standard_basis_to_basis_of_last_barycentric_frame;
    return {position_map(through_degrees_of_freedom.position()),
            velocity_map(through_degrees_of_freedom.velocity())};
  };

  return transforms;
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>>>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::DummyForTesting() {
  return make_not_null_unique<Transforms>();
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::first(
    Mobile const& mobile,
    LazyTrajectory<FromFrame> const& from_trajectory) {
  typename Trajectory<FromFrame>::template Transform<ThroughFrame> const first =
      std::bind(first_, from_trajectory, _1, _2, _3);
  return (mobile.*from_trajectory)().first_with_transform(first);
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::first_on_or_after(
    Mobile const& mobile,
    LazyTrajectory<FromFrame> const& from_trajectory,
    Instant const& time) {
  typename Trajectory<FromFrame>::template Transform<ThroughFrame> const first =
      std::bind(first_, from_trajectory, _1, _2, _3);
  return (mobile.*from_trajectory)().on_or_after_with_transform(time, first);
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<ThroughFrame>::template TransformingIterator<ToFrame>
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::second(
    Trajectory<ThroughFrame> const& through_trajectory) {
  return through_trajectory.first_with_transform(second_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
FrameField<ToFrame>
Transforms<FromFrame, ThroughFrame, ToFrame>::coordinate_frame() const {
  return coordinate_frame_;
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
template<typename Frame1, typename Frame2>
bool
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::Cache<Frame1, Frame2>::
Lookup(not_null<Trajectory<Frame1> const*> const trajectory,
    Instant const& time,
    not_null<DegreesOfFreedom<Frame2>**> degrees_of_freedom) {
  bool found = false;
  ++number_of_lookups_[trajectory];
  auto const it = map_.find(std::make_pair(trajectory, time));
  if (it != map_.end()) {
    ++number_of_hits_[trajectory];
    *degrees_of_freedom = &it->second;
    found = true;
  }
  VLOG_EVERY_N(1, 1000) << "Hit ratio for trajectory " << trajectory << " is "
                        << static_cast<double>(number_of_hits_[trajectory]) /
                           number_of_lookups_[trajectory];
  return found;
}

template<typename Mobile,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
template<typename Frame1, typename Frame2>
void
Transforms<Mobile, FromFrame, ThroughFrame, ToFrame>::Cache<Frame1, Frame2>::
Insert(not_null<Trajectory<Frame1> const*> const trajectory,
    Instant const& time,
    DegreesOfFreedom<Frame2> const& degrees_of_freedom) {
  map_.emplace(std::make_pair(trajectory, time), degrees_of_freedom);
}

}  // namespace physics
}  // namespace principia
