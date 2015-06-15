#pragma once

#include <algorithm>
#include <list>

#include "physics/transformz.hpp"

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
// corresponding angular velocity.  |barycentre_degrees_of_freedom| must be a
// convex combination of the two other degrees of freedom.
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

template<typename ThroughFrame, typename ToFrame>
void FromStandardBasisToBasisOfLastBarycentricFrame(
    Instant const& last,
    MassiveBody const& primary,
    ContinuousTrajectory<ToFrame> const& to_primary_trajectory,
    typename std::list<
        typename ContinuousTrajectory<ToFrame>::Hint>::iterator const
        to_primary_hint,
    MassiveBody const& secondary,
    ContinuousTrajectory<ToFrame> const& to_secondary_trajectory,
    typename std::list<
        typename ContinuousTrajectory<ToFrame>::Hint>::iterator const
        to_secondary_hint,
    not_null<Rotation<ThroughFrame, ToFrame>*> const rotation,
    not_null<DegreesOfFreedom<ToFrame>*> const
        last_barycentre_degrees_of_freedom) {
  // TODO(phl): Add hinting.
  DegreesOfFreedom<ToFrame> const& last_primary_degrees_of_freedom =
      to_primary_trajectory.EvaluateDegreesOfFreedom(
          last, &*to_primary_hint);
  DegreesOfFreedom<ToFrame> const& last_secondary_degrees_of_freedom =
      to_secondary_trajectory.EvaluateDegreesOfFreedom(
          last, &*to_secondary_hint);
  *last_barycentre_degrees_of_freedom =
      Barycentre<ToFrame, GravitationalParameter>(
          {last_primary_degrees_of_freedom,
           last_secondary_degrees_of_freedom},
          {primary.gravitational_parameter(),
           secondary.gravitational_parameter()});
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

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transformz<FromFrame, ThroughFrame, ToFrame>>>
Transformz<FromFrame, ThroughFrame, ToFrame>::BodyCentredNonRotating(
    MassiveBody const& centre,
    ContinuousTrajectory<FromFrame> const& from_centre_trajectory,
    ContinuousTrajectory<ToFrame> const& to_centre_trajectory) {
  not_null<std::unique_ptr<Transformz>> transforms =
      make_not_null_unique<Transformz>();
  transforms->from_hints_.emplace_front();
  transforms->to_hints_.emplace_front();
  auto const from_centre_hint = transforms->from_hints_.begin();
  auto const to_centre_hint = transforms->to_hints_.begin();

  transforms->coordinate_frame_ =
      [](Instant const& last,
         Position<ToFrame> const& q) {
    return CoordinateFrame<ToFrame>()(q);
  };

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  not_null<Transformz*> that = transforms.get();
  transforms->first_ =
      [&from_centre_trajectory, from_centre_hint, that](
          bool const cacheable,
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom,
          not_null<Trajectory<FromFrame> const*> const trajectory) ->
      DegreesOfFreedom<ThroughFrame> {
    //TODO(phl): Not for first()!
    // First check if the result is cached.
    DegreesOfFreedom<ThroughFrame>* cached_through_degrees_of_freedom = nullptr;
    if (cacheable &&
        that->first_cache_.Lookup(trajectory, t,
                                  &cached_through_degrees_of_freedom)) {
      return *cached_through_degrees_of_freedom;
    }

    DegreesOfFreedom<FromFrame> const& centre_degrees_of_freedom =
        from_centre_trajectory.EvaluateDegreesOfFreedom(t, &*from_centre_hint);

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
    if (cacheable) {
      that->first_cache_.Insert(trajectory, t, through_degrees_of_freedom);
    }
    return through_degrees_of_freedom;
  };

  transforms->second_ =
      [&to_centre_trajectory, to_centre_hint](
          Instant const last,
          Instant const& t,
          DegreesOfFreedom<ThroughFrame> const& through_degrees_of_freedom,
          Trajectory<ThroughFrame> const* trajectory) ->
      DegreesOfFreedom<ToFrame> {
    DegreesOfFreedom<ToFrame> const& last_centre_degrees_of_freedom =
        to_centre_trajectory.EvaluateDegreesOfFreedom(last, &*to_centre_hint);

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

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transformz<FromFrame, ThroughFrame, ToFrame>>>
Transformz<FromFrame, ThroughFrame, ToFrame>::BarycentricRotating(
    MassiveBody const& primary,
    ContinuousTrajectory<FromFrame> const& from_primary_trajectory,
    ContinuousTrajectory<ToFrame> const& to_primary_trajectory,
    MassiveBody const& secondary,
    ContinuousTrajectory<FromFrame> const& from_secondary_trajectory,
    ContinuousTrajectory<ToFrame> const& to_secondary_trajectory) {
  not_null<std::unique_ptr<Transformz>> transforms =
      make_not_null_unique<Transformz>();
  transforms->from_hints_.emplace_front();
  transforms->to_hints_.emplace_front();
  auto const from_primary_hint = transforms->from_hints_.begin();
  auto const to_primary_hint = transforms->to_hints_.begin();
  transforms->from_hints_.emplace_front();
  transforms->to_hints_.emplace_front();
  auto const from_secondary_hint = transforms->from_hints_.begin();
  auto const to_secondary_hint = transforms->to_hints_.begin();

  transforms->coordinate_frame_ =
      [&primary, &to_primary_trajectory, to_primary_hint,
       &secondary, &to_secondary_trajectory, to_secondary_hint](
          Instant const& last,
          Position<ToFrame> const& q) ->
      Rotation<ToFrame, ToFrame> {
    Rotation<ThroughFrame, ToFrame>
        from_standard_basis_to_basis_of_last_barycentric_frame =
            Rotation<ThroughFrame, ToFrame>::Identity();
    DegreesOfFreedom<ToFrame> dummy = {ToFrame::origin, Velocity<ToFrame>()};
    FromStandardBasisToBasisOfLastBarycentricFrame<ThroughFrame, ToFrame>(
        last,
        primary,
        to_primary_trajectory,
        to_primary_hint,
        secondary,
        to_secondary_trajectory,
        to_secondary_hint,
        &from_standard_basis_to_basis_of_last_barycentric_frame,
        &dummy);

    return from_standard_basis_to_basis_of_last_barycentric_frame *
               Rotation<ToFrame, ThroughFrame>::Identity();
  };

  // From the perspective of the lambda the following variable is really |this|,
  // hence the name.
  not_null<Transformz*> that = transforms.get();
  transforms->first_ =
      [&primary, &from_primary_trajectory, from_primary_hint,
       &secondary, &from_secondary_trajectory, from_secondary_hint, that](
          bool const cacheable,
          Instant const& t,
          DegreesOfFreedom<FromFrame> const& from_degrees_of_freedom,
          not_null<Trajectory<FromFrame> const*> const trajectory) ->
      DegreesOfFreedom<ThroughFrame> {
    // First check if the result is cached.
    DegreesOfFreedom<ThroughFrame>* cached_through_degrees_of_freedom = nullptr;
    if (cacheable &&
        that->first_cache_.Lookup(trajectory, t,
                                  &cached_through_degrees_of_freedom)) {
      return *cached_through_degrees_of_freedom;
    }

    DegreesOfFreedom<FromFrame> const& primary_degrees_of_freedom =
        from_primary_trajectory.EvaluateDegreesOfFreedom(
            t, &*from_primary_hint);
    DegreesOfFreedom<FromFrame> const& secondary_degrees_of_freedom =
        from_secondary_trajectory.EvaluateDegreesOfFreedom(
            t, &*from_secondary_hint);
    DegreesOfFreedom<FromFrame> const barycentre_degrees_of_freedom =
        Barycentre<FromFrame, GravitationalParameter>(
            {primary_degrees_of_freedom,
             secondary_degrees_of_freedom},
            {primary.gravitational_parameter(),
             secondary.gravitational_parameter()});
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
    if (cacheable) {
      that->first_cache_.Insert(trajectory, t, through_degrees_of_freedom);
    }
    return through_degrees_of_freedom;
  };

  transforms->second_ =
      [&primary, &to_primary_trajectory, to_primary_hint,
       &secondary, &to_secondary_trajectory, to_secondary_hint](
          Instant const& last,
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
        last,
        primary,
        to_primary_trajectory,
        to_primary_hint,
        secondary,
        to_secondary_trajectory,
        to_secondary_hint,
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

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
not_null<std::unique_ptr<Transformz<FromFrame, ThroughFrame, ToFrame>>>
Transformz<FromFrame, ThroughFrame, ToFrame>::DummyForTesting() {
  return make_not_null_unique<Transformz>();
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
Transformz<FromFrame, ThroughFrame, ToFrame>::first(
    Trajectory<FromFrame> const& from_trajectory) {
  typename Trajectory<FromFrame>::template Transform<ThroughFrame> const first =
      std::bind(first_, false /*cacheable*/, _1, _2, _3);
  return from_trajectory.first_with_transform(first);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<FromFrame>::template TransformingIterator<ThroughFrame>
Transformz<FromFrame, ThroughFrame, ToFrame>::first_with_caching(
    not_null<Trajectory<FromFrame> const*> const from_trajectory) {
  typename Trajectory<FromFrame>::template Transform<ThroughFrame> const first =
      std::bind(first_, true /*cacheable*/, _1, _2, _3);
  return from_trajectory->first_with_transform(first);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
typename Trajectory<ThroughFrame>::template TransformingIterator<ToFrame>
Transformz<FromFrame, ThroughFrame, ToFrame>::second(
    Instant const& last,
    Trajectory<ThroughFrame> const& through_trajectory) {
  typename Trajectory<ThroughFrame>::template Transform<ToFrame> second =
      std::bind(second_, last, _1, _2, _3);
  return through_trajectory.first_with_transform(second);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
FrameField<ToFrame>
Transformz<FromFrame, ThroughFrame, ToFrame>::coordinate_frame(
    Instant const& last) const {
  return std::bind(coordinate_frame_, last, _1);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
template<typename Frame1, typename Frame2>
bool
Transformz<FromFrame, ThroughFrame, ToFrame>::Cache<Frame1, Frame2>::
Lookup(not_null<Trajectory<Frame1> const*> const trajectory,
       Instant const& time,
       not_null<DegreesOfFreedom<Frame2>**> degrees_of_freedom) {
  bool found = false;
  ++number_of_lookups_[trajectory];
  auto const it1 = cache_.find(trajectory);
  if (it1 != cache_.end()) {
    auto const it2 = it1->second.find(time);
    if (it2 != it1->second.end()) {
      ++number_of_hits_[trajectory];
      *degrees_of_freedom = &it2->second;
      found = true;
    }
  }
  VLOG_EVERY_N(1, 1000) << "Hit ratio for trajectory " << trajectory << " is "
                        << static_cast<double>(number_of_hits_[trajectory]) /
                           number_of_lookups_[trajectory];
  return found;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
template<typename Frame1, typename Frame2>
void
Transformz<FromFrame, ThroughFrame, ToFrame>::Cache<Frame1, Frame2>::Insert(
    not_null<Trajectory<Frame1> const*> const trajectory,
       Instant const& time,
       DegreesOfFreedom<Frame2> const& degrees_of_freedom) {
  cache_[trajectory].emplace(std::make_pair(time, degrees_of_freedom));
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
template<typename Frame1, typename Frame2>
void
Transformz<FromFrame, ThroughFrame, ToFrame>::Cache<Frame1, Frame2>::Delete(
    not_null<Trajectory<Frame1> const*> const trajectory) {
  cache_.erase[trajectory];
}

}  // namespace physics
}  // namespace principia
