
#pragma once

#include <functional>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace astronomy {
namespace internal_frames {

using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Degree;

// A reference frame with a basis.
// The frame is the International Celestial Reference Frame.
// The basis is defined from the orbit of the Earth at J2000.0 as follows:
// The xy plane is the plane of the Earth's orbit at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at J2000.0.
// The basis is right-handed and orthonormal.
using ICRFJ2000Ecliptic =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_ECLIPTIC, true>;

// The xy plane is the plane of the Earth's mean equator at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is along the Earth's mean north pole at J2000.0.
// The basis is right-handed and orthonormal.
// Note that |ICRFJ2000Equator| and |ICRFJ2000Ecliptic| share their x axis.
using ICRFJ2000Equator =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_EQUATOR, true>;

//TODO(phl):comment
using Trappist =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::TRAPPIST, true>;

// Rotation around the common x axis mapping equatorial coordinates to ecliptic
// coordinates.  The angle is the one defined by the XVIth General Assembly of
// the International Astronomical Union.
geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic> const
    ICRFJ200EquatorialToEcliptic =
        geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic>(
            23 * Degree + 26 * ArcMinute + 21.448 * ArcSecond,
            geometry::Bivector<double, ICRFJ2000Equator>({1, 0, 0}),
            geometry::DefinesFrame<ICRFJ2000Ecliptic>{});

geometry::Position<ICRFJ2000Ecliptic> const SolarSystemBarycentreEcliptic;
geometry::Position<ICRFJ2000Equator> const SolarSystemBarycentreEquator;

}  // namespace internal_frames

using internal_frames::ICRFJ2000Ecliptic;
using internal_frames::ICRFJ2000Equator;
using internal_frames::ICRFJ200EquatorialToEcliptic;
using internal_frames::SolarSystemBarycentreEcliptic;
using internal_frames::SolarSystemBarycentreEquator;

}  // namespace astronomy
}  // namespace principia
