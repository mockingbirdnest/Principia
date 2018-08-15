
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

// The International Celestial Reference System.
// The origin is the barycentre of the solar system.  The axes are fixed with
// respect to the quasars.
// The xy plane is close to the mean equator at J2000.0.
// The x axis is close to the dynamical equinox at J2000.
// The z axis is perpendicular to the xy plane and points to the northern side
// of the xy plane.
using ICRS = Frame<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ICRS,
                   /*frame_is_inertial=*/true>;

// The Geocentric Celestial Reference System.
// The origin is the centre of mass of the whole Earth system, including oceans
// and atmosphere; the axes are those of the ICRS.
// To be honest:
// - The BCRS and GCRS are relativistic reference systems.
// - The transformation between BCRS and GCRS spatial coordinates has no
//   rotational component, i.e., the GCRS is kinematically non-rotating with
//   respect to the BCRS; as a result, it rotates dynamically with respect to
//   the BCRS; in particular, the equations of motion in the GCRS include de
//   Sitter and Lense–Thirring precession.
// - The time coordinate of the BCRS is TCB, and the time coordinate of the GCRS
//   is TCG; TT is a linear scaling of TCG, and TDB is a linear scaling of TDB;
//   for practical purposes TT and TDB are within 2 ms of each other;
//   Principia's |Instant| is TT.
using GCRS = Frame<serialization::Frame::SolarSystemTag,
                   serialization::Frame::GCRS,
                   /*frame_is_inertial=*/false>;

// The International Terrestrial Reference System.
// The origin is the centre of mass of the whole Earth system, including oceans
// and atmosphere; there is no global residual rotation with respect to
// horizontal motions at the earth's surface.
// Practical notes:
// - This is identical to WGS84 at 1 m level, and to recent realizations of WGS
//   84 at 10 cm level; the xz plane is a bit over a hundred metres from the
//   Royal Observatory in Greenwich.
//   Caveat: using |SphericalCoordinates<Length>| on ITRS coordinates is
//   inconsistent with WGS 84 latitudes; the WGS 84 geoid should be used.
// - The ITRS is related to the European Terrestrial Reference System (ETRS89)
//   by a time-dependent rigid transformation.
// To be honest:
// - There are many realizations of the ITRS; as of this writing, ICRFs 89, 90,
//   91, 92, 93, 94, 96, 97, 2000, 2005, 2008, and 2014.  They are related by
//   linearly time-dependent linearized similarity transformations with nonzero
//   scaling, which the current abstractions of Principia cannot support; see
//   Transformation parameters from ITRF2014 to past ITRFs,
//   http://itrf.ensg.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt.  Identifying
//   recent ITRFs (2005 and later) results in errors of a few millimetres at the
//   geocentre, scale errors on the order of 1 part-per-billion, and angular
//   errors below 1 μas.  Identifying earlier ITRFs results in errors may result
//   in errors of tens of centimetres at the geocentre, scale errors on the
//   order of ten parts per billion, and arcsecond-level angular errors.
// - ETRFyy is always related to the corresponding ITRFyy by a linearized rigid
//   transformation depending linearly on time, however it is related to other
//   ITRFs by time-dependent similarity with time-dependent nonzero scaling; see
//   the tables in EUREF Technical Note 1,
//   http://etrs89.ensg.ign.fr/pub/EUREF-TN-1.pdf
using ITRS = Frame<serialization::Frame::SolarSystemTag,
                   serialization::Frame::ITRS,
                   /*frame_is_inertial=*/false>;

// The following two reference systems are both geocentric, and they share the
// same z axis, the Celestial Intermediate Pole.  Their common xy plane is the
// intermediate equator.
// The rotation around the z axis relating the Terrestrial and Celestial
// Intermediate Reference Systems, i.e., the angle between the Terrestrial and
// Celestial Intermediate Origins, is the Earth Rotation Angle (ERA); the Earth
// Rotation Angle is an affine function of UT1 defined by equation 5.14 of the
// IERS conventions (2010).  Note that the rate is faster than 2π radians per
// UT1 day, as UT1 is notionally mean solar days, whereas revolutions of the ERA
// are stellar days.

// The glossary of the Nomenclature for Fundamental Astronomy notes for both
// that "since the acronym for this system is close to another acronym (namely
// [ITRS, respectively ICRS]), it is suggested that wherever possible the
// complete name be used".

// The Terrestrial Intermediate Reference System.
// The x axis is the Terrestrial Intermediate Origin.
// The ITRS is related to the Terrestrial Intermediate Reference System by polar
// motion, as described in sections 5.4.1 and 5.5.1 of the IERS conventions
// (2010).
// TODO(egg): identifying this with the ITRS, besides leading to acronym
// confusion, results in arcsecond-level errors.
using TerrestrialIntermediateReferenceSystem = ITRS;

// The origin is geocentric; the basis is right-handed and orthonormal, the xy
// plane is the same as that of the Terrestrial Intermediate Reference System,
// i.e., the intermediate equator.
// The x axis is the Celestial Intermediate Origin.
// The GCRS is related to the Celestial Intermediate Reference System by
// precession-nutation and frame biases, see sections 5.4.4 and 5.5.4 of the
// IERS conventions (2010).
// TODO(egg): identifying this with the GCRS leads to errors on the order of 10
// arcminutes at the time of this writing, and on the order of a degree over a
// century.
using CelestialIntermediateReferenceSystem = GCRS;

// A reference frame for transit exoplanetology.
// The xy plane is the plane of the sky.  The origin is at the barycentre of the
// extrasolar system.  The z axis goes from the extrasolar system to the Earth.
// The xz plane is (approximately) the plane of the extrasolar system.
using Sky =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::SKY, true>;

}  // namespace internal_frames

using internal_frames::ICRS;
using internal_frames::GCRS;
using internal_frames::ITRS;
using internal_frames::Sky;

}  // namespace astronomy
}  // namespace principia
