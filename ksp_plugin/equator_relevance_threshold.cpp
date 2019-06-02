#include "ksp_plugin/equator_relevance_threshold.hpp"

#include <algorithm>

#include "physics/geopotential.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_equator_relevance_threshold {

using physics::Geopotential;
using physics::OblateBody;
using quantities::Cbrt;
using quantities::Pow;
using quantities::si::Metre;
using quantities::si::Radian;

Length EquatorRelevanceThreshold(RotatingBody<Barycentric> const& body) {
  OblateBody<Barycentric> const* oblate_body =
      dynamic_cast<OblateBody<Barycentric> const*>(&body);
  Length const j2_threshold =
      oblate_body == nullptr
          ? 0 * Metre
          : Geopotential<Barycentric>(
                oblate_body,
                /*tolerance=*/0x1p-24).degree_damping()[2].inner_threshold();
  Length const supersynchronous_threshold =
      Cbrt(body.gravitational_parameter() /
           Pow<2>(body.angular_frequency() / (2 * Radian)));
  return std::max(j2_threshold, supersynchronous_threshold);
}

}  // namespace internal_equator_relevance_threshold
}  // namespace ksp_plugin
}  // namespace principia
