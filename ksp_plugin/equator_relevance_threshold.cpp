#include "ksp_plugin/equator_relevance_threshold.hpp"

#include <algorithm>

#include "numerics/elementary_functions.hpp"
#include "physics/geopotential.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _equator_relevance_threshold {
namespace internal {

using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_geopotential;
using namespace principia::physics::_oblate_body;
using namespace principia::quantities::_si;

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

}  // namespace internal
}  // namespace _equator_relevance_threshold
}  // namespace ksp_plugin
}  // namespace principia
