
#include "ksp_plugin/interface.hpp"

#include "base/array.hpp"
#include "base/status.hpp"
#include "base/status_or.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
namespace interface {

using base::Error;
using base::StatusOr;
using base::UniqueBytes;

namespace {

int error_code(Error error) {
  return static_cast<int>(error);
}

char const* release_string(const std::string& s) {
  UniqueBytes allocated_string(s.size() + 1);
  std::memcpy(allocated_string.data.get(), s.data(), s.size() + 1);
  return reinterpret_cast<char const*>(allocated_string.data.release());
}

}  // namespace

int principia__ExternalFlowBodyCentric(
    Plugin const* const plugin,
    int const central_body_index,
    QP const world_body_centric_initial_degrees_of_freedom,
    double const t_initial,
    double const t_final,
    QP* const world_body_centric_final_degrees_of_freedom,
    char const** const error_message) {
  journal::Method<journal::ExternalFlowBodyCentric> m{
      {plugin,
       central_body_index,
       world_body_centric_initial_degrees_of_freedom,
       t_initial,
       t_final},
      {world_body_centric_final_degrees_of_freedom, error_message}};
  if (plugin == nullptr) {
    *error_message = release_string("|plugin| must not be null");
    return m.Return(error_code(Error::INVALID_ARGUMENT));
  }
  return error_code(Error::UNIMPLEMENTED);
}

int principia__ExternalGetNearestPlannedCoastDegreesOfFreedom(
    Plugin const* const plugin,
    int const central_body_index,
    char const* const vessel_guid,
    int manoeuvre_index,
    XYZ world_body_centric_reference_position,
    QP* world_body_centric_nearest_degrees_of_freedom,
    char const** const error_message) {
  journal::Method<journal::ExternalGetNearestPlannedCoastDegreesOfFreedom> m{
      {plugin,
       central_body_index,
       vessel_guid,
       manoeuvre_index,
       world_body_centric_reference_position},
      {world_body_centric_nearest_degrees_of_freedom, error_message}};
  return m.Return(error_code(Error::UNIMPLEMENTED));
}

}  // namespace interface
}  // namespace principia
