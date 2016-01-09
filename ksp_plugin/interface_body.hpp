#pragma once

#include "ksp_plugin/interface.hpp"

namespace principia {
namespace interface {

inline bool operator==(Burn const& left, Burn const& right) {
  return left.thrust == right.thrust &&
         left.specific_impulse == right.specific_impulse &&
         left.frame == right.frame &&
         left.initial_time == right.initial_time &&
         left.delta_v == right.delta_v;
}

inline bool operator==(NavigationFrameParameters const& left,
                       NavigationFrameParameters const& right) {
  return left.extension == right.extension &&
         left.centre_index == right.centre_index &&
         left.primary_index == right.primary_index &&
         left.secondary_index == right.secondary_index;
}

inline bool operator==(NavigationManoeuvre const& left,
                       NavigationManoeuvre const& right) {
  return left.burn == right.burn &&
         left.initial_mass == right.initial_mass &&
         left.final_mass == right.final_mass &&
         left.mass_flow == right.mass_flow &&
         left.duration == right.duration &&
         left.final_time == right.final_time &&
         left.time_of_half_delta_v == right.time_of_half_delta_v &&
         left.time_to_half_delta_v == right.time_to_half_delta_v &&
         left.direction == right.direction;
}

inline bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

inline bool operator==(WXYZ const& left, WXYZ const& right) {
  return left.w == right.w && left.x == right.x &&
         left.y == right.y && left.z == right.z;
}

inline bool operator==(XYZ const& left, XYZ const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

inline bool operator==(XYZSegment const& left, XYZSegment const& right) {
  return left.begin == right.begin && left.end == right.end;
}

}  // namespace interface
}  // namespace principia
