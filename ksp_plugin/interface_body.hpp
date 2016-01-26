
#pragma once

#include "ksp_plugin/interface.hpp"

namespace principia {
namespace interface {

inline bool operator==(Burn const& left, Burn const& right) {
  return left.thrust_in_kilonewtons == right.thrust_in_kilonewtons &&
         left.specific_impulse_in_seconds_g0 ==
             right.specific_impulse_in_seconds_g0 &&
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
         left.initial_mass_in_tonnes == right.initial_mass_in_tonnes &&
         left.final_mass_in_tonnes == right.final_mass_in_tonnes &&
         left.mass_flow == right.mass_flow &&
         left.duration == right.duration &&
         left.final_time == right.final_time &&
         left.time_of_half_delta_v == right.time_of_half_delta_v &&
         left.time_to_half_delta_v == right.time_to_half_delta_v &&
         left.inertial_direction == right.inertial_direction;
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

inline R3Element<double> ToR3Element(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

inline WXYZ ToWXYZ(Quaternion const& quaternion) {
  return {quaternion.real_part(),
          quaternion.imaginary_part().x,
          quaternion.imaginary_part().y,
          quaternion.imaginary_part().z};
}

inline XYZ ToXYZ(R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

}  // namespace interface
}  // namespace principia
