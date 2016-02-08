
#pragma once

#include "ksp_plugin/interface.hpp"

#include <cmath>

namespace principia {
namespace interface {

inline bool NaNIndependentEq(double const left, double const right) {
  return (left == right) || (std::isnan(left) && std::isnan(right));
}

inline bool operator==(Burn const& left, Burn const& right) {
  return NaNIndependentEq(left.thrust_in_kilonewtons,
                          right.thrust_in_kilonewtons) &&
         NaNIndependentEq(left.specific_impulse_in_seconds_g0,
                          right.specific_impulse_in_seconds_g0) &&
         left.frame == right.frame &&
         NaNIndependentEq(left.initial_time, right.initial_time) &&
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
         NaNIndependentEq(left.initial_mass_in_tonnes,
                          right.initial_mass_in_tonnes) &&
         NaNIndependentEq(left.final_mass_in_tonnes,
                          right.final_mass_in_tonnes) &&
         NaNIndependentEq(left.mass_flow, right.mass_flow) &&
         NaNIndependentEq(left.duration, right.duration) &&
         NaNIndependentEq(left.final_time, right.final_time) &&
         NaNIndependentEq(left.time_of_half_delta_v,
                          right.time_of_half_delta_v) &&
         NaNIndependentEq(left.time_to_half_delta_v,
                          right.time_to_half_delta_v) &&
         left.inertial_direction == right.inertial_direction &&
         left.binormal == right.binormal &&
         left.normal == right.normal &&
         left.tangent == right.tangent;
}

inline bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

inline bool operator==(WXYZ const& left, WXYZ const& right) {
  return NaNIndependentEq(left.w, right.w) &&
         NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
}

inline bool operator==(XYZ const& left, XYZ const& right) {
  return NaNIndependentEq(left.x, right.x) &&
         NaNIndependentEq(left.y, right.y) &&
         NaNIndependentEq(left.z, right.z);
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

inline not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const* const plugin,
    NavigationFrameParameters const& parameters) {
  switch (parameters.extension) {
    case serialization::BarycentricRotatingDynamicFrame::
        kBarycentricRotatingDynamicFrameFieldNumber:
      return CHECK_NOTNULL(plugin)->NewBarycentricRotatingNavigationFrame(
          parameters.primary_index, parameters.secondary_index);
    case serialization::BodyCentredNonRotatingDynamicFrame::
        kBodyCentredNonRotatingDynamicFrameFieldNumber:
      return CHECK_NOTNULL(plugin)->NewBodyCentredNonRotatingNavigationFrame(
          parameters.centre_index);
    default:
      LOG(FATAL) << "Unexpected extension " << parameters.extension;
      base::noreturn();
  }
}

}  // namespace interface
}  // namespace principia
