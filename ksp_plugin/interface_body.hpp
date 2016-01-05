#pragma once

#include "ksp_plugin/interface.hpp"

namespace principia {
namespace ksp_plugin {

inline bool operator==(XYZ const& left, XYZ const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

inline bool operator==(XYZSegment const& left, XYZSegment const& right) {
  return left.begin == right.begin && left.end == right.end;
}

inline bool operator==(WXYZ const& left, WXYZ const& right) {
  return left.w == right.w && left.x == right.x &&
         left.y == right.y && left.z == right.z;
}

inline bool operator==(QP const& left, QP const& right) {
  return left.q == right.q && left.p == right.p;
}

inline bool operator==(NavigationFrameParameters const& left,
                       NavigationFrameParameters const& right) {
  return left.extension == right.extension &&
         left.centre_index == right.centre_index &&
         left.primary_index == right.primary_index &&
         left.secondary_index == right.secondary_index;
}

}  // namespace ksp_plugin
}  // namespace principia
