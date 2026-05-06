#include "ksp_plugin/interface.hpp"

#include <utility>

#include "absl/log/log.h"
#include "astronomy/лидов.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.

namespace principia {
namespace interface {

using namespace principia::astronomy::_лидов;

double __cdecl principia__GraphLidovFrozenLine(double const c₂) {
  journal::Method<journal::GraphLidovFrozenLine> m({c₂});
  return m.Return(ЛидовFrozenLine(c₂));
}

double __cdecl principia__GraphLidovMaximalEccentricityLine(double const e,
                                                            double const c₂) {
  journal::Method<journal::GraphLidovMaximalEccentricityLine> m({e, c₂});
  return m.Return(ЛидовMaximalEccentricityLine(e, c₂));
}

Interval __cdecl principia__GraphLidovMaximalEccentricityLineC2Range(
    double const e) {
  journal::Method<journal::GraphLidovMaximalEccentricityLineC2Range> m({e});
  auto const [c₂_min, c₂_max] = ЛидовMaximalEccentricityLineC₂Range(e);
  return m.Return({c₂_min, c₂_max});
}

double __cdecl principia__GraphLidovMaximalInclinationLine(
    double const inclination_in_degrees,
    double const c₂) {
  journal::Method<journal::GraphLidovMaximalInclinationLine> m(
      {inclination_in_degrees, c₂});
  Angle const i = inclination_in_degrees * Degree;
  return m.Return(ЛидовMaximalInclinationLine(i, c₂));
}

Interval __cdecl principia__GraphLidovMaximalInclinationLineC2Range(
    double const inclination_in_degrees) {
  journal::Method<journal::GraphLidovMaximalInclinationLineC2Range> m(
      {inclination_in_degrees});
  Angle const i = inclination_in_degrees * Degree;
  auto const [c₂_min, c₂_max] = ЛидовMaximalInclinationLineC₂Range(i);
  return m.Return({c₂_min, c₂_max});
}

double __cdecl principia__GraphLidovMinimalInclinationLine(
    double const inclination_in_degrees,
    double const c₂) {
  journal::Method<journal::GraphLidovMinimalInclinationLine> m(
      {inclination_in_degrees, c₂});
  Angle const i = inclination_in_degrees * Degree;
  return m.Return(ЛидовMinimalInclinationLine(i, c₂));
}

Interval __cdecl principia__GraphLidovMinimalInclinationLineC2Range(
    double const inclination_in_degrees) {
  journal::Method<journal::GraphLidovMinimalInclinationLineC2Range> m(
      {inclination_in_degrees});
  Angle const i = inclination_in_degrees * Degree;
  auto const [c₂_min, c₂_max] = ЛидовMinimalInclinationLineC₂Range(i);
  return m.Return({c₂_min, c₂_max});
}

double __cdecl principia__GraphLidovMinimalEccentricityLeftLine(
    double const e,
    double const c₂) {
  journal::Method<journal::GraphLidovMinimalEccentricityLeftLine> m({e, c₂});
  return m.Return(ЛидовMinimalEccentricityLeftLine(e, c₂));
}

Interval __cdecl principia__GraphLidovMinimalEccentricityLeftLineC2Range(
    double const e) {
  journal::Method<journal::GraphLidovMinimalEccentricityLeftLineC2Range> m({e});
  auto const [c₂_min, c₂_max] = ЛидовMinimalEccentricityLeftLineC₂Range(e);
  return m.Return({c₂_min, c₂_max});
}

void __cdecl principia__GraphLidovMinimalEccentricityRightLineC2AndC1Max(
    double const e,
    double* const c₂,
    double* const c₁_max) {
  journal::Method<journal::GraphLidovMinimalEccentricityRightLineC2AndC1Max> m(
      {e}, {c₂, c₁_max});
  CHECK(c₂ != nullptr);
  CHECK(c₁_max != nullptr);
  *c₂ = ЛидовMinimalEccentricityRightLineC₂(e);
  *c₁_max = ЛидовMinimalEccentricityRightLineC₁Max(e);
  return m.Return();
}

}  // namespace interface
}  // namespace principia
