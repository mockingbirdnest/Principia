#include "ksp_plugin/interface.hpp"

#include <utility>

#include "absl/log/log.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace interface {

using namespace principia::journal::_method;
using namespace principia::quantities::_si;
using namespace principia::quantities::_quantities;
using namespace principia::numerics::_elementary_functions;

namespace {

static const Angle i_critical = ArcCos(Sqrt(3.0 / 5.0));

}  // namespace

double __cdecl principia__GraphLidovFrozenLine(double const c₂) {
  journal::Method<journal::GraphLidovFrozenLine> m({c₂});
  return m.Return(3.0 / 5.0 - 2 * Sqrt(-3.0 / 5.0 * c₂) - c₂);
}

double __cdecl principia__GraphLidovMaximalEccentricityLine(double const e,
                                                            double const c₂) {
  journal::Method<journal::GraphLidovMaximalEccentricityLine> m({e, c₂});
  double const e² = Pow<2>(e);
  return m.Return(3.0 / 5.0 - c₂ + c₂ / e² - 3 * e² / 5);
}

Interval __cdecl principia__GraphLidovMaximalEccentricityLineC2Range(
    double const e) {
  journal::Method<journal::GraphLidovMaximalEccentricityLineC2Range> m({e});
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return m.Return({-3 * e⁴ / 5, 2 * e² / 5});
}

double __cdecl principia__GraphLidovMaximalInclinationLine(
    double const inclination_in_degrees,
    double const c₂) {
  journal::Method<journal::GraphLidovMaximalInclinationLine> m(
      {inclination_in_degrees, c₂});
  Angle const i = inclination_in_degrees * Degree;
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return m.Return(c₂ < 0 ? cos²_i * (5 * cos²_i - 5 * c₂ - 3) / (5 * cos²_i - 3)
                         : (2 - 5 * c₂) * cos²_i / 2);
}

Interval __cdecl principia__GraphLidovMaximalInclinationLineC2Range(
    double const inclination_in_degrees) {
  journal::Method<journal::GraphLidovMaximalInclinationLineC2Range> m(
      {inclination_in_degrees});
  Angle const i = inclination_in_degrees * Degree;
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return m.Return(
      {i > i_critical ? -Pow<2>(1 - 5 * Cos(2 * i)) / 60 : 0, 2.0 / 5.0});
}

double __cdecl principia__GraphLidovMinimalInclinationLine(
    double const inclination_in_degrees,
    double const c₂) {
  journal::Method<journal::GraphLidovMinimalInclinationLine> m(
      {inclination_in_degrees, c₂});
  Angle const i = inclination_in_degrees * Degree;
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return m.Return(cos²_i * (5 * cos²_i - 5 * c₂ - 3) / (5 * cos²_i - 3));
}

Interval __cdecl principia__GraphLidovMinimalInclinationLineC2Range(
    double const inclination_in_degrees) {
  journal::Method<journal::GraphLidovMinimalInclinationLineC2Range> m(
      {inclination_in_degrees});
  Angle const i = inclination_in_degrees * Degree;
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return m.Return(i > i_critical ? Interval{cos²_i - 3.0 / 5.0,
                                            -Pow<2>(1 - 5 * Cos(2 * i)) / 60}
                                 : Interval{0, cos²_i - 3.0 / 5.0});
}

double __cdecl principia__GraphLidovMinimalEccentricityLeftLine(double const e,
                                                            double const c₂) {
  journal::Method<journal::GraphLidovMinimalEccentricityLeftLine> m({e, c₂});
  double const e² = Pow<2>(e);
  return m.Return(3.0 / 5.0 - c₂ + c₂ / e² - 3 * e² / 5);
}

Interval __cdecl principia__GraphLidovMinimalEccentricityLeftLineC2Range(
    double const e) {
  journal::Method<journal::GraphLidovMinimalEccentricityLeftLineC2Range> m({e});
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return m.Return({-3 * e² / 5, -3 * e⁴ / 5});
}

void __cdecl principia__GraphLidovMinimalEccentricityRightLineC2AndC1Max(
    double const e,
    double* const c₂,
    double* const c₁_max) {
  journal::Method<journal::GraphLidovMinimalEccentricityRightLineC2AndC1Max> m(
      {e}, {c₂, c₁_max});
  double const e² = Pow<2>(e);
  *c₂ = 2 * e² / 5;
  *c₁_max = 1 - e²;
  return m.Return();
}

}  // namespace interface
}  // namespace principia
