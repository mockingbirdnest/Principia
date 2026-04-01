#include "ksp_plugin/interface.hpp"

#include <utility>

#include "glog/logging.h"
#include "journal/method.hpp"
#include "journal/profiles.hpp"  // 🧙 For generated profiles.

namespace principia {
namespace interface {

using namespace principia::quantities::_si;
using namespace principia::quantities::_quantities;
using namespace principia::numerics::_elementary_functions;

static const Angle i_critical = ArcCos(Sqrt(3.0 / 5.0));

double __cdecl principia__GraphLidovMinimalInclinationLine(
    double const inclination_in_degrees,
    double const c₂) {
  journal::Method<journal::GraphLidovMinimalInclinationLine> m(
      {inclination_in_degrees, c₂});
  Angle const i = inclination_in_degrees * Degree;
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return m.Return(cos²_i * (5 * cos²_i - 3 - 5 * c₂) / (5 * cos²_i - 3));
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

}  // namespace interface
}  // namespace principia