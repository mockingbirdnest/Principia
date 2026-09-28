#include "astronomy/лидов.hpp"

namespace principia {
namespace astronomy {
namespace _лидов {
namespace internal {

// All functions in this file refer to an orbit perturbed as in the analysis of
// [Лид61].  The parameters c₁ and c₂ are as defined there.

Angle const i_critical = ArcCos(Sqrt(3.0 / 5.0));

// Returns c₁ such that an orbit with these values of c₁ and c₂ has no
// eccentricity-inclination exchange.
double FrozenLine(double const c₂) {
  CHECK_LE(c₂, 0);
  return 3.0 / 5.0 - 2 * Sqrt(-3.0 / 5.0 * c₂) - c₂;
}

// Returns c₁ such that the upper bound of eccentricity for an orbit with these
// values of c₁ and c₂ is e.
double MaximalEccentricityLine(double const e, double const c₂) {
  double const e² = Pow<2>(e);
  return 3.0 / 5.0 - c₂ + c₂ / e² - 3 * e² / 5.0;
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of eccentricity for an orbit with these values of c₁ and c₂
// is e.
Interval<double> MaximalEccentricityLineC₂Range(double const e) {
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return {-3.0 * e⁴ / 5.0, 2.0 * e² / 5.0};
}

// Returns c₁ such that the upper bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
double MaximalInclinationLine(Angle const i, double const c₂) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return c₂ < 0
             ? cos²_i * (5.0 * cos²_i - 5.0 * c₂ - 3.0) / (5.0 * cos²_i - 3.0)
             : (2.0 - 5.0 * c₂) * cos²_i / 2.0;
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// upper bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> MaximalInclinationLineC₂Range(Angle const i) {
  double const cos_i = Cos(i);
  return {i > i_critical ? -Pow<2>(1.0 - 5.0 * Cos(2.0 * i)) / 60.0 : 0,
          2.0 / 5.0};
}

// Returns c₁ such that the lower bound of inclination for an orbit with these
// values of c₁ and c₂ is i.
double MinimalInclinationLine(Angle const i, double const c₂) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return cos²_i * (5.0 * cos²_i - 5.0 * c₂ - 3.0) / (5.0 * cos²_i - 3.0);
}

// Returns the range of values of c₂ such that there exists a c₁ such that the
// lower bound of inclination for an orbit with these values of c₁ and c₂
// is i.
Interval<double> MinimalInclinationLineC₂Range(Angle const i) {
  double const cos_i = Cos(i);
  double const cos²_i = Pow<2>(cos_i);
  return i > i_critical
             ? Interval<double>{cos²_i - 3.0 / 5.0,
                                -Pow<2>(1.0 - 5.0 * Cos(2 * i)) / 60.0}
             : Interval<double>{0, cos²_i - 3.0 / 5.0};
}

// Returns the value of c₁ such that the lower bound of eccentricity for an
// orbit with these values of c₁ and c₂ is e.
double MinimalEccentricityLeftLine(double const e, double const c₂) {
  double const e² = Pow<2>(e);
  return 3.0 / 5.0 - c₂ + c₂ / e² - 3.0 * e² / 5.0;
}

// Returns the range of negative values of c₂ such that there exists a c₁ such
// that the lower bound of eccentricity for an orbit with these values of c₁
// and c₂ is e.
Interval<double> MinimalEccentricityLeftLineC₂Range(double const e) {
  double const e² = Pow<2>(e);
  double const e⁴ = Pow<4>(e);
  return {-3.0 * e² / 5.0, -3.0 * e⁴ / 5.0};
}

// Returns the positive value of c₂ for which the lower bound of eccentricity is
// e.
double MinimalEccentricityRightLineC₂(double const e) {
  double const e² = Pow<2>(e);
  return 2.0 * e² / 5;
}

// Returns the maximal possible value of c₁ that can be attained when c₂ has the
// positive value for which the lower bound of eccentricity is e.
double MinimalEccentricityRightLineC₁Max(double const e) {
  double const e² = Pow<2>(e);
  return 1.0 - e²;
}

Graph<double, double> ЛидовGraph(OrbitalElements const& elements,
                                 std::int64_t const width,
                                 std::int64_t const height,
                                 RGBA32 const background,
                                 RGB24 const region_boundary_colour,
                                 RGB24 const inclination_colour,
                                 RGB24 const eccentricity_colour,
                                 RGB24 const лидов_parameter_colour,
                                 ЛидовGrid const grid) {
  Graph<double, double> graph(
      width, height, {-3.0 / 5.0, 2.0 / 5.0}, {0, 1}, background);
  graph.PlotVerticalLine(0, region_boundary_colour);
  graph.Plot(FrozenLine, {-3.0 / 5.0, 0}, region_boundary_colour);
  switch (grid) {
    case ЛидовGrid::None:
      graph.PlotHorizontalLine(0, region_boundary_colour);
      graph.Plot(
          [](double const c₂) {
            return MaximalInclinationLine(0 * Radian, c₂);
          },
          {0, 2.0 / 5.0},
          region_boundary_colour);
      break;
    case ЛидовGrid::MaxEccentricityMinInclination:
      for (int ten_e_max = 1; ten_e_max <= 10; ++ten_e_max) {
        double const e_max = ten_e_max / 10.0;
        graph.Plot(
            [e_max](double const c₂) {
              return MaximalEccentricityLine(e_max, c₂);
            },
            MaximalEccentricityLineC₂Range(e_max),
            eccentricity_colour);
      }
      for (int i_min_degrees = 0; i_min_degrees <= 80; i_min_degrees += 10) {
        Angle const i_min = i_min_degrees * Degree;
        graph.Plot(
            [i_min](double const c₂) {
              return MinimalInclinationLine(i_min, c₂);
            },
            MinimalInclinationLineC₂Range(i_min),
            inclination_colour);
      }
      for (int i_min_degrees = 10; i_min_degrees <= 30; i_min_degrees += 10) {
        Interval const c₂ =
            MinimalInclinationLineC₂Range(i_min_degrees * Degree);
        graph.AddLabel(std::pair{c₂.max, 0.0},
                       absl::StrCat(i_min_degrees, "°"),
                       inclination_colour,
                       Label::TextPlacement::Below);
      }
      for (int i_min_degrees = 40; i_min_degrees <= 80; i_min_degrees += 10) {
        Interval const c₂ =
            MinimalInclinationLineC₂Range(i_min_degrees * Degree);
        graph.AddLabel(std::pair{c₂.min, 0.0},
                       absl::StrCat(i_min_degrees, "°"),
                       inclination_colour,
                       Label::TextPlacement::Below);
      }
      break;
    case ЛидовGrid::MinEccentricityMaxInclination:
      for (int ten_e_min = 1; ten_e_min <= 9; ++ten_e_min) {
        double const e_min = ten_e_min / 10.0;
        graph.Plot(
            [e_min](double const c₂) {
              return MinimalEccentricityLeftLine(e_min, c₂);
            },
            MinimalEccentricityLeftLineC₂Range(e_min),
            eccentricity_colour);
        graph.PlotVerticalLine(
            MinimalEccentricityRightLineC₂(e_min),
            eccentricity_colour,
            {{0, MinimalEccentricityRightLineC₁Max(e_min)}});
      }
      for (int i_max_degrees = 0; i_max_degrees <= 90; i_max_degrees += 10) {
        Angle const i_max = i_max_degrees * Degree;
        graph.Plot(
            [i_max](double const c₂) {
              return MaximalInclinationLine(i_max, c₂);
            },
            MaximalInclinationLineC₂Range(i_max),
            inclination_colour);
      }
      for (int ten_e_min = 4; ten_e_min <= 9; ++ten_e_min) {
        double e_min = ten_e_min / 10.0;
        Interval c₂ = MinimalEccentricityLeftLineC₂Range(e_min);
        graph.AddLabel(std::pair{c₂.min, 0.0},
                       absl::StrCat(".", ten_e_min),
                       eccentricity_colour,
                       Label::TextPlacement::Below);
      }
      for (int ten_e_min = 6; ten_e_min <= 9; ++ten_e_min) {
        double const e_min = ten_e_min / 10.0;
        double const c₂ = MinimalEccentricityRightLineC₂(e_min);
        graph.AddLabel(std::pair{c₂, 0.0},
                       absl::StrCat(".", ten_e_min),
                       eccentricity_colour,
                       Label::TextPlacement::Below);
      }
      break;
  }
  // The inclination labels on the frozen curve are the same for both the max
  // and min lines (because the inclination is frozen there); it is easiest to
  // position them based on the maximal inclination lines (because they are
  // then at the minimal c₂ for the line, instead of at either end of the c₂
  // interval depending on the inclination).  Likewise the eccentricity labels
  // on the equatorial curve are the for both max and min e.
  if (grid != ЛидовGrid::None) {
    for (int i_degrees = 10; i_degrees <= 60; i_degrees += 10) {
      Angle const i = i_degrees * Degree;
      Interval c₂ =
          MaximalInclinationLineC₂Range(i);
      double c₁ =
          MaximalInclinationLine(i, c₂.min);
      graph.AddLabel(std::pair{c₂.min, c₁},
                     absl::StrCat(i_degrees, "°"),
                     inclination_colour,
                     Label::TextPlacement::Left);
    }
    for (int ten_e = 2; ten_e <= 9; ++ten_e) {
      double e = ten_e / 10.0;
      double const c₂ = MinimalEccentricityRightLineC₂(e);
      double const c₁ = MinimalEccentricityRightLineC₁Max(e);
      graph.AddLabel(std::pair{c₂, c₁},
                     absl::StrCat(" .", ten_e),
                     eccentricity_colour,
                     Label::TextPlacement::Right);
    }
  }
  graph.ListPointPlot(
      elements.mean_elements() |
          std::ranges::views::transform(
              [](OrbitalElements::ClassicalElements const& elements) {
                auto const [sin_i, cos_i] = SinCos(elements.inclination);
                auto const sin_ω = Sin(elements.argument_of_periapsis);
                double const sin²_i = Pow<2>(sin_i);
                double const cos²_i = Pow<2>(cos_i);
                double const& e = elements.eccentricity;
                double const e² = Pow<2>(e);
                double const sin²_ω = Pow<2>(sin_ω);
                double const c₂ = e² * (2.0 / 5.0 - sin²_i * sin²_ω);
                double const c₁ = (1 - e²) * cos²_i;
                return std::pair{c₂, c₁};
              }),
      лидов_parameter_colour);
  return graph;
}

}  // namespace internal
}  // namespace _лидов
}  // namespace astronomy
}  // namespace principia
