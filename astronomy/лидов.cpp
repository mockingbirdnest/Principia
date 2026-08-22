#include "astronomy/лидов.hpp"

namespace principia {
namespace astronomy {
namespace _лидов {
namespace internal {

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
  graph.Plot(ЛидовFrozenLine, {-3.0 / 5.0, 0}, region_boundary_colour);
  switch (grid) {
    case ЛидовGrid::None:
      graph.PlotHorizontalLine(0, region_boundary_colour);
      graph.Plot(
          [](double const c₂) {
            return ЛидовMaximalInclinationLine(0 * Radian, c₂);
          },
          {0, 2.0 / 5.0},
          region_boundary_colour);
      break;
    case ЛидовGrid::MaxEccentricityMinInclination:
      for (int ten_e_max = 1; ten_e_max <= 10; ++ten_e_max) {
        double const e_max = ten_e_max / 10.0;
        graph.Plot(
            [e_max](double const c₂) {
              return ЛидовMaximalEccentricityLine(e_max, c₂);
            },
            ЛидовMaximalEccentricityLineC₂Range(e_max),
            eccentricity_colour);
      }
      for (int i_min_degrees = 0; i_min_degrees <= 80; i_min_degrees += 10) {
        Angle const i_min = i_min_degrees * Degree;
        graph.Plot(
            [i_min](double const c₂) {
              return ЛидовMinimalInclinationLine(i_min, c₂);
            },
            ЛидовMinimalInclinationLineC₂Range(i_min),
            inclination_colour);
      }
      break;
    case ЛидовGrid::MinEccentricityMaxInclination:
      for (int ten_e_min = 1; ten_e_min <= 9; ++ten_e_min) {
        double const e_min = ten_e_min / 10.0;
        graph.Plot(
            [e_min](double const c₂) {
              return ЛидовMinimalEccentricityLeftLine(e_min, c₂);
            },
            ЛидовMinimalEccentricityLeftLineC₂Range(e_min),
            eccentricity_colour);
        graph.PlotVerticalLine(
            ЛидовMinimalEccentricityRightLineC₂(e_min),
            eccentricity_colour,
            {{0, ЛидовMinimalEccentricityRightLineC₁Max(e_min)}});
      }
      for (int i_max_degrees = 0; i_max_degrees <= 90; i_max_degrees += 10) {
        Angle const i_max = i_max_degrees * Degree;
        graph.Plot(
            [i_max](double const c₂) {
              return ЛидовMaximalInclinationLine(i_max, c₂);
            },
            ЛидовMaximalInclinationLineC₂Range(i_max),
            inclination_colour);
      }
      break;
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
