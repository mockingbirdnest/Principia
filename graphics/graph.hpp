#pragma once

#include <optional>
#include <ranges>
#include <vector>

#include "base/algebra.hpp"
#include "geometry/interval.hpp"
#include "graphics/colours.hpp"

namespace principia {
namespace graphics {
namespace _graph {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_interval;
using namespace principia::graphics::_colours;

template<typename Point, typename Abscissa, typename Ordinate>
concept graph_point =
    (std::tuple_size_v<Point> == 2) &&
    std::convertible_to<std::tuple_element_t<0, Point>, Abscissa> &&
    std::convertible_to<std::tuple_element_t<1, Point>, Ordinate>;

template<affine Abscissa, affine Ordinate>
class Graph {
 public:
  Graph(std::int64_t width,
        std::int64_t height,
        Interval<Abscissa> const& x_range,
        Interval<Ordinate> const& y_range,
        RGBA32 background);

  template<std::ranges::range Points>
    requires graph_point<std::ranges::range_value_t<Points>, Abscissa, Ordinate>
  void ListPointPlot(Points const& points, RGB24 colour);

  void Plot(std::function<Ordinate(Abscissa)> const& f,
            Interval<Abscissa> const& range,
            RGB24 colour);
  void PlotVerticalLine(
      Abscissa x,
      RGB24 colour,
      std::optional<Interval<Ordinate>> y_range = std::nullopt);
  void PlotHorizontalLine(Ordinate y, RGB24 colour);

  std::vector<RGBA32> const& pixels() const;
  std::int64_t width() const;
  std::int64_t height() const;

 private:
  std::int64_t abscissa_to_pixel(Abscissa x) const;
  Interval<Abscissa> pixel_to_abscissa(std::int64_t i) const;
  // This function is decreasing: `y_range_.max` maps to 0.
  std::int64_t ordinate_to_pixel(Ordinate y) const;

  void SetPixel(std::int64_t column, std::int64_t row, RGBA32 pixel);


  Interval<Abscissa> const x_range_;
  Interval<Ordinate> const y_range_;

  std::int64_t const width_;
  std::int64_t const height_;
  Difference<Abscissa> const pixel_width_ = x_range_.measure() / width_;
  Difference<Ordinate> const pixel_height_ = y_range_.measure() / height_;
  Inverse<Difference<Abscissa>> const inverse_pixel_width_ = 1 / pixel_width_;
  Inverse<Difference<Ordinate>> const inverse_pixel_height_ = 1 / pixel_height_;
  std::vector<RGBA32> pixels_;
};

}  // namespace internal

using internal::Graph;
using internal::RGBA32;
using internal::RGB24;

}  // namespace _graph
}  // namespace graphics
}  // namespace principia

#include "graphics/graph_body.hpp"
