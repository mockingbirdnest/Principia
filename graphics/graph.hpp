#pragma once

#include <optional>
#include <ranges>

#include "base/algebra.hpp"
#include "geometry/interval.hpp"

namespace principia {
namespace graphics {
namespace _graph {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_interval;

struct RGB24 {
  std::uint8_t red;
  std::uint8_t green;
  std::uint8_t blue;
};

struct RGBA32 {
  RGB24 colour;
  std::uint8_t alpha;
};

template<affine Abscissa, affine Ordinate>
class Graph {
 public:
  Graph(int width,
        int height,
        Interval<Abscissa> x_range,
        Interval<Ordinate> y_range,
        RGBA32 background);

  template<std::ranges::range Points>
  void ListPointPlot(Points const& points, RGB24 colour);

  void Plot(std::function<Ordinate(Abscissa)> const& f,
            Interval<Abscissa> range,
            RGB24 colour);
  void PlotVerticalLine(
      Abscissa x,
      RGB24 colour,
      std::optional<Interval<Ordinate>> y_range = std::nullopt);
  void PlotHorizontalLine(Ordinate y, RGB24 colour);

  std::vector<RGBA32> const& pixels() const;
  int width() const;
  int height() const;

 private:
  int AbscissaToPixel(Abscissa x) const;
  Interval<Abscissa> PixelToAbscissa(int i) const;
  int OrdinateToPixel(Ordinate y) const;

  void SetPixel(int column, int row, RGBA32 pixel);


  Interval<Abscissa> x_range_;
  Interval<Ordinate> y_range_;

  int width_;
  int height_;
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