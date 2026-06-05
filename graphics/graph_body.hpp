#pragma once
#include "graphics/graph.hpp"

#include <optional>

namespace principia {
namespace graphics {
namespace _graph {
namespace internal {

template<affine Abscissa, affine Ordinate>
Graph<Abscissa, Ordinate>::Graph(int const width,
                                 int const height,
                                 Interval<Abscissa> const x_range,
                                 Interval<Ordinate> const y_range,
                                 RGBA32 const background)
    : x_range_(x_range),
      y_range_(y_range),
      width_(width),
      height_(height),
      pixels_(width * height, background) {}

template<affine Abscissa, affine Ordinate>
template<std::ranges::range Points>
void Graph<Abscissa, Ordinate>::ListPointPlot(Points const& points,
                                             RGB24 colour) {
  for (auto const& [x, y] : points) {
    SetPixel(AbscissaToPixel(x),
             OrdinateToPixel(y),
             {.colour = colour, .alpha = 255});
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::Plot(std::function<Ordinate(Abscissa)> const& f,
                                     Interval<Abscissa> range,
                                     RGB24 colour) {
  for (int i = AbscissaToPixel(range.min); i <= AbscissaToPixel(range.max);
       ++i) {
    Interval pixel_range = PixelToAbscissa(i).IntersectedWith(range);
    double f_x_min = f(pixel_range.min);
    double f_x_max = f(pixel_range.max);
    double f_min;
    double f_max;
    if (f_x_min <= f_x_max) {
      f_min = f_x_min;
      f_max = f_x_max;
    } else {
      f_min = f_x_max;
      f_max = f_x_min;
    }
    for (int j = OrdinateToPixel(f_max); j <= OrdinateToPixel(f_min); ++j) {
      SetPixel(i, j, {.colour = colour, .alpha = 255});
    }
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotVerticalLine(
    Abscissa x,
    RGB24 colour,
    std::optional<Interval<Ordinate>> y_range) {
  for (int j = y_range.has_value() ? OrdinateToPixel(y_range->max) : 0;
       j < (y_range.has_value() ? OrdinateToPixel(y_range->min) : height_);
       ++j) {
    SetPixel(AbscissaToPixel(x), j, {.colour = colour, .alpha = 255});
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotHorizontalLine(Ordinate y,
                                                          RGB24 colour) {
  for (int i = 0; i < width_; ++i) {
    SetPixel(i, OrdinateToPixel(y), {.colour = colour, .alpha = 255});
  }
}

template<affine Abscissa, affine Ordinate>
std::vector<RGBA32> const& Graph<Abscissa, Ordinate>::pixels() const {
  return pixels_;
}

template<affine Abscissa, affine Ordinate>
int Graph<Abscissa, Ordinate>::width() const {
  return width_;
}

template<affine Abscissa, affine Ordinate>
int Graph<Abscissa, Ordinate>::height() const {
  return height_;
}

template<affine Abscissa, affine Ordinate>
int Graph<Abscissa, Ordinate>::AbscissaToPixel(Abscissa x) const {
  return static_cast<int>(width_ * (x - x_range_.min) / x_range_.measure());
}

template<affine Abscissa, affine Ordinate>
Interval<Abscissa> Graph<Abscissa, Ordinate>::PixelToAbscissa(int i) const {
  return {i * x_range_.measure() / width_ + x_range_.min,
          (i + 1) * x_range_.measure() / width_ + x_range_.min};
}

template<affine Abscissa, affine Ordinate>
int Graph<Abscissa, Ordinate>::OrdinateToPixel(Ordinate y) const {
  return static_cast<int>(height_ * (y_range_.max - y) /
                          y_range_.measure());
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::SetPixel(int const column,
                                         int const row,
                                         RGBA32 const pixel) {
  if (row < 0 || column < 0 || row >= height_ || column >= width_) {
    return;
  }
  pixels_[row * width_ + column] = pixel;
}

}
}  // namespace _graph
}  // namespace graphics
}  // namespace principia
