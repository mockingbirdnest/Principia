#pragma once

#include "graphics/graph.hpp"

#include <optional>
#include <vector>

namespace principia {
namespace graphics {
namespace _graph {
namespace internal {

template<affine Abscissa, affine Ordinate>
Graph<Abscissa, Ordinate>::Graph(std::int64_t const width,
                                 std::int64_t const height,
                                 Interval<Abscissa> const& x_range,
                                 Interval<Ordinate> const& y_range,
                                 RGBA32 const background)
    : x_range_(x_range),
      y_range_(y_range),
      width_(width),
      height_(height),
      pixels_(width * height, background) {}

template<affine Abscissa, affine Ordinate>
template<std::ranges::range Points>
void Graph<Abscissa, Ordinate>::ListPointPlot(Points const& points,
                                              RGB24 const colour) {
  for (auto const& [x, y] : points) {
    SetPixel(abscissa_to_pixel(x), ordinate_to_pixel(y), Opaque(colour));
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::Plot(std::function<Ordinate(Abscissa)> const& f,
                                     Interval<Abscissa> const& range,
                                     RGB24 const colour) {
  for (std::int64_t i = abscissa_to_pixel(range.min);
       i <= abscissa_to_pixel(range.max);
       ++i) {
    Interval const pixel_range = Intersection(pixel_to_abscissa(i), range);
    double const f_x_min = f(pixel_range.min);
    double const f_x_max = f(pixel_range.max);
    double f_min;
    double f_max;
    if (f_x_min <= f_x_max) {
      f_min = f_x_min;
      f_max = f_x_max;
    } else {
      f_min = f_x_max;
      f_max = f_x_min;
    }
    for (std::int64_t j = ordinate_to_pixel(f_max);
         j <= ordinate_to_pixel(f_min);
         ++j) {
      SetPixel(i, j, Opaque(colour));
    }
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotVerticalLine(
    Abscissa const x,
    RGB24 const colour,
    std::optional<Interval<Ordinate>> const y_range) {
  for (std::int64_t j = y_range.has_value() ? ordinate_to_pixel(y_range->max)
                                            : 0;
       j < (y_range.has_value() ? ordinate_to_pixel(y_range->min) : height_);
       ++j) {
    SetPixel(abscissa_to_pixel(x), j, Opaque(colour));
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotHorizontalLine(Ordinate y,
                                                          RGB24 colour) {
  for (std::int64_t i = 0; i < width_; ++i) {
    SetPixel(i, ordinate_to_pixel(y), Opaque(colour));
  }
}

template<affine Abscissa, affine Ordinate>
std::vector<RGBA32> const& Graph<Abscissa, Ordinate>::pixels() const {
  return pixels_;
}

template<affine Abscissa, affine Ordinate>
std::int64_t Graph<Abscissa, Ordinate>::width() const {
  return width_;
}

template<affine Abscissa, affine Ordinate>
std::int64_t Graph<Abscissa, Ordinate>::height() const {
  return height_;
}

template<affine Abscissa, affine Ordinate>
std::int64_t Graph<Abscissa, Ordinate>::abscissa_to_pixel(
    Abscissa const x) const {
  return static_cast<std::int64_t>((x - x_range_.min) * inverse_pixel_width_);
}

template<affine Abscissa, affine Ordinate>
Interval<Abscissa> Graph<Abscissa, Ordinate>::pixel_to_abscissa(
    std::int64_t i) const {
  return {i * pixel_width_ + x_range_.min,
          (i + 1) * pixel_width_ + x_range_.min};
}

template<affine Abscissa, affine Ordinate>
std::int64_t Graph<Abscissa, Ordinate>::ordinate_to_pixel(
    Ordinate const y) const {
  return static_cast<std::int64_t>((y_range_.max - y) * inverse_pixel_height_);
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::SetPixel(std::int64_t const column,
                                         std::int64_t const row,
                                         RGBA32 const pixel) {
  if (row < 0 || column < 0 || row >= height_ || column >= width_) {
    return;
  }
  pixels_[row * width_ + column] = pixel;
}

}  // namespace internal
}  // namespace _graph
}  // namespace graphics
}  // namespace principia
