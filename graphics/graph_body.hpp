#pragma once

#include "graphics/graph.hpp"

#include <algorithm>
#include <optional>
#include <string>
#include <utility>
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
                                 RGBA32 const background,
                                 Graph::Options const options)
    : x_range_(x_range),
      y_range_(y_range),
      width_(width),
      height_(height),
      anti_aliasing_(options.anti_aliasing),
      samples_(samples_width_ * samples_height_, background) {}

template<affine Abscissa, affine Ordinate>
template<std::ranges::range Points>
  requires graph_point<std::ranges::range_value_t<Points>, Abscissa, Ordinate>
void Graph<Abscissa, Ordinate>::ListPointPlot(
    Points const& points,
    RGB24 const colour,
    ListPointPlotOptions const options) {
  InvalidateAntiAliasing();
  if (options.joined) {
    std::optional<Abscissa> xᵢ₋₁;
    std::optional<Ordinate> yᵢ₋₁;
    for (auto const& [xᵢ, yᵢ] : points) {
      if (xᵢ₋₁.has_value()) {
        if (*xᵢ₋₁ == xᵢ) {
          PlotVerticalLine(
              xᵢ, colour, {{std::min(*yᵢ₋₁, yᵢ), std::max(*yᵢ₋₁, yᵢ)}});
        } else {
          auto const slope = (yᵢ - *yᵢ₋₁) / (xᵢ - *xᵢ₋₁);
          Plot([xᵢ₋₁ = *xᵢ₋₁, yᵢ₋₁ = *yᵢ₋₁, slope](
                   Abscissa const& x) { return yᵢ₋₁ + (x - xᵢ₋₁) * slope; },
               {std::min(*xᵢ₋₁, xᵢ), std::max(*xᵢ₋₁, xᵢ)},
               colour);
        }
      }
      xᵢ₋₁ = xᵢ;
      yᵢ₋₁ = yᵢ;
    }
  } else {
    for (auto const& [x, y] : points) {
      DrawRectangle(x, y, pixel_width_, pixel_height_, colour);
    }
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::Plot(std::function<Ordinate(Abscissa)> const& f,
                                     Interval<Abscissa> const& range,
                                     RGB24 const colour) {
  InvalidateAntiAliasing();
  for (std::int64_t i = abscissa_to_sample(range.min);
       i <= abscissa_to_sample(range.max);
       ++i) {
    Interval const sample_range = Intersection(sample_to_abscissa_range(i), range);
    Ordinate const f_x_min = f(sample_range.min);
    Ordinate const f_x_max = f(sample_range.max);
    Ordinate f_min;
    Ordinate f_max;
    if (f_x_min <= f_x_max) {
      f_min = f_x_min;
      f_max = f_x_max;
    } else {
      f_min = f_x_max;
      f_max = f_x_min;
    }
    for (auto y = f_min; y <= f_max; y += sample_height_) {
      DrawRectangle(sample_range.midpoint(), y, pixel_width_, pixel_height_, colour);
    }
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotVerticalLine(
    Abscissa const x,
    RGB24 const colour,
    std::optional<Interval<Ordinate>> const y_range) {
  InvalidateAntiAliasing();
  auto const& range = y_range.value_or(y_range_);
  DrawRectangle(x,
                range.midpoint(),
                pixel_width_,
                range.measure() + pixel_height_,
                colour);
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::PlotHorizontalLine(Ordinate const y,
                                                   RGB24 const colour) {
  InvalidateAntiAliasing();
  DrawRectangle(x_range_.midpoint(),
                y,
                x_range_.measure() + pixel_width_,
                pixel_height_,
                colour);
}

template<affine Abscissa, affine Ordinate>
template<graph_point<Abscissa, Ordinate> Point>
void Graph<Abscissa, Ordinate>::AddLabel(Point const& point,
                                         std::string text,
                                         RGB24 const colour,
                                         Label::TextPlacement const placement) {
  auto const& [x, y] = point;
  labels_.emplace_back(abscissa_to_sample(x) / anti_aliasing_.supersampling,
                       ordinate_to_sample(y) / anti_aliasing_.supersampling,
                       std::move(text),
                       colour,
                       placement);
}

template<affine Abscissa, affine Ordinate>
inline void Graph<Abscissa, Ordinate>::AntiAlias() {
  if (anti_aliasing_.supersampling == 1) {
    return;
  }
  pixels_.emplace(width_ * height_);
  std::int64_t const sample_count =
      anti_aliasing_.supersampling * anti_aliasing_.supersampling;
  for (std::int64_t row = 0; row < height_; ++row) {
    for (std::int64_t column = 0; column < width_; ++column) {
      // TODO(egg): This is bad colour theory.
      std::int64_t r_total = 0;
      std::int64_t g_total = 0;
      std::int64_t b_total = 0;
      std::int64_t a_total = 0;
      for (std::int64_t sub_row = 0; sub_row < anti_aliasing_.supersampling;
           ++sub_row) {
        for (std::int64_t sub_column = 0;
             sub_column < anti_aliasing_.supersampling;
             ++sub_column) {
          RGBA32 const& sample =
              samples_[(row * anti_aliasing_.supersampling + sub_row) *
                           samples_width_ +
                       column * anti_aliasing_.supersampling + sub_column];
          r_total += sample.colour.red;
          g_total += sample.colour.green;
          b_total += sample.colour.blue;
          a_total += sample.alpha;
        }
      }
      (*pixels_)[row * width_ + column] = {
          {static_cast<std::uint8_t>(r_total / sample_count),
           static_cast<std::uint8_t>(g_total / sample_count),
           static_cast<std::uint8_t>(b_total / sample_count)},
          static_cast<std::uint8_t>(a_total / sample_count)};
    }
  }
}

template<affine Abscissa, affine Ordinate>
std::vector<RGBA32> const& Graph<Abscissa, Ordinate>::pixels() const {
  if (anti_aliasing_.supersampling == 1) {
    return samples_;
  }
  return *pixels_;
}

template<affine Abscissa, affine Ordinate>
std::vector<Label> const& Graph<Abscissa, Ordinate>::labels() const {
  return labels_;
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
std::int64_t Graph<Abscissa, Ordinate>::abscissa_to_sample(
    Abscissa const x) const {
  return static_cast<std::int64_t>((x - x_range_.min) * inverse_sample_width_ +
                                   0.5 * anti_aliasing_.supersampling);
}

template<affine Abscissa, affine Ordinate>
Interval<Abscissa> Graph<Abscissa, Ordinate>::sample_to_abscissa_range(
    std::int64_t i) const {
  return {
      (i - 0.5 * anti_aliasing_.supersampling) * sample_width_ + x_range_.min,
      (i + 0.5 * anti_aliasing_.supersampling) * sample_width_ + x_range_.min};
}

template<affine Abscissa, affine Ordinate>
std::int64_t Graph<Abscissa, Ordinate>::ordinate_to_sample(
    Ordinate const y) const {
  return static_cast<std::int64_t>((y_range_.max - y) * inverse_sample_height_ +
                                   0.5 * anti_aliasing_.supersampling);
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::DrawRectangle(Abscissa const x,
                                          Ordinate const y,
                     Difference<Abscissa> width,
                     Difference<Ordinate> height,
                     RGB24 colour) {
  std::int64_t const first_row =
      abscissa_to_sample(x - (width - sample_width_) * 0.5);
  std::int64_t const last_row =
      abscissa_to_sample(x + (width - sample_width_) * 0.5);
  // Note the signs: `ordinate_to_sample` is decreasing.
  std::int64_t const first_column =
      ordinate_to_sample(y + (height - sample_height_) * 0.5);
  std::int64_t const last_column =
      ordinate_to_sample(y - (height - sample_height_) * 0.5);
  for (std::int64_t row = first_row; row <= last_row; ++row) {
    for (std::int64_t column = first_column; column <= last_column; ++column) {
      SetSample(row, column, Opaque(colour));
    }
  }
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::DrawDisk(Abscissa x,
                                         Ordinate y,
                                         double r,
                                         RGB24 colour) {}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::SetSample(std::int64_t const column,
                                          std::int64_t const row,
                                          RGBA32 const sample) {
  if (row < 0 || column < 0 || row >= samples_height_ ||
      column >= samples_width_) {
    return;
  }
  samples_[row * samples_width_ + column] = sample;
}

template<affine Abscissa, affine Ordinate>
void Graph<Abscissa, Ordinate>::InvalidateAntiAliasing() {
  pixels_.reset();
}

}  // namespace internal
}  // namespace _graph
}  // namespace graphics
}  // namespace principia
