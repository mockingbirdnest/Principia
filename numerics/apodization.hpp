#pragma once

#include "geometry/instant.hpp"
#include "numerics/poisson_series.hpp"

// The order and terminology in this file follows
// https://en.wikipedia.org/wiki/Window_function.

namespace principia {
namespace numerics {
namespace _apodization {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::numerics::_poisson_series;

// ISO 18431-2:2004, section 5.4.
PoissonSeries<double, 0, 0> Dirichlet(Instant const& t_min,
                                      Instant const& t_max);

PoissonSeries<double, 0, 0> Sine(Instant const& t_min,
                                 Instant const& t_max);

// ISO 18431-2:2004, section 5.2.
PoissonSeries<double, 0, 0> Hann(Instant const& t_min,
                                 Instant const& t_max);

PoissonSeries<double, 0, 0> Hamming(Instant const& t_min,
                                    Instant const& t_max);

PoissonSeries<double, 0, 0> Blackman(Instant const& t_min,
                                     Instant const& t_max);

PoissonSeries<double, 0, 0> ExactBlackman(Instant const& t_min,
                                          Instant const& t_max);

PoissonSeries<double, 0, 0> Nuttall(Instant const& t_min,
                                    Instant const& t_max);

PoissonSeries<double, 0, 0> BlackmanNuttall(Instant const& t_min,
                                            Instant const& t_max);

PoissonSeries<double, 0, 0> BlackmanHarris(Instant const& t_min,
                                           Instant const& t_max);

// The flat-top window in Wikipedia is not normalized and comes from Matlab (?).
// We use the normalized function from ISO 18431-2:2004/Cor.1:2008, section 5.3
// instead.
PoissonSeries<double, 0, 0> ISO18431_2(Instant const& t_min,
                                       Instant const& t_max);

}  // namespace internal

using internal::Blackman;
using internal::BlackmanHarris;
using internal::BlackmanNuttall;
using internal::Dirichlet;
using internal::ExactBlackman;
using internal::Hann;
using internal::Hamming;
using internal::ISO18431_2;
using internal::Nuttall;
using internal::Sine;

}  // namespace _apodization
}  // namespace numerics
}  // namespace principia
