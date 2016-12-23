
#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using geometry::Point;
using geometry::Vector;
using quantities::Quantity;

std::string Apply(std::string const& function,
                  std::vector<std::string> const& arguments);

template<typename T>
std::string Option(std::string const& name, T const& right);

template<typename T>
std::string Assign(std::string const& name, T const& right);

std::string Export(std::string const& file, std::string const& expression);

template<typename T, typename U>
std::string PlottableDataset(std::vector<T> const& x, std::vector<U> const& y);

template<typename T>
std::string ToMathematica(std::vector<T> const& list);

std::string ToMathematica(double const& real);

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity);

template<typename S, typename F>
std::string ToMathematica(Vector<S, F> const& vector);

template<typename V>
std::string ToMathematica(Point<V> const& point);

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple);

// Returns its argument.
std::string ToMathematica(std::string const& str);

// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string Escape(std::string const& str);

}  // namespace internal_mathematica

using internal_mathematica::Apply;
using internal_mathematica::Assign;
using internal_mathematica::Escape;
using internal_mathematica::Export;
using internal_mathematica::Option;
using internal_mathematica::PlottableDataset;
using internal_mathematica::ToMathematica;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
