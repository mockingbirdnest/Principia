#pragma once

#include <string>
#include <vector>

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Quantity;

namespace mathematica {

std::string Apply(std::string const& function,
                                std::vector<std::string> const& arguments);

template<typename T>
std::string Option(std::string const& name, T const& right);

template<typename T>
std::string Assign(std::string const& name, T const& right);

std::string Export(std::string const& file, std::string const& expression);

template<typename T, typename U>
std::string PlottableDataset(std::vector<T> x, std::vector<U> y);

template<typename T>
std::string ToMathematica(std::vector<T> list);
std::string ToMathematica(double const& real);
template<typename D>
std::string ToMathematica(Quantity<D> const& quantity);
// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string Escape(std::string const& str);
// Returns its argument.
std::string ToMathematica(std::string const& str);

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
