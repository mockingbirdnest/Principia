#pragma once

#include <string>
#include <vector>

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Quantity;

namespace testing_utilities {

template<typename T>
std::string MathematicaFunction(std::string const& function,
                                std::vector<T> const& arguments);

std::string MathematicaFunction(std::string const& function,
                                std::vector<std::string> const& arguments);

template<typename T>
std::string MathematicaOption(std::string const& name, T const& right);

template<typename T>
std::string MathematicaAssign(std::string const& name, T const& right);

template<typename T, typename U>
std::string MathematicaPlottableDataset(std::vector<T> x, std::vector<U> y);

template<typename T>
std::string ToMathematica(std::vector<T> list);
std::string ToMathematica(double const& real);
template<typename D>
std::string ToMathematica(Quantity<D> const& quantity);
// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string MathematicaEscape(std::string const& str);
// Returns its argument.
std::string ToMathematica(std::string const& str);

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/mathematica_body.hpp"
