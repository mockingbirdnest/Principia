
#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using geometry::Point;
using geometry::R3Element;
using geometry::Vector;
using numerics::FixedVector;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Quantity;
using quantities::Quotient;

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

template<typename It>
std::string ToMathematica(It begin, It end);

std::string ToMathematica(double const& real);

template<typename T, int size>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector);

template<typename T>
std::string ToMathematica(R3Element<T> const& r3_element);

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity);

template<typename S, typename F>
std::string ToMathematica(Vector<S, F> const& vector);

template<typename V>
std::string ToMathematica(Point<V> const& point);

template<typename F>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom);

template<typename... Types>
std::string ToMathematica(std::tuple<Types...> const& tuple);

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>>
std::string ToMathematica(R ref);

std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements);

// Returns its argument.
std::string ToMathematica(std::string const& str);

// Wraps the string in quotes.
// TODO(egg): escape things properly.
std::string Escape(std::string const& str);

// TODO(phl): This doesn't work well for complex structures like orbital
// elements or trajectory iterators.  Surely we can do better.
template<typename T>
struct RemoveUnit;

template<typename T>
typename RemoveUnit<T>::Unitless ExpressIn(
    typename RemoveUnit<T>::Unit const& unit,
    T const& value);

template<typename V>
typename RemoveUnit<Point<V>>::Unitless ExpressIn(
    typename RemoveUnit<Point<V>>::Unit const& unit,
    Point<V> const& value);

template<typename T>
typename RemoveUnit<std::vector<T>>::Unitless ExpressIn(
    typename RemoveUnit<std::vector<T>>::Unit const& unit,
    std::vector<T> const& values);

}  // namespace internal_mathematica

using internal_mathematica::Apply;
using internal_mathematica::Assign;
using internal_mathematica::Escape;
using internal_mathematica::Export;
using internal_mathematica::ExpressIn;
using internal_mathematica::Option;
using internal_mathematica::PlottableDataset;
using internal_mathematica::ToMathematica;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
