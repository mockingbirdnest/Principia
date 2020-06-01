
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "base/file.hpp"
#include "base/traits.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/fixed_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using base::all_different_v;
using base::OFStream;
using geometry::Bivector;
using geometry::Point;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::Vector;
using numerics::FixedVector;
using physics::DegreesOfFreedom;
using quantities::Amount;
using quantities::Angle;
using quantities::Current;
using quantities::Exponentiation;
using quantities::is_tuple_v;
using quantities::Length;
using quantities::LuminousIntensity;
using quantities::Mass;
using quantities::Pow;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Temperature;
using quantities::Time;

// Define this value to 1 to force the logger to append "_new" to the file
// names, which is useful for regression testing of the logger.
#define PRINCIPIA_MATHEMATICA_LOGGER_REGRESSION_TEST 1

// A helper class for type erasure of quantities.  It may be used with the
// functions in this file to remove the dimensions of quantities (we know that
// Mathematica is sluggish when processing quantities).  Usage:
//
//   ToMathematica(... , ExpressIn(Metre, Second, Degree));
//
// The construction parameters must be values of distinct SI base quantities
// (or Angle). They define a system of units.  They may be in any order.  If the
// other arguments of the functions contain quantities that are not spanned by
// that system of units, the call is ill-formed.
template<typename... Qs>
class ExpressIn {
 public:
  // Check that only SI base quantities or Angle are specified.
  static_assert(
      ((std::is_same_v<Qs, Length> || std::is_same_v<Qs, Mass> ||
        std::is_same_v<Qs, Time> || std::is_same_v<Qs, Current> ||
        std::is_same_v<Qs, Temperature> || std::is_same_v<Qs, Amount> ||
        std::is_same_v<Qs, LuminousIntensity> || std::is_same_v<Qs, Angle>) &&
       ...), "Must instantiate with SI base quantities or Angle");

  // Check that all quantities are different.
  static_assert(all_different_v<Qs...>, "Must use different quantities");

  ExpressIn(Qs const&... qs);  // NOLINT(runtime/explicit)

  template<typename Q>
  double operator()(Q const& q) const;

 private:
  template<std::int64_t exponent, typename Q1, typename Q2>
  Quotient<Q2, Exponentiation<Q1, exponent>> Divide(Q2 const& q2) const;

  std::tuple<Qs...> units_;
};

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string Option(std::string const& name,
                   T const& right,
                   OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string Assign(std::string const& name,
                   T const& right,
                   OptionalExpressIn express_in = std::nullopt);

template<typename T, typename U, typename OptionalExpressIn = std::nullopt_t>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::vector<T> const& list,
                          OptionalExpressIn express_in = std::nullopt);

template<typename It, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(It begin, It end,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(bool b,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T,
         typename = std::enable_if_t<std::is_integral_v<T>>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(T integer,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T,
         typename = std::enable_if_t<std::is_floating_point_v<T>>,
         typename OptionalExpressIn = std::nullopt_t,
         typename = void>
std::string ToMathematica(T real,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, int size, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R3Element<T> const& r3_element,
                          OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R3x3Matrix<T> const& r3x3_matrix,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Quaternion const& quaternion,
                          OptionalExpressIn express_in = std::nullopt);

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in);

template<typename D, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Quantity<D> const& quantity,
                          std::nullopt_t express_in = std::nullopt);

template<typename S, typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Vector<S, F> const& vector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename S, typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          OptionalExpressIn express_in = std::nullopt);

template<typename V, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Point<V> const& point,
                          OptionalExpressIn express_in = std::nullopt);

template<typename S,
         typename F,
         template<typename, typename> typename M,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(SymmetricBilinearForm<S, F, M> const& form,
                          OptionalExpressIn express_in = std::nullopt);

template<typename F, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          OptionalExpressIn express_in = std::nullopt);

template<typename Tuple,
         typename = std::enable_if_t<is_tuple_v<Tuple>>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(Tuple const& tuple,
                          OptionalExpressIn express_in = std::nullopt);

template<typename R,
         typename = std::void_t<decltype(std::declval<R>().time)>,
         typename = std::void_t<decltype(std::declval<R>().degrees_of_freedom)>,
         typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(R ref,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    OptionalExpressIn express_in = std::nullopt);

template<typename T, typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::optional<T> const& opt,
                          OptionalExpressIn express_in = std::nullopt);

template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(char const* str,
                          OptionalExpressIn express_in = std::nullopt);
template<typename OptionalExpressIn = std::nullopt_t>
std::string ToMathematica(std::string const& str,
                          OptionalExpressIn express_in = std::nullopt);

// An RAII object to help with Mathematica logging.
class Logger final {
 public:
  // Creates a logger object that will, at destruction, write to the given file.
  // If make_unique is true, a unique id is inserted before the file extension
  // to identify different loggers.
  Logger(std::filesystem::path const& path,
         bool make_unique = true);
  ~Logger();

  // Appends an element to the list of values for the List variable |name|.  The
  // |args...| are passed verbatim to ToMathematica for stringification.  When
  // this object is destroyed, an assignment is generated for each of the
  // variables named in a call to Append.
  template<typename... Args>
  void Append(std::string const& name, Args... args);

  // Sets an element as the single value for the variable |name|.  The
  // |args...| are passed verbatim to ToMathematica for stringification.  When
  // this object is destroyed, an assignment is generated for each of the
  // variables named in a call to Set.
  template<typename... Args>
  void Set(std::string const& name, Args... args);

 private:
  OFStream file_;
  std::map<std::string, std::vector<std::string>> name_and_multiple_values_;
  std::map<std::string, std::string> name_and_single_value_;

  static std::atomic_uint64_t id_;
};

}  // namespace internal_mathematica

using internal_mathematica::Assign;
using internal_mathematica::ExpressIn;
using internal_mathematica::Logger;
using internal_mathematica::Option;
using internal_mathematica::PlottableDataset;
using internal_mathematica::ToMathematica;

}  // namespace mathematica
}  // namespace principia

#include "mathematica/mathematica_body.hpp"
