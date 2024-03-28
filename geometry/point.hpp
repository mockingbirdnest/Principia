#pragma once

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _point {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;

// Point<Vector> is an affine space on the vector space Vector.
template<typename Vector>
class Point final {
  // This cannot be a constraint, as it would lead to recursive instantation of
  // Position, which is used in Frame.
  static_assert(additive_group<Vector>);

 public:
  constexpr Point();

#if PRINCIPIA_COMPILER_MSVC && !__INTELLISENSE__
  // Explicitly define constexpr default copy and move constructors because
  // otherwise MSVC fails to initialize constant expressions.  In addition,
  // Intellisense gets confused by these (because of course MSVC and
  // Intellisense are different compilers and have bugs in different places).
  constexpr Point(Point const& other);
  constexpr Point(Point&& other);
  Point& operator=(Point const& other) = default;
  Point& operator=(Point&& other) = default;
#endif

  constexpr friend bool operator==(Point const& left,
                                   Point const& right) = default;
  constexpr friend bool operator!=(Point const& left,
                                   Point const& right) = default;
  constexpr friend auto operator<=>(Point const& left, Point const& right)
    requires convertible_to_quantity<Vector> = default;

  constexpr Vector operator-(Point const& from) const;

  constexpr Point operator+(Vector const& translation) const;
  constexpr Point operator-(Vector const& translation) const;

  Point& operator+=(Vector const& translation);
  Point& operator-=(Vector const& translation);

  void WriteToMessage(not_null<serialization::Point*> message) const;
  static Point ReadFromMessage(serialization::Point const& message)
    requires serializable<Vector>;

 private:
  // This constructor allows for C++11 functional constexpr operators, and
  // possibly move-magic.
  // TODO(egg): We may want to reconsider this after we truly have C++14.
  constexpr explicit Point(Vector const& coordinates);

  Vector coordinates_;

  template<typename V>
  friend constexpr Point<V> operator+(V const& translation,
                                      Point<V> const& point);
  template<typename L, typename R>
  friend Point<Product<L, R>> FusedMultiplyAdd(L const& a, R const& b,
                                               Point<Product<L, R>> const& c);
  template<typename L, typename R>
  friend Point<Product<L, R>> FusedNegatedMultiplyAdd(
      L const& a, R const& b, Point<Product<L, R>> const& c);

  template<typename V>
    requires convertible_to_quantity<V>
  friend constexpr Point<V> NextUp(Point<V> x);
  template<typename V>
    requires convertible_to_quantity<V>
  friend constexpr Point<V> NextDown(Point<V> x);

  template<typename V>
  friend std::string DebugString(Point<V> const& point);
};

template<typename Vector>
constexpr Point<Vector> operator+(Vector const& translation,
                                  Point<Vector> const& point);

template<typename L, typename R>
Point<Product<L, R>> FusedMultiplyAdd(L const& a,
                                      R const& b,
                                      Point<Product<L, R>> const& c);
template<typename L, typename R>
Point<Product<L, R>> FusedNegatedMultiplyAdd(L const& a,
                                             R const& b,
                                             Point<Product<L, R>> const& c);

template<typename Vector>
  requires convertible_to_quantity<Vector>
constexpr Point<Vector> NextUp(Point<Vector> x);
template<typename Vector>
  requires convertible_to_quantity<Vector>
constexpr Point<Vector> NextDown(Point<Vector> x);

template<typename Vector>
std::string DebugString(Point<Vector> const& point);

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Point<Vector> const& point);

}  // namespace internal

using internal::Point;

}  // namespace _point
}  // namespace geometry
}  // namespace principia

#include "geometry/point_body.hpp"
