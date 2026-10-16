#pragma once

#include <immintrin.h>

#include "base/algebra.hpp"
#include "base/macros.hpp"
#include "base/traits.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/concepts.hpp"

namespace principia {
namespace numerics {
namespace _wide {
namespace internal {

using namespace base::_algebra;
using namespace base::_traits;
using namespace geometry::_grassmann;
using namespace geometry::_r3_element;
using namespace quantities::_concepts;

#if PRINCIPIA_USE_AVX()
using UnderlyingType = __m256d;
#else
using UnderlyingType = __m128d;
#endif

template<typename T>
concept r3_element_coordinates = requires(T v) {
  requires instance<std::remove_cvref_t<decltype(v.coordinates())>, R3Element>;
};

// For a type `Scalar` represented by a `double`, a `Wide<Scalar>` behaves like
// a `Scalar`, but is broadcast to all elements of a register of the appropriate
// width for fast multiplications with `R3Element` or types wrapping it; this
// can be used to factor out broadcast instructions in performance-critical
// code.  Since the point of using `Wide` is to minimize the cost incurred by
// the implicit widening in scalar-vector multiplications widening, `Wide` is
// not implicitly constructible from the narrow scalar, and operations between
// wide scalars require that both arguments be `Wide`: `2.0 * Wide(2.0)` is
// ill-formed, and must be written `Wide(2.0) * Wide(2.0)` instead.
template<real_affine_space Scalar>
  requires(real_dimension<Scalar> == 1)
class Wide {
 public:
  Wide() = default;

  explicit Wide(Scalar narrow);

  static Wide Identity()
    requires ring<Scalar>;

  Scalar narrow() const;

  template<additive_group T>
  Wide<T>& operator+=(Wide<T>& left, Wide<T> right);
  template<additive_group T>
  Wide<T>& operator-=(Wide<T>& left, Wide<T> right);

  template<typename T, std::integral U>
  Wide<T>& operator*=(Wide<T>& left, U right);
  template<ring T>
  Wide<T>& operator*=(Wide<T>& left, Wide<T> right);

 private:
  explicit Wide(UnderlyingType value);

  UnderlyingType value_;

  template<quantity T>
  friend Wide<T> Widen(T narrow);

  template<additive_group T>
  friend Wide<T> operator-(Wide<T> right);
  template<affine T, affine U>
  friend Wide<Sum<T, U>> operator+(Wide<T> left, Wide<U> right);
  template<affine T, affine U>
  friend Wide<Difference<T, U>> operator-(Wide<T> left, Wide<U> right);
  template<typename T, typename U>
  friend Wide<Product<T, U>> operator*(Wide<T> left, Wide<U> right);
  template<typename T, typename U>
  friend Wide<Quotient<T, U>> operator/(Wide<T> left, Wide<U> right);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> operator*(R3Element<T> const& left,
                                            Wide<U> right);
  template<typename T, instance<Multivector> U>
  friend R3Element<Product<T, U>> operator*(Wide<T> left,
                                            R3Element<U> const& right);
  template<typename T, typename U>
  friend R3Element<Quotient<T, U>> operator/(R3Element<T> const& left,
                                             Wide<U> right);
  template<typename T>
  friend R3Element<T>& operator*=(R3Element<T>& left, Wide right)
    requires std::same_as<Scalar, double>;
  template<typename T>
  friend R3Element<T>& operator/=(R3Element<T>& left, Wide right)
    requires std::same_as<Scalar, double>;
};

template<additive_group T>
Wide<T> operator+(Wide<T> right);

template<additive_group T>
Wide<T> operator-(Wide<T> right);

template<affine T, affine U>
Wide<Sum<T, U>> operator+(Wide<T> left, Wide<U> right);

template<affine T, affine U>
Wide<Difference<T, U>> operator-(Wide<T> left, Wide<U> right);

template<typename T, typename U>
Wide<Product<T, U>> operator*(Wide<T> left, Wide<U> right);

template<typename T, typename U>
Wide<Quotient<T, U>> operator/(Wide<T> left, Wide<U> right);

template<typename T, std::integral U>
Wide<T> operator*(Wide<T> left, U right);

template<typename T, std::integral U>
Wide<T> operator*(U left, Wide<T> right);

// The following are the operations that are optimized by the use of `Wide`.

template<r3_element_coordinates T, typename U>
Product<T, U> operator*(T const& left, Wide<U> right);

template<typename T, r3_element_coordinates U>
Product<T, U> operator*(Wide<T> left, U const& right);

template<r3_element_coordinates T, typename U>
Quotient<T, U> operator/(T left, Wide<U> right);

template<r3_element_coordinates T>
T& operator*=(T& left, Wide<double> right);

template<r3_element_coordinates T>
T& operator/=(T& left, Wide<double> right);

template<typename T, typename U>
R3Element<Product<T, U>> operator*(R3Element<T> const& left, Wide<U> right);

template<typename T, instance<Multivector> U>
R3Element<Product<T, U>> operator*(Wide<T> left, R3Element<U> const& right);

template<typename T, typename U>
R3Element<Quotient<T, U>> operator/(R3Element<T> const& left, Wide<U> right);

template<typename T>
R3Element<T>& operator*=(R3Element<T>& left, Wide<double> right);

template<typename T>
R3Element<T>& operator/=(R3Element<T>& left, Wide<double> right);

}  // namespace internal

using internal::Wide;

}  // namespace _wide
}  // namespace numerics
}  // namespace principia

#include "numerics/wide_body.hpp"
