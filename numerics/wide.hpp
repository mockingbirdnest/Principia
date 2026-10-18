#pragma once

#include <immintrin.h>

#include "base/algebra.hpp"
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

// Do not export this type.
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
  requires (real_dimension<Scalar> == 1)  // NOLINT
class Wide {
 public:
  Wide() = default;

  explicit Wide(Scalar narrow) requires quantity<Scalar>;
  explicit Wide(Scalar narrow)
    requires (!quantity<Scalar>) && quantity<Difference<Scalar>>;  // NOLINT

  static Wide Identity() requires ring<Scalar>;

  explicit operator Scalar() const;

  Wide& operator+=(Wide<Difference<Scalar>> right);
  Wide& operator-=(Wide<Difference<Scalar>> right);

  template<std::integral U>
  Wide& operator*=(U right) requires additive_group<Scalar>;
  template<ring R>
  Wide& operator*=(Wide<R> right) requires module<Scalar, R>;
  template<field K>
  Wide& operator/=(Wide<K> right) requires vector_space<Scalar, K>;

 private:
  explicit Wide(UnderlyingType value);

  UnderlyingType value_;

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

  friend Wide<Product<T, U>> FusedMultiplyAdd(
      Wide<T> a,
      Wide<U> b,
      Wide<Product<T, U>> c);
  template<typename T, typename U>
  friend Wide<Product<T, U>> FusedMultiplySubtract(
      Wide<T> a,
      Wide<U> b,
      Wide<Product<T, U>> c);
  template<typename T, typename U>
  friend Wide<Product<T, U>> FusedNegatedMultiplyAdd(
      Wide<T> a,
      Wide<U> b,
      Wide<Product<T, U>> c);
  template<typename T, typename U>
  friend Wide<Product<T, U>> FusedNegatedMultiplySubtract(
      Wide<T> a,
      Wide<U> b,
      Wide<Product<T, U>> c);

  template<typename T, typename U>
  friend R3Element<Product<T, U>> operator*(R3Element<T> const& left,
                                            Wide<U> right);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> operator*(Wide<T> left,
                                            R3Element<U> const& right);
  template<typename T, typename U>
  friend R3Element<Quotient<T, U>> operator/(R3Element<T> const& left,
                                             Wide<U> right);
  template<typename T, typename U>
  friend R3Element<T>& operator*=(R3Element<T>& left, Wide<U> right)
    requires std::same_as<U, double>;
  template<typename T, typename U>
  friend R3Element<T>& operator/=(R3Element<T>& left, Wide<U> right)
    requires std::same_as<U, double>;

  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedMultiplyAdd(
      R3Element<T> const& a,
      Wide<U> b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedMultiplyAdd(
      Wide<T> a,
      R3Element<U> const& b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedMultiplySubtract(
      R3Element<T> const& a,
      Wide<U> b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedMultiplySubtract(
      Wide<T> a,
      R3Element<U> const& b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
      R3Element<T> const& a,
      Wide<U> b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
      Wide<T> a,
      R3Element<U> const& b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
      R3Element<T> const& a,
      Wide<U> b,
      R3Element<Product<T, U>> const& c);
  template<typename T, typename U>
  friend R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
      Wide<T> a,
      R3Element<U> const& b,
      R3Element<Product<T, U>> const& c);
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

template<typename T, typename U>
Wide<Product<T, U>> FusedMultiplyAdd(Wide<T> a,
                                     Wide<U> b,
                                     Wide<Product<T, U>> c);
template<typename T, typename U>
Wide<Product<T, U>> FusedMultiplySubtract(Wide<T> a,
                                          Wide<U> b,
                                          Wide<Product<T, U>> c);
template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplyAdd(Wide<T> a,
                                            Wide<U> b,
                                            Wide<Product<T, U>> c);
template<typename T, typename U>
Wide<Product<T, U>> FusedNegatedMultiplySubtract(Wide<T> a,
                                                 Wide<U> b,
                                                 Wide<Product<T, U>> c);

// The following are the operations that are optimized by the use of `Wide`.

template<typename T, typename U>
R3Element<Product<T, U>> operator*(R3Element<T> const& left, Wide<U> right);

template<typename T, typename U>
R3Element<Product<T, U>> operator*(Wide<T> left, R3Element<U> const& right);

template<typename T, typename U>
R3Element<Quotient<T, U>> operator/(R3Element<T> const& left, Wide<U> right);

template<typename T, typename U>
R3Element<T>& operator*=(R3Element<T>& left, Wide<U> right)
  requires std::same_as<U, double>;

template<typename T, typename U>
R3Element<T>& operator/=(R3Element<T>& left, Wide<U> right)
  requires std::same_as<U, double>;

template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(R3Element<T> const& a,
                                          Wide<U> b,
                                          R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplyAdd(Wide<T> a,
                                          R3Element<U> const& b,
                                          R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedMultiplySubtract(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplyAdd(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    R3Element<T> const& a,
    Wide<U> b,
    R3Element<Product<T, U>> const& c);
template<typename T, typename U>
R3Element<Product<T, U>> FusedNegatedMultiplySubtract(
    Wide<T> a,
    R3Element<U> const& b,
    R3Element<Product<T, U>> const& c);

// The following don't need to be friend, they access the `coordinates` member.

template<r3_element_coordinates T, typename U>
Product<T, U> operator*(T const& left, Wide<U> right);

template<typename T, r3_element_coordinates U>
Product<T, U> operator*(Wide<T> left, U const& right);

template<r3_element_coordinates T, typename U>
Quotient<T, U> operator/(T const& left, Wide<U> right);

template<r3_element_coordinates T>
T& operator*=(T& left, Wide<double> right);

template<r3_element_coordinates T>
T& operator/=(T& left, Wide<double> right);

template<r3_element_coordinates T, typename U>
Product<T, U> FusedMultiplyAdd(T const& a,
                               Wide<U> b,
                               Product<T, U> const& c);
template<typename T, r3_element_coordinates U>
Product<T, U> FusedMultiplyAdd(Wide<T> a,
                               U const& b,
                               Product<T, U> const& c);
template<r3_element_coordinates T, typename U>
Product<T, U> FusedMultiplySubtract(T const& a,
                                    Wide<U> b,
                                    Product<T, U> const& c);
template<typename T, r3_element_coordinates U>
Product<T, U> FusedMultiplySubtract(Wide<T> a,
                                    U const& b,
                                    Product<T, U> const& c);
template<r3_element_coordinates T, typename U>
Product<T, U> FusedNegatedMultiplyAdd(T const& a,
                                      Wide<U> b,
                                      Product<T, U> const& c);
template<typename T, r3_element_coordinates U>
Product<T, U> FusedNegatedMultiplyAdd(Wide<T> a,
                                      U const& b,
                                      Product<T, U> const& c);
template<r3_element_coordinates T, typename U>
Product<T, U> FusedNegatedMultiplySubtract(T const& a,
                                           Wide<U> b,
                                           Product<T, U> const& c);
template<typename T, r3_element_coordinates U>
Product<T, U> FusedNegatedMultiplySubtract(Wide<T> a,
                                           U const& b,
                                           Product<T, U> const& c);

}  // namespace internal

using internal::Wide;

}  // namespace _wide
}  // namespace numerics
}  // namespace principia

#include "numerics/wide_body.hpp"
