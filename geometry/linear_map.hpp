
#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_linear_map {

using base::not_null;

template<typename FromFrame, typename ToFrame>
class LinearMap {
 public:
  LinearMap() = default;
  virtual ~LinearMap() = default;

  virtual Sign Determinant() const = 0;

// The following is the contract that must be implemented by subclasses.
// Apologies for the commented-out code, but we cannot write this in real C++
// because templates cannot be virtual and because the return type is not
// covariant in inheritance.
//
//   virtual LinearMap<ToFrame, FromFrame> Inverse() const = 0;
//
//   template<typename Scalar>
//   virtual Vector<Scalar, ToFrame> operator()(
//       Vector<Scalar, FromFrame> const& vector) const = 0;
//
//   template<typename Scalar>
//   virtual Bivector<Scalar, ToFrame> operator()(
//       Bivector<Scalar, FromFrame> const& bivector) const = 0;
//
//   template<typename Scalar>
//   virtual Trivector<Scalar, ToFrame> operator()(
//       Trivector<Scalar, FromFrame> const& trivector) const = 0;
//
//   template<typename Scalar>
//   SymmetricBilinearForm<Scalar, ToFrame> operator()(
//       SymmetricBilinearForm<Scalar, FromFrame> const& form) const = 0;
//
//   template<typename T>
//   typename base::Mappable<LinearMap, T>::type operator()(T const& t) const;
//
 protected:
  // Serialization of the frames.  These are just helper functions for
  // implementing the subclasses, they don't dispatch to the subclasses.
  static void WriteToMessage(not_null<serialization::LinearMap*> message);
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static void ReadFromMessage(serialization::LinearMap const& message);

//   template<typename Scalar>
//   virtual R3Element<Scalar> operator()(
//       R3Element<Scalar> const& r3_element) const = 0;
};

}  // namespace internal_linear_map

using internal_linear_map::LinearMap;

}  // namespace geometry
}  // namespace principia

#include "geometry/linear_map_body.hpp"
