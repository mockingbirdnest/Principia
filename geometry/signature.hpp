#pragma once

#include "geometry/frame.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

template<typename FromFrame, typename ToFrame>
class Signature {
 public:
  static constexpr Sign determinant =
      FromFrame::handedness == ToFrame::handedness ? Sign::Positive()
                                                   : Sign::Negative();
  struct deduce_sign_t {};
  static constexpr deduce_sign_t deduce_sign_from_handedness{};

  constexpr Signature(Sign x, Sign y, deduce_sign_t z);
  constexpr Signature(Sign x, deduce_sign_t y, Sign z);
  constexpr Signature(deduce_sign_t x, Sign y, Sign z);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness == T::handedness>>
  static constexpr Signature Identity();

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness != T::handedness>>
  static constexpr Signature CentralReflection();

  constexpr Signature<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename Scalar>
  SymmetricBilinearForm<Scalar, ToFrame> operator()(
      SymmetricBilinearForm<Scalar, FromFrame> const& form) const;

  // TODO(egg): Serialization.

 private:
  constexpr Signature(Sign x, Sign y, Sign z);

  Sign x_;
  Sign y_;
  Sign z_;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Signature<FromFrame, ToFrame> operator*(
    Signature<ThroughFrame, ToFrame> const& left,
    Signature<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Signature<FromFrame, ToFrame> const& signature);

}  // namespace internal_signature

using internal_signature::Signature;

}  // namespace geometry
}  // namespace principia

#include "geometry/signature_body.hpp"
