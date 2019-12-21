#pragma once

#include "geometry/frame.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

using base::not_constructible;

struct deduce_sign_direct_t final {};
struct deduce_sign_indirect_t final {};

constexpr deduce_sign_direct_t deduce_sign_for_direct_isometry{};
constexpr deduce_sign_indirect_t deduce_sign_for_indirect_isometry{};

template<typename FromFrame, typename ToFrame>
class Signature {
 public:
  static constexpr Sign determinant =
      FromFrame::handedness == ToFrame::handedness ? Sign::Positive()
                                                   : Sign::Negative();
  using deduce_sign_t = std::conditional_t<determinant.Positive(),
                                           deduce_sign_direct_t,
                                           deduce_sign_indirect_t>;

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
  static constexpr Signature CentralInversion();

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

  template<typename From, typename To>
  friend class Signature;

  template<typename From, typename Through, typename To>
  friend Signature<From, To> operator*(Signature<Through, To> const& left,
                                       Signature<From, Through> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Signature<FromFrame, ToFrame> operator*(
    Signature<ThroughFrame, ToFrame> const& left,
    Signature<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Signature<FromFrame, ToFrame> const& signature);

}  // namespace internal_signature

using internal_signature::deduce_sign_from_handedness;
using internal_signature::Signature;

}  // namespace geometry
}  // namespace principia

#include "geometry/signature_body.hpp"
