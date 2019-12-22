#pragma once

#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

using base::not_constructible;

struct DeduceSignPreservingOrientation final {};
struct DeduceSignReversingOrientation final {};

template<typename FromFrame, typename ToFrame>
class Signature {
 public:
  static constexpr Sign determinant =
      FromFrame::handedness == ToFrame::handedness ? Sign::Positive()
                                                   : Sign::Negative();
  using DeduceSign = std::conditional_t<determinant.is_positive(),
                                        DeduceSignPreservingOrientation,
                                        DeduceSignReversingOrientation>;

  constexpr Signature(Sign x, Sign y, DeduceSign z);
  constexpr Signature(Sign x, DeduceSign y, Sign z);
  constexpr Signature(DeduceSign x, Sign y, Sign z);

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

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

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

  template<typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Signature<From, To> const& signature);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Signature<FromFrame, ToFrame> operator*(
    Signature<ThroughFrame, ToFrame> const& left,
    Signature<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Signature<FromFrame, ToFrame> const& signature);

}  // namespace internal_signature

using internal_signature::DeduceSignPreservingOrientation;
using internal_signature::DeduceSignReversingOrientation;
using internal_signature::Signature;

}  // namespace geometry
}  // namespace principia

#include "geometry/signature_body.hpp"
