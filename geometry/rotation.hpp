
#pragma once

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class OrthogonalMap;

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation);

// An orientation-preserving orthogonal map between the inner product spaces
// |FromFrame| and |ToFrame|, as well as the induced maps on the exterior
// algebra.
template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
 public:
  explicit Rotation(Quaternion const& quaternion);
  template<typename Scalar>
  Rotation(quantities::Angle const& angle,
           Bivector<Scalar, FromFrame> const& axis);

  // Construct a rotation from the axes of |ToFrame|, expressed in |FromFrame|.
  template<int rank_x, int rank_y, int rank_z,
           typename = std::enable_if_t<
                          (rank_x + rank_y + rank_z) % 2 == 0 &&
                          !std::is_same<FromFrame, ToFrame>::value>>
  Rotation(Multivector<double, FromFrame, rank_x> x_to_frame,
           Multivector<double, FromFrame, rank_y> y_to_frame,
           Multivector<double, FromFrame, rank_z> z_to_frame);
  // Construct a rotation from the axes of |FromFrame|, expressed in |ToFrame|.
  // The |typename = void| template parameter gives it a different signature
  // from the above when |FromFrame| and |ToFrame| are the same (otherwise this
  // is a problem when instantiating the class, even though they are both
  // disabled by |enable_if|).
  template<int rank_x, int rank_y, int rank_z,
           typename = std::enable_if_t<
                          (rank_x + rank_y + rank_z) % 2 == 0 &&
                          !std::is_same<FromFrame, ToFrame>::value>,
           typename = void>
  Rotation(Multivector<double, ToFrame, rank_x> x_from_frame,
           Multivector<double, ToFrame, rank_y> y_from_frame,
           Multivector<double, ToFrame, rank_z> z_from_frame);

  ~Rotation() override = default;

  Sign Determinant() const override;

  Rotation<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename T>
  typename base::Mappable<Rotation, T>::type operator()(T const& t) const;

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

  static Rotation Identity();

  Quaternion const& quaternion() const;

  void WriteToMessage(not_null<serialization::LinearMap*> const message) const;
  static Rotation ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Rotation*> const message) const;
  static Rotation ReadFromMessage(serialization::Rotation const& message);

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  Quaternion quaternion_;

  // For constructing a rotation using a quaternion.
  template<typename From, typename To>
  friend class Permutation;

  template<typename From, typename Through, typename To>
  friend Rotation<From, To> operator*(Rotation<Through, To> const& left,
                                      Rotation<From, Through> const& right);

  friend std::ostream& operator<<<>(std::ostream& out,  // NOLINT
                                    Rotation const& rotation);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/rotation_body.hpp"
