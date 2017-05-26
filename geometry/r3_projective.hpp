
#pragma once

#include "geometry/r3_element.hpp"

namespace principia {
namespace geometry {

template<typename Scalar>
class R3Projective {
 public:
  explicit R3Projective(R3Element<Scalar> const& coordinates);

  bool at_infinity() const;
  Scalar slope() const;

  Scalar const x() const;
  Scalar const y() const;

 private:
  R3Element<Scalar> coordinates_;

  friend bool operator==(R3Projective<Scalar> const& left,
                         R3Projective<Scalar> const& right);
  friend bool operator!=(R3Projective<Scalar> const& left,
                         R3Projective<Scalar> const& right);
};

template<typename Scalar>
bool operator==(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/projective_body.hpp"
