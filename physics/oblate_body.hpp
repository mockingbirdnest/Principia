
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#include "physics/rotating_body.hpp"
#endif  // PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#ifndef PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
#define PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_

#include <vector>

#include "geometry/grassmann.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_oblate_body {

using base::not_null;
using geometry::Vector;
using numerics::FixedLowerTriangularMatrix;
using quantities::Degree2SphericalHarmonicCoefficient;
using quantities::Degree3SphericalHarmonicCoefficient;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Quotient;

template<typename Frame>
class OblateBody : public RotatingBody<Frame> {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  static constexpr int max_geopotential_degree = 50;
  using GeopotentialCoefficients =
      FixedLowerTriangularMatrix<double, max_geopotential_degree + 1>;

  class Parameters final {
   public:
    Parameters(double j2,
               Length const& reference_radius);

    static Parameters ReadFromMessage(
        serialization::OblateBody::Geopotential const& message,
        Length const& reference_radius);

    void WriteToMessage(
        not_null<serialization::OblateBody::Geopotential*> message) const;

   private:
    // Only for use when building from a geopotential.
    explicit Parameters(Length const& reference_radius);

    Length reference_radius_;

    double j2_;
    Quotient<Degree2SphericalHarmonicCoefficient, GravitationalParameter>
        j2_over_μ_;
    GeopotentialCoefficients cos_;
    GeopotentialCoefficients sin_;
    int degree_;
    bool is_zonal_;

    template<typename F>
    friend class OblateBody;
  };

  OblateBody(MassiveBody::Parameters const& massive_body_parameters,
             typename RotatingBody<Frame>::Parameters const&
                 rotating_body_parameters,
             Parameters const& parameters);

  // These parameters are unnormalized.
  double j2() const;
  Quotient<Degree2SphericalHarmonicCoefficient,
           GravitationalParameter> const& j2_over_μ() const;

  // These parameters are normalized.
  GeopotentialCoefficients const& cos() const;
  GeopotentialCoefficients const& sin() const;
  int geopotential_degree() const;

  // Returns true iff the geopotential only contains zonal terms.
  bool is_zonal() const;

  Length const& reference_radius() const;

  // Returns false.
  bool is_massless() const override;

  // Returns true.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const override;

  static not_null<std::unique_ptr<OblateBody<Frame>>> ReadFromMessage(
      serialization::OblateBody const& message,
      MassiveBody::Parameters const& massive_body_parameters,
      typename RotatingBody<Frame>::Parameters const& rotating_body_parameters);

 private:
  Parameters parameters_;
};

}  // namespace internal_oblate_body

using internal_oblate_body::OblateBody;

}  // namespace physics
}  // namespace principia

#include "physics/oblate_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
