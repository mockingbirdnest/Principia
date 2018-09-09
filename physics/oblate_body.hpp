
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#include "physics/rotating_body.hpp"
#endif  // PRINCIPIA_PHYSICS_ROTATING_BODY_HPP_
#ifndef PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
#define PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_

#include <optional>
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

// TODO(phl): c22, s22 and j3 in this class are temporary until we get the
// general geopotential model to work at least as well as c22/s22/j3.
template<typename Frame>
class OblateBody : public RotatingBody<Frame> {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  static constexpr int max_geopotential_degree = 5;
  using GeopotentialCoefficients =
      FixedLowerTriangularMatrix<double, max_geopotential_degree + 1>;

  class PHYSICS_DLL Parameters final {
   public:
    Parameters(double j2,
               Length const& reference_radius);

    Parameters(double j2,
               double c22,
               double s22,
               Length const& reference_radius);

    Parameters(double j2,
               double c22,
               double s22,
               double j3,
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

    std::optional<double> j2_;
    std::optional<Degree2SphericalHarmonicCoefficient> pre_descartes_j2_;
    std::optional<Quotient<Degree2SphericalHarmonicCoefficient,
                           GravitationalParameter>> j2_over_μ_;
    std::optional<Degree2SphericalHarmonicCoefficient> c22_;
    std::optional<Quotient<Degree2SphericalHarmonicCoefficient,
                           GravitationalParameter>> c22_over_μ_;
    std::optional<Degree2SphericalHarmonicCoefficient> s22_;
    std::optional<Quotient<Degree2SphericalHarmonicCoefficient,
                           GravitationalParameter>> s22_over_μ_;
    std::optional<Degree3SphericalHarmonicCoefficient> j3_;
    std::optional<Quotient<Degree3SphericalHarmonicCoefficient,
                           GravitationalParameter>> j3_over_μ_;
    std::optional<GeopotentialCoefficients> cos_;
    std::optional<GeopotentialCoefficients> sin_;
    std::optional<int> degree_;

    template<typename F>
    friend class OblateBody;
  };

  OblateBody(MassiveBody::Parameters const& massive_body_parameters,
             typename RotatingBody<Frame>::Parameters const&
                 rotating_body_parameters,
             Parameters const& parameters);

  // Selectors for the various spherical harmonics coefficients.

  // These parameters are unnormalized.
  double j2() const;
  Quotient<Degree2SphericalHarmonicCoefficient,
           GravitationalParameter> const& j2_over_μ() const;
  Degree2SphericalHarmonicCoefficient const c22() const;
  Quotient<Degree2SphericalHarmonicCoefficient,
           GravitationalParameter> const c22_over_μ() const;
  Degree2SphericalHarmonicCoefficient const s22() const;
  Quotient<Degree2SphericalHarmonicCoefficient,
           GravitationalParameter> const s22_over_μ() const;
  Degree3SphericalHarmonicCoefficient const j3() const;
  Quotient<Degree3SphericalHarmonicCoefficient,
           GravitationalParameter> const j3_over_μ() const;

  // These parameters are normalized.
  GeopotentialCoefficients const& cos() const;
  GeopotentialCoefficients const& sin() const;

  Length const& reference_radius() const;

  // Whether this body has a c22, s22, j3 or geopotential.
  bool has_c22() const;
  bool has_s22() const;
  bool has_j3() const;
  bool has_geopotential() const;

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

#if !PHYSICS_DLL_IMPORT
#include "physics/oblate_body_body.hpp"
#endif

#endif  // PRINCIPIA_PHYSICS_OBLATE_BODY_HPP_
