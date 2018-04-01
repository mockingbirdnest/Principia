
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
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_oblate_body {

using base::not_null;
using geometry::Vector;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Order2ZonalCoefficient;
using quantities::Quotient;

template<typename Frame>
class OblateBody : public RotatingBody<Frame> {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  class PHYSICS_DLL Parameters final {
   public:
    explicit Parameters(Order2ZonalCoefficient const& j2);
    Parameters(double const j2,
               Length const& reference_radius);

   private:
    std::optional<Order2ZonalCoefficient> j2_;
    std::optional<
        Quotient<Order2ZonalCoefficient, GravitationalParameter>> j2_over_μ_;
    template<typename F>
    friend class OblateBody;
  };

  OblateBody(MassiveBody::Parameters const& massive_body_parameters,
             typename RotatingBody<Frame>::Parameters const&
                 rotating_body_parameters,
             Parameters const& parameters);

  // Returns the j2 coefficient.
  Order2ZonalCoefficient const& j2() const;

  // Returns |j2 / μ|.
  Quotient<Order2ZonalCoefficient,
           GravitationalParameter> const& j2_over_μ() const;

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
