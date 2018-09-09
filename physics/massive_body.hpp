
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#endif  // PRINCIPIA_PHYSICS_BODY_HPP_
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_

#include <string>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_massive_body {

using base::not_null;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Mass;

class PHYSICS_DLL MassiveBody : public Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G.  The parameter must not
  // be zero.
  class PHYSICS_DLL Parameters final {
   public:
    // The constructors are implicit on purpose.
    Parameters(
        GravitationalParameter const& gravitational_parameter);  // NOLINT
    Parameters(std::string const& name,
               GravitationalParameter const& gravitational_parameter);
    Parameters(Mass const& mass);  // NOLINT(runtime/explicit)
    Parameters(std::string const& name,
               Mass const& mass);

    // For pre-Διόφαντος compatibility.
    GravitationalParameter const& gravitational_parameter() const;

   private:
    std::string const name_;
    GravitationalParameter const gravitational_parameter_;
    Mass const mass_;
    friend class MassiveBody;
  };

  explicit MassiveBody(Parameters const& parameters);

  // Returns the construction parameter.
  std::string const& name() const;
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns zero.
  virtual Length mean_radius() const;

  // Returns false.
  bool is_massless() const override;

  // Returns false.
  bool is_oblate() const override;

  // Call the following |WriteToMessage|, which dispatches to the proper
  // subclass.
  void WriteToMessage(not_null<serialization::Body*> message) const override;

  // Must be overridden by each subclass and first call the same method of the
  // superclass.
  virtual void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const;

  // |message.has_massive_body()| must be true.
  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::Body const& message);

  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::MassiveBody const& message);

 private:
  Parameters const parameters_;
};

}  // namespace internal_massive_body

using internal_massive_body::MassiveBody;

}  // namespace physics
}  // namespace principia

#if !PHYSICS_DLL_IMPORT
#include "physics/massive_body_body.hpp"
#endif

#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
