#pragma once

#include "physics/jacobi_coordinates.hpp"

#include <experimental/optional>
#include <iterator>

namespace principia {

namespace physics {

template<typename Frame>
JacobiCoordinates<Frame>::JacobiCoordinates(MassiveBody const& primary) {
  DegreesOfFreedom<PrimocentricFrame> motionless_origin = {
      PrimocentricFrame::origin, Velocity<PrimocentricFrame>()};
  primocentric_dof_.emplace_back(motionless_origin);
  system_barycentre_.Add(primocentric_dof_.back(),
                         primary.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    RelativeDegreesOfFreedom<Frame> const& dof_wrt_system) {
  primocentric_dof_.emplace_back(system_barycentre_.Get() +
                                 id_fp_(dof_wrt_system));
  system_barycentre_.Add(primocentric_dof_.back(),
                         body.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    KeplerianElements<Frame> const& osculating_elements_wrt_system) {
  Instant const epoch;
  Add(body,
      KeplerOrbit<Frame>(/*primary=*/System(),
                         /*secondary=*/body,
                         osculating_elements_wrt_system,
                         epoch).StateVectors(epoch));
}

template<typename Frame>
MassiveBody JacobiCoordinates<Frame>::System() const {
  return MassiveBody(system_barycentre_.weight());
}

template<typename Frame>
std::vector<RelativeDegreesOfFreedom<Frame>>
JacobiCoordinates<Frame>::BarycentricCoordinates() const {
  std::vector<RelativeDegreesOfFreedom<Frame>> result;
  DegreesOfFreedom<PrimocentricFrame> system_barycentre =
      system_barycentre_.Get();
  for (auto const& dof : primocentric_dof_) {
    result.emplace_back(id_pf_(dof - system_barycentre));
  }
  return result;
}

template<typename Frame>
Identity<typename JacobiCoordinates<Frame>::PrimocentricFrame, Frame> const
    JacobiCoordinates<Frame>::id_pf_;
template<typename Frame>
Identity<Frame, typename JacobiCoordinates<Frame>::PrimocentricFrame> const
    JacobiCoordinates<Frame>::id_fp_;

template<typename Frame>
HierarchicalSystem<Frame>::HierarchicalSystem(
    not_null<std::unique_ptr<MassiveBody const>> primary)
    : system_(std::move(primary)) {
  subsystems_[system_.primary.get()] = &system_;
}

template<typename Frame>
void HierarchicalSystem<Frame>::Add(
    not_null<std::unique_ptr<MassiveBody const>> body,
    not_null<MassiveBody const*> const parent,
    KeplerianElements<Frame> const& jacobi_osculating_elements) {
  System& parent_system = *subsystems_[parent];
  parent_system.satellites.emplace_back(
      make_not_null_unique<Subsystem>(std::move(body)));
  not_null<Subsystem*> inserted_system = parent_system.satellites.back().get();
  {
    not_null<MassiveBody const*> body = inserted_system->primary.get();
    subsystems_[body] = inserted_system;
    inserted_system->jacobi_osculating_elements = jacobi_osculating_elements;
  }
}

template<typename Frame>
typename HierarchicalSystem<Frame>::BarycentricSystem
HierarchicalSystem<Frame>::Get() {
  auto const semimajor_axis_less_than = [](
      not_null<std::unique_ptr<Subsystem>> const& left,
      not_null<std::unique_ptr<Subsystem>> const& right) -> bool {
    return left->jacobi_osculating_elements.semimajor_axis <
           right->jacobi_osculating_elements.semimajor_axis;
  };

  struct BarycentricSubystem {
    std::unique_ptr<MassiveBody> equivalent_body;
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    std::vector<RelativeDegreesOfFreedom<Frame>> barycentric_degrees_of_freedom;
  };

  std::function<BarycentricSubystem(System&)> to_barycentric =
      [&semimajor_axis_less_than, &to_barycentric](System& system) {
        std::sort(system.satellites.begin(),
                  system.satellites.end(),
                  semimajor_axis_less_than);

        BarycentricSubystem result;

        // Jacobi coordinates for |system|, with satellite subsystems treated
        // as point masses at their barycentres.
        JacobiCoordinates<Frame> jacobi_coordinates(*system.primary);
        result.bodies.emplace_back(std::move(system.primary));

        // The |n|th element of |satellite_degrees_of_freedom| contains the list
        // of degrees of freedom of the bodies in the |n|th satellite subsystem
        // with respect to their barycentre.
        std::vector<std::vector<RelativeDegreesOfFreedom<Frame>>>
            satellite_degrees_of_freedom;

        // Fill |satellite_degrees_of_freedom|, |jacobi_coordinates|, and
        // |result.bodies|.
        for (auto const& subsystem : system.satellites) {
          BarycentricSubystem barycentric_satellite_subsystem =
              to_barycentric(*subsystem);
          satellite_degrees_of_freedom.emplace_back(std::move(
              barycentric_satellite_subsystem.barycentric_degrees_of_freedom));
          jacobi_coordinates.Add(
              *barycentric_satellite_subsystem.equivalent_body,
              subsystem->jacobi_osculating_elements);
          std::move(barycentric_satellite_subsystem.bodies.begin(),
                    barycentric_satellite_subsystem.bodies.end(),
                    std::back_inserter(result.bodies));
        }

        // Fill |result.barycentric_degrees_of_freedom|.
        std::vector<RelativeDegreesOfFreedom<Frame>> const
            barycentres_of_subsystems =
                jacobi_coordinates.BarycentricCoordinates();
        // The primary.
        result.barycentric_degrees_of_freedom.emplace_back(
            barycentres_of_subsystems.front());
        for (int n = 0; n < satellite_degrees_of_freedom.size(); ++n) {
          // |n + 1| because the primary is at |barycentres_of_subsystems[0]|.
          RelativeDegreesOfFreedom<Frame> const subsystem_barycentre =
              barycentres_of_subsystems[n + 1];
          for (auto const& body_dof_wrt_subsystem_barycentre :
               satellite_degrees_of_freedom[n]) {
            result.barycentric_degrees_of_freedom.emplace_back(
                body_dof_wrt_subsystem_barycentre + subsystem_barycentre);
          }
        }

        result.equivalent_body =
            std::make_unique<MassiveBody>(jacobi_coordinates.System());
        return std::move(result);
      };

  to_barycentric(system_);


  return BarycentricSystem();
}

}  // namespace physics
}  // namespace principia
