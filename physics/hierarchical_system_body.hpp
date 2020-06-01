
#pragma once

#include "physics/hierarchical_system.hpp"

#include <algorithm>
#include <iterator>
#include <vector>

#include "google/protobuf/repeated_field.h"

namespace principia {
namespace physics {
namespace internal_hierarchical_system {

using base::make_not_null_unique;
using geometry::Identity;
using geometry::Velocity;

template<typename Frame>
HierarchicalSystem<Frame>::HierarchicalSystem(
    not_null<std::unique_ptr<MassiveBody const>> primary)
    : system_(std::move(primary)) {
  systems_[system_.primary.get()] = &system_;
}

template<typename Frame>
void HierarchicalSystem<Frame>::Add(
    not_null<std::unique_ptr<MassiveBody const>> body,
    not_null<MassiveBody const*> const parent,
    KeplerianElements<Frame> const& jacobi_osculating_elements) {
  not_null<MassiveBody const*> unowned_body = body.get();
  System& parent_system = *systems_[parent];
  parent_system.satellites.emplace_back(
      make_not_null_unique<Subsystem>(std::move(body)));
  {  // Hide the moved-from |body|.
    not_null<MassiveBody const*> body = unowned_body;
    not_null<Subsystem*> inserted_system =
        parent_system.satellites.back().get();
    systems_[body] = inserted_system;
    inserted_system->jacobi_osculating_elements = jacobi_osculating_elements;
  }
}

template<typename Frame>
typename HierarchicalSystem<Frame>::BarycentricSystem
HierarchicalSystem<Frame>::ConsumeBarycentricSystem() {
  BarycentricSystem result;
  auto barycentric_result = ToBarycentric(system_);
  result.bodies = std::move(barycentric_result.bodies);
  static DegreesOfFreedom<Frame> const system_barycentre = {Frame::origin,
                                                            Frame::unmoving};
  for (auto const& barycentric_dof :
       barycentric_result.barycentric_degrees_of_freedom) {
    result.degrees_of_freedom.emplace_back(system_barycentre + barycentric_dof);
  }
  return std::move(result);
}

template<typename Frame>
void HierarchicalSystem<Frame>::WriteToMessage(
    not_null<serialization::HierarchicalSystem*> const message) const {
  system_.primary->WriteToMessage(message->mutable_system()->mutable_primary());
  WriteToMessage(system_.satellites,
                 message->mutable_system()->mutable_satellites());
}

template<typename Frame>
typename HierarchicalSystem<Frame>::BarycentricSubsystem
HierarchicalSystem<Frame>::ToBarycentric(System& system) {
  auto const semimajor_axis_less_than = [](
      not_null<std::unique_ptr<Subsystem>> const& left,
      not_null<std::unique_ptr<Subsystem>> const& right) -> bool {
    bool const has_semimajor_axes =
        left->jacobi_osculating_elements.semimajor_axis &&
        right->jacobi_osculating_elements.semimajor_axis;
    bool const has_mean_motions =
        left->jacobi_osculating_elements.mean_motion &&
        right->jacobi_osculating_elements.mean_motion;
    bool const has_periods =
        left->jacobi_osculating_elements.period &&
        right->jacobi_osculating_elements.period;
    if (has_semimajor_axes) {
      return left->jacobi_osculating_elements.semimajor_axis <
             right->jacobi_osculating_elements.semimajor_axis;
    }
    if (has_mean_motions) {
      return left->jacobi_osculating_elements.mean_motion >
             right->jacobi_osculating_elements.mean_motion;
    }
    if (has_periods) {
      return left->jacobi_osculating_elements.period <
             right->jacobi_osculating_elements.period;
    }
    LOG(FATAL) << "improperly initialized elements";
    base::noreturn();
  };

  std::sort(system.satellites.begin(),
            system.satellites.end(),
            semimajor_axis_less_than);

  BarycentricSubsystem result;

  // A reference frame wherein the barycentre of |system| is motionless at the
  // origin.
  using SystemBarycentre = geometry::Frame<enum class SystemBarycentreTag>;
  static DegreesOfFreedom<SystemBarycentre> const system_barycentre = {
      SystemBarycentre::origin, SystemBarycentre::unmoving};
  static Identity<SystemBarycentre, Frame> const id_bf;
  static Identity<Frame, SystemBarycentre> const id_fb;

  // Jacobi coordinates for |system|, with satellite subsystems treated as point
  // masses at their barycentres.
  JacobiCoordinates<Frame> jacobi_coordinates(*system.primary);
  // Add the primary first (preorder).
  result.bodies.emplace_back(std::move(system.primary));

  // The |n|th element of |satellite_degrees_of_freedom| contains the list of
  // degrees of freedom of the bodies in the |n|th satellite subsystem with
  // respect to their barycentre.
  std::vector<std::vector<RelativeDegreesOfFreedom<Frame>>>
      satellite_degrees_of_freedom;

  // Fill |satellite_degrees_of_freedom|, |jacobi_coordinates|, and
  // |result.bodies|.
  for (auto const& subsystem : system.satellites) {
    BarycentricSubsystem barycentric_satellite_subsystem =
        ToBarycentric(*subsystem);
    satellite_degrees_of_freedom.emplace_back(std::move(
        barycentric_satellite_subsystem.barycentric_degrees_of_freedom));
    jacobi_coordinates.Add(*barycentric_satellite_subsystem.equivalent_body,
                           subsystem->jacobi_osculating_elements);
    std::move(barycentric_satellite_subsystem.bodies.begin(),
              barycentric_satellite_subsystem.bodies.end(),
              std::back_inserter(result.bodies));
  }

  std::vector<DegreesOfFreedom<SystemBarycentre>> barycentres_of_subsystems;
  {
    // TODO(egg): should BarycentricDegreesOfFreedom return DegreesOfFreedom
    // instead of RelativeDegreesOfFreedom?  In what frame?
    auto const barycentric_dof =
        jacobi_coordinates.BarycentricDegreesOfFreedom();
    for (auto const& dof : barycentric_dof) {
      barycentres_of_subsystems.push_back(system_barycentre + id_fb(dof));
    }
  }

  // Fill |result.barycentric_degrees_of_freedom|.
  // The primary.
  result.barycentric_degrees_of_freedom.emplace_back(
      id_bf(barycentres_of_subsystems.front() - system_barycentre));
  // The bodies in satellite subsystems.
  for (int n = 0; n < satellite_degrees_of_freedom.size(); ++n) {
    // |n + 1| because the primary is at |barycentres_of_subsystems[0]|.
    DegreesOfFreedom<SystemBarycentre> const subsystem_barycentre =
        barycentres_of_subsystems[n + 1];
    for (auto const& body_dof_wrt_subsystem_barycentre :
         satellite_degrees_of_freedom[n]) {
      DegreesOfFreedom<SystemBarycentre> const body_dof =
          subsystem_barycentre + id_fb(body_dof_wrt_subsystem_barycentre);
      result.barycentric_degrees_of_freedom.emplace_back(
          id_bf(body_dof - system_barycentre));
    }
  }

  result.equivalent_body =
      std::make_unique<MassiveBody>(jacobi_coordinates.System());
  return std::move(result);
}

template<typename Frame>
void HierarchicalSystem<Frame>::WriteToMessage(
    std::vector<not_null<std::unique_ptr<Subsystem>>> const& subsystems,
    google::protobuf::RepeatedPtrField<
        serialization::HierarchicalSystem::Subsystem>* const messages) {
  // Sort the subsystems by name to ensure stability of the serialization.
  std::vector<Subsystem const*> sorted_subsystems;
  for (auto const& subsystem : subsystems) {
    sorted_subsystems.push_back(subsystem.get());
  }
  std::sort(sorted_subsystems.begin(),
            sorted_subsystems.end(),
            [](Subsystem const* const lhs, Subsystem const* const rhs) {
              return lhs->primary->name() < rhs->primary->name();
            });

  for (auto const& subsystem : sorted_subsystems) {
    serialization::HierarchicalSystem::Subsystem* const message =
        messages->Add();
    subsystem->primary->WriteToMessage(message->mutable_primary());
    subsystem->jacobi_osculating_elements.WriteToMessage(
        message->mutable_jacobi_osculating_elements());
    WriteToMessage(subsystem->satellites, message->mutable_satellites());
  }
}

}  // namespace internal_hierarchical_system
}  // namespace physics
}  // namespace principia
