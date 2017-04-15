
#pragma once

#include <typeindex>
#include <type_traits>

#include "base/macros.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"

namespace principia {
namespace interface {

// This is used for interfacing, and should only appear in C++ code in tests
// and generated code; we allow ourselves to pollute the |interface| namespace
// with convenience |using|s.

using base::not_null;
using base::PullSerializer;
using base::PushDeserializer;
using geometry::Instant;
using ksp_plugin::Barycentric;
using ksp_plugin::NavigationFrame;
using ksp_plugin::PileUp;
using ksp_plugin::Plugin;
using ksp_plugin::Vessel;
using ksp_plugin::World;
using physics::DiscreteTrajectory;

// A wrapper for a container and an iterator into that container.
class Iterator {
 public:
  virtual ~Iterator() = default;

  virtual bool AtEnd() const = 0;
  virtual void Increment() = 0;
  virtual int Size() const = 0;
};

// A concrete, typed subclass of |Iterator| which holds a |Container|.
template<typename Container>
class TypedIterator : public Iterator {
 public:
  explicit TypedIterator(Container container);

  // Obtains the element denoted by this iterator and converts it to some
  // |Interchange| type using |convert|.
  template<typename Interchange>
  Interchange Get(
      std::function<Interchange(typename Container::value_type const&)> const&
          convert) const;

  bool AtEnd() const override;
  void Increment() override;
  int Size() const override;

 private:
  Container container_;
  typename Container::const_iterator iterator_;
};

// A specialization for |DiscreteTrajectory<World>|.
template<>
class TypedIterator<DiscreteTrajectory<World>> : public Iterator {
 public:
  TypedIterator(not_null<std::unique_ptr<DiscreteTrajectory<World>>> trajectory,
                not_null<Plugin const*> plugin);

  // Obtains the element denoted by this iterator and converts it to some
  // |Interchange| type using |convert|.
  template<typename Interchange>
  Interchange Get(
      std::function<Interchange(
          DiscreteTrajectory<World>::Iterator const&)> const& convert) const;

  bool AtEnd() const override;
  void Increment() override;
  int Size() const override;

  not_null<Plugin const*> plugin() const;

 private:
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> trajectory_;
  DiscreteTrajectory<World>::Iterator iterator_;
  not_null<Plugin const*> plugin_;
};

// Takes ownership of |**pointer| and returns it to the caller.  Nulls
// |*pointer|.  |pointer| must not be null.  No transfer of ownership of
// |*pointer|.
template<typename T>
std::unique_ptr<T> TakeOwnership(T** pointer);
template<typename T>
std::unique_ptr<T[]> TakeOwnershipArray(T** pointer);

#include "ksp_plugin/interface.generated.h"

extern "C" PRINCIPIA_DLL
void CDECL principia__ActivateRecorder(bool activate);

extern "C" PRINCIPIA_DLL
void CDECL principia__InitGoogleLogging();

bool operator==(AdaptiveStepParameters const& left,
                AdaptiveStepParameters const& right);
bool operator==(Burn const& left, Burn const& right);
bool operator==(NavigationFrameParameters const& left,
                NavigationFrameParameters const& right);
bool operator==(NavigationManoeuvre const& left,
                NavigationManoeuvre const& right);
bool operator==(NavigationManoeuvreFrenetTrihedron const& left,
                NavigationManoeuvreFrenetTrihedron const& right);
bool operator==(QP const& left, QP const& right);
bool operator==(WXYZ const& left, WXYZ const& right);
bool operator==(XYZ const& left, XYZ const& right);

physics::Ephemeris<Barycentric>::AdaptiveStepParameters
FromAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters);
physics::KeplerianElements<Barycentric> FromKeplerianElements(
    KeplerianElements const& keplerian_elements);
geometry::R3Element<double> FromXYZ(XYZ const& xyz);

AdaptiveStepParameters ToAdaptiveStepParameters(
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters);
KeplerianElements ToKeplerianElements(
    physics::KeplerianElements<Barycentric> const& keplerian_elements);
WXYZ ToWXYZ(geometry::Quaternion const& quaternion);
XYZ ToXYZ(geometry::R3Element<double> const& r3_element);

// TODO(phl): These utilities should maybe go into a separate file.
Instant FromGameTime(Plugin const& plugin, double t);
double ToGameTime(Plugin const& plugin, Instant const& t);

// TODO(phl): Return a reference.
not_null<Vessel*> GetVessel(Plugin const& plugin, char const* vessel_guid);

// A factory for NavigationFrame objects.
not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const& plugin,
    NavigationFrameParameters const& parameters);

}  // namespace interface
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
