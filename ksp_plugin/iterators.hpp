
#pragma once

#include "base/not_null.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {

using base::not_null;
using physics::DiscreteTrajectory;

// A wrapper for a container and an iterator into that container.
class Iterator {
 public:
  virtual ~Iterator() = default;

  virtual bool AtEnd() const = 0;
  virtual void Increment() = 0;
  virtual void Reset() = 0;
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
  void Reset() override;
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
  void Reset() override;
  int Size() const override;

  DiscreteTrajectory<World>::Iterator iterator() const;
  not_null<Plugin const*> plugin() const;

 private:
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> trajectory_;
  DiscreteTrajectory<World>::Iterator iterator_;
  not_null<Plugin const*> plugin_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/iterators_body.hpp"
