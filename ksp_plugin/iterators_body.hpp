#pragma once

#include <string>
#include <utility>

#include "ksp_plugin/iterators.hpp"

namespace principia {
namespace ksp_plugin {
namespace _iterators {
namespace internal {

template<typename Container>
TypedIterator<Container>::TypedIterator(Container container,
                                        Plugin const* const plugin)
    : container_(std::move(container)),
      iterator_(container_.begin()),
      plugin_(plugin) {}

template<typename Container>
template<typename Interchange>
Interchange TypedIterator<Container>::Get(
    std::function<Interchange(typename Container::value_type const&)> const&
        convert) const {
  CHECK(iterator_ != container_.end());
  return convert(*iterator_);
}

template<typename Container>
bool TypedIterator<Container>::AtEnd() const {
  return iterator_ == container_.end();
}

template<typename Container>
void TypedIterator<Container>::Increment() {
  ++iterator_;
}

template<typename Container>
void TypedIterator<Container>::Reset() {
  iterator_ = container_.begin();
}

template<typename Container>
int TypedIterator<Container>::Size() const {
  return container_.size();
}

template<typename Container>
Plugin const* TypedIterator<Container>::plugin() const {
  return plugin_;
}

inline TypedIterator<DiscreteTrajectory<World>>::TypedIterator(
    DiscreteTrajectory<World> trajectory,
    not_null<Plugin const*> const plugin)
    : trajectory_(std::move(trajectory)),
      iterator_(trajectory_.begin()),
      plugin_(plugin) {}

template<typename Interchange>
Interchange TypedIterator<DiscreteTrajectory<World>>::Get(
    std::function<Interchange(
        DiscreteTrajectory<World>::iterator const&)> const& convert) const {
  CHECK(iterator_ != trajectory_.end());
  return convert(iterator_);
}

inline bool TypedIterator<DiscreteTrajectory<World>>::AtEnd() const {
  return iterator_ == trajectory_.end();
}

inline void TypedIterator<DiscreteTrajectory<World>>::Increment() {
  ++iterator_;
}

inline void TypedIterator<DiscreteTrajectory<World>>::Reset() {
  iterator_ = trajectory_.begin();
}

inline int TypedIterator<DiscreteTrajectory<World>>::Size() const {
  return trajectory_.size();
}

inline DiscreteTrajectory<World>::iterator TypedIterator<
    DiscreteTrajectory<World>>::iterator() const {
  return iterator_;
}

inline not_null<Plugin const*> TypedIterator<
    DiscreteTrajectory<World>>::plugin() const {
  return plugin_;
}

}  // namespace internal
}  // namespace _iterators
}  // namespace ksp_plugin
}  // namespace principia
