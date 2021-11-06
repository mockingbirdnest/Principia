
#pragma once

#include <string>

#include "ksp_plugin/iterators.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Container>
TypedIterator<Container>::TypedIterator(Container container)
    : container_(std::move(container)),
      iterator_(container_.begin()) {}

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

inline TypedIterator<DiscreteTraject0ry<World>>::TypedIterator(
    DiscreteTraject0ry<World> trajectory,
    not_null<Plugin const*> const plugin)
    : trajectory_(std::move(trajectory)),
      iterator_(trajectory_.begin()),
      plugin_(plugin) {}

template<typename Interchange>
Interchange TypedIterator<DiscreteTraject0ry<World>>::Get(
    std::function<Interchange(
        DiscreteTraject0ry<World>::iterator const&)> const& convert) const {
  CHECK(iterator_ != trajectory_.end());
  return convert(iterator_);
}

inline bool TypedIterator<DiscreteTraject0ry<World>>::AtEnd() const {
  return iterator_ == trajectory_.end();
}

inline void TypedIterator<DiscreteTraject0ry<World>>::Increment() {
  ++iterator_;
}

inline void TypedIterator<DiscreteTraject0ry<World>>::Reset() {
  iterator_ = trajectory_.begin();
}

inline int TypedIterator<DiscreteTraject0ry<World>>::Size() const {
  return trajectory_.size();
}

inline DiscreteTraject0ry<World>::iterator TypedIterator<
    DiscreteTraject0ry<World>>::iterator() const {
  return iterator_;
}

inline not_null<Plugin const*> TypedIterator<
    DiscreteTraject0ry<World>>::plugin() const {
  return plugin_;
}

}  // namespace ksp_plugin
}  // namespace principia
