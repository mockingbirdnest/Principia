
#pragma once

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

inline TypedIterator<DiscreteTrajectory<World>>::TypedIterator(
    not_null<std::unique_ptr<DiscreteTrajectory<World>>> trajectory,
    not_null<Plugin const*> const plugin)
    : trajectory_(std::move(trajectory)),
      iterator_(trajectory_->begin()),
      plugin_(plugin) {
  CHECK(trajectory_->is_root());
}

template<typename Interchange>
Interchange TypedIterator<DiscreteTrajectory<World>>::Get(
    std::function<Interchange(
        DiscreteTrajectory<World>::Iterator const&)> const& convert) const {
  CHECK(iterator_ != trajectory_->end());
  return convert(iterator_);
}

inline bool TypedIterator<DiscreteTrajectory<World>>::AtEnd() const {
  return iterator_ == trajectory_->end();
}

inline void TypedIterator<DiscreteTrajectory<World>>::Increment() {
  ++iterator_;
}

inline void TypedIterator<DiscreteTrajectory<World>>::Reset() {
  iterator_ = trajectory_->begin();
}

inline int TypedIterator<DiscreteTrajectory<World>>::Size() const {
  return trajectory_->Size();
}

inline DiscreteTrajectory<World>::Iterator TypedIterator<
    DiscreteTrajectory<World>>::iterator() const {
  return iterator_;
}

inline not_null<Plugin const*> TypedIterator<
    DiscreteTrajectory<World>>::plugin() const {
  return plugin_;
}

}  // namespace ksp_plugin
}  // namespace principia
