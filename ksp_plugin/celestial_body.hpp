#pragma once

#include "ksp_plugin/celestial.hpp"

namespace principia {
namespace ksp_plugin {

inline Celestial::Celestial(not_null<MassiveBody const*> body)
    : body_(body),
      current_time_hint_(
          make_not_null_unique<ContinuousTrajectory<Barycentric>::Hint>) {}

inline bool Celestial::is_initialized() const {
  return trajectory_ != nullptr;
}

inline ContinuousTrajectory<Barycentric> const& Celestial::trajectory() const {
  return *trajectory_;
}

not_null<ContinuousTrajectory<Barycentric>::Hint*>
Celestial::current_time_hint() const {
  return current_time_hint_.get();
}

inline MassiveBody const& Celestial::body() const {
  return *body_;
}

inline bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

inline Celestial const* Celestial::parent() const {
  return parent_;
}

inline void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline void Celestial::WriteToMessage(
    not_null<serialization::Celestial*> const message) const {
  CHECK(is_initialized());
  body_->WriteToMessage(message->mutable_body());
  // TODO(phl): implement.
}

inline std::unique_ptr<Celestial> Celestial::ReadFromMessage(
    serialization::Celestial const& message) {
  // TODO(phl): implement.
  LOG(FATAL) << "Not yet implemented";
}

}  // namespace ksp_plugin
}  // namespace principia
