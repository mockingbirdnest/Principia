#pragma once

#include "physics/clientele.hpp"

namespace principia {
namespace physics {
namespace _clientele {
namespace internal {

template<typename Key>
Clientele<Key>::Clientele(Key const& default_value)
    : default_value_(default_value) {}

template<typename Key>
void Clientele<Key>::Join(Key const& key) {
  absl::MutexLock l(&lock_);
  clients_.insert(key);
}

template<typename Key>
void Clientele<Key>::Leave(Key const& key) {
  absl::MutexLock l(&lock_);
  auto const it = clients_.find(key);
  CHECK(it != clients_.end()) << key;
  clients_.erase(it);
}

template<typename Key>
Key const& Clientele<Key>::first() const {
  absl::MutexLock l(&lock_);
  if (auto const it = clients_.begin(); it == clients_.end()) {
    return default_value_;
  } else {
    return *it;
  }
}

template<typename Key>
Client<Key>::Client(Key const& key, Clientele<Key>& clientele)
    : t_(key), clientele_(&clientele) {
  clientele_->Join(t_);
}

template<typename Key>
Client<Key>::~Client() {
  clientele_->Leave(t_);
}

}  // namespace internal
}  // namespace _clientele
}  // namespace physics
}  // namespace principia
