#pragma once

#include "absl/container/btree_set.h"
#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"

namespace principia {
namespace physics {
namespace internal_clientele {

using namespace principia::base::_not_null;

// A helper class to manage a set of clients to a service.  Clients are ordered
// based on the value of type |Key| passed when they join the |Clientele|.  A
// |Clientele| may contain multiple clients with the same key.
template<typename Key>
class Clientele {
 public:
  // The |default_key| is returned by |first| when the object is empty.
  explicit Clientele(Key const& default_key);

  // Indicates that a client with the given |key| is joining or leaving the
  // clientele.
  void Join(Key const& key);
  void Leave(Key const& key);

  // Returns the smallest key of all clients currently in this object.
  Key const& first() const;

 private:
  Key const default_value_;

  mutable absl::Mutex lock_;
  absl::btree_multiset<Key> clients_;
};

// An RAII object to manage clients that join and leave a |Clientele|.
template<typename Key>
class Client {
 public:
  Client(Key const& key, Clientele<Key>& clientele);
  ~Client();

 private:
  Key const t_;
  not_null<Clientele<Key>*> const clientele_;
};

}  // namespace internal_clientele

using internal_clientele::Client;
using internal_clientele::Clientele;

}  // namespace physics
}  // namespace principia

#include "physics/clientele_body.hpp"
