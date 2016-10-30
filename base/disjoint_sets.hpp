#pragma once

#include <list>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// For the purposes of this class, |T| represents the set of its values, and
// a single globally unique partition is being built.  If |MakeSingleton| is
// called on an element |e| of type |T|, all properties of the subset previously
// containing |e| are invalidated.

template<typename T>
class Subset;

template<typename T>
Subset<T> MakeSingleton(not_null<T*> element);
// The arguments are invalidated; the result may be used to get information
// about the united subset.
template<typename T>
Subset<T> Unite(Subset<T> left, Subset<T> right);
// Returns the subset containing |element|.
template<typename T>
Subset<T> Find(T& element);

// Any properties about a subset of |T| that can be efficiently maintained when
// merging (e.g. a list of elements) should be kept in |SubsetProperties<T>|;
// specialize it as needed.
template<typename T>
class SubsetProperties {
 public:
  void MergeWith(SubsetProperties& other) {}
};

template<typename T>
class SubsetNode;

// A subset of |T|.
template<typename T>
class Subset {
 public:
  SubsetProperties<T> const& properties();

 private:
  explicit Subset(not_null<SubsetNode<T>*> node);

  not_null<SubsetNode<T>*> const node_;

  friend Subset<T> MakeSingleton<>(not_null<T*> element);
  friend Subset<T> Unite<>(Subset<T> left, Subset<T> right);
  friend Subset<T> Find<>(T& element);

  friend bool operator==(Subset left, Subset right) {
    return left.node_ == right.node_;
  }
};

template<typename T>
class SubsetNode {
 public:
  SubsetNode();
  SubsetProperties<T> properties_;
 private:
  not_null<SubsetNode<T>*> Root();

  // TODO(egg): maybe they should be mutable? those intrusive structures are
  // confusing...
  not_null<SubsetNode<T>*> parent_;
  int rank_ = 0;
  friend Subset<T> MakeSingleton<>(not_null<T*> element);
  friend Subset<T> Unite<>(Subset<T> left, Subset<T> right);
  friend Subset<T> Find<>(T& element);
};

// Specialize to return a |SubsetNode| owned by |element| (in constant time).
template<typename T>
not_null<SubsetNode<T>*> GetSubsetNode(T& element) {
  static_assert(false, "GetSubsetNode must be specialized");
}

}  // namespace base
}  // namespace principia

#include "base/disjoint_sets_body.hpp"
