
#pragma once

#include <optional>
#include <list>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// For the purposes of this class, |T| represents the set of its values, and
// a single globally unique partition is being built.  If |MakeSingleton| is
// called on an element |e| of type |T|, all properties of the subset previously
// containing |e| are invalidated.
// To use an union-find algorithm on elements of |T|, specialize
// |Subset<T>::Node::Get|, run |Subset<T>::MakeSingleton| on all elements
// involved, and proceed with calls to |Subset<T>::Unite| and |Subset<T>::Find|.

// A subset of |T|.
template<typename T>
class Subset final {
 public:
  // Any properties about a subset of |T| that can be efficiently maintained
  // when merging (e.g. a list of elements) should be kept in
  // |Subset<T>::Properties|; specialize it as needed.
  class Properties final {
   public:
    void MergeWith(Properties& other) {}
  };

  // The |SubsetPropertiesArgs| are forwarded to the constructor of
  // |Properties|.
  template<typename... SubsetPropertiesArgs>
  static Subset MakeSingleton(
      T& element,
      SubsetPropertiesArgs... subset_properties_args);

  // The arguments are invalidated; the result may be used to get information
  // about the united subset.
  static Subset Unite(Subset left, Subset right);
  // Returns the subset containing |element|.
  static Subset Find(T& element);

  Properties const& properties() const;
  Properties& mutable_properties();

  class Node final {
   public:
    Node();

   private:
    // Specialize to return a |Node| owned by |element| (in constant time).  The
    // compiler will warn about returning from a non-void function if this is
    // not specialized.
    static not_null<typename Subset::Node*> Get(T& element) {}

    not_null<Node*> Root();

    not_null<Node*> parent_;
    int rank_ = 0;

    // Do not require a default constructor for |Properties|.
    std::optional<Properties> properties_;

    friend class Subset<T>;
  };

  friend bool operator==(Subset left, Subset right) {
    return left.node_ == right.node_;
  }

 private:
  explicit Subset(not_null<Node*> node);

  not_null<Node*> const node_;
};

}  // namespace base
}  // namespace principia

#include "base/disjoint_sets_body.hpp"
