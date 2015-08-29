#pragma once

#include <list>
#include <map>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {

using geometry::Instant;

namespace physics {

//TODO(phl): Fix all the comments.
template<typename Tr4jectory, typename TimelineConstIterator_>
class Forkable {
 public:
  using TimelineConstIterator = TimelineConstIterator_;

  Forkable();

  not_null<Tr4jectory*> NewFork(Instant const& time);

  void DeleteFork(not_null<Forkable**> const forkable);

  bool is_root() const;

  not_null<Tr4jectory const*> root() const;
  not_null<Tr4jectory*> root();

  Instant const* ForkTime() const;  // optional

  class Iterator {
   public:
    bool operator==(Iterator const& right) const;
    bool operator!=(Iterator const& right) const;

    Iterator& operator++();

    TimelineConstIterator current() const;

   private:
    bool at_end() const;

    // |ancestry_| has one more element than |forks_|.  The first element in
    // |ancestry_| is the root.  There is no element in |forks_| for the root.
    // It is therefore empty for a root trajectory.
    TimelineConstIterator current_;
    std::list<not_null<Tr4jectory const*>> ancestry_;  // Pointers not owned.

    template<typename Tr4jectory, typename TimelineConstIterator_>
    friend class Forkable;
  };

  Iterator End() const;

  Iterator Find(Instant const& time) const;

  Iterator Wrap(
      not_null<const Tr4jectory*> const ancestor,
      TimelineConstIterator const position_in_ancestor_timeline) const;

 protected:
  virtual not_null<Tr4jectory*> that() = 0;
  virtual not_null<Tr4jectory const*> that() const = 0;

  virtual TimelineConstIterator timeline_begin() const = 0;
  virtual TimelineConstIterator timeline_end() const = 0;
  virtual TimelineConstIterator timeline_find(Instant const& time) const = 0;
  virtual void timeline_insert(TimelineConstIterator begin,
                               TimelineConstIterator end) = 0;
  virtual bool timeline_empty() const = 0;

 private:
  // There may be several forks starting from the same time, hence the multimap.
  using Children = std::multimap<Instant, Tr4jectory>;

  // Null for a root.
  Tr4jectory* parent_;

  // This iterator is only at |end()| for a root.
  //TODO(phl): const? (3x)
  typename Children::const_iterator position_in_parent_children_;

  // This iterator is at |end()| if the parent's timeline is empty, or if this
  // object is a root.
  TimelineConstIterator position_in_parent_timeline_;

  Children children_;
};

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
