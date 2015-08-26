#pragma once

#include <list>
#include <map>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {

using geometry::Instant;

namespace physics {

//TODO(phl): Fix all the comments.
template<typename TimelineConstIterator_>
class Forkable {
 public:
  using TimelineConstIterator = TimelineConstIterator_;

  not_null<Forkable*> NewFork(Instant const& time);

  void DeleteFork(not_null<Forkable**> const forkable);

  bool is_root() const;

  not_null<Forkable const*> root() const;
  not_null<Forkable*> root();

  Instant const* ForkTime() const;  // optional

  class Iterator {
   public:
    static Iterator New(not_null<Forkable*> const forkable,
                        Instant const& time);
    static Iterator New(not_null<Forkable*> const forkable,
                        not_null<Forkable*> const ancestor,
                        TimelineConstIterator const 
                            position_in_ancestor_timeline);

    Iterator& operator++();

    TimelineConstIterator current() const;

   private:
    // |ancestry_| has one more element than |forks_|.  The first element in
    // |ancestry_| is the root.  There is no element in |forks_| for the root.
    // It is therefore empty for a root trajectory.
    TimelineConstIterator current_;
    std::list<not_null<Forkable const*>> ancestry_;  // Pointers not owned.
  };

 protected:
  virtual TimelineConstIterator timeline_end() const = 0;
  virtual TimelineConstIterator timeline_find(Instant const& time) const = 0;
  virtual void timeline_insert(TimelineConstIterator begin,
                               TimelineConstIterator end) = 0;
  virtual bool timeline_empty() const = 0;

 private:
  // There may be several forks starting from the same time, hence the multimap.
  using Children = std::multimap<Instant, Forkable>;

  // A constructor for creating a child object during forking.
  Forkable(not_null<Forkable*> const parent,
           typename Children::const_iterator position_in_parent_children,
           TimelineConstIterator position_in_parent_timeline);

  // Null for a root.
  Forkable* const parent_;

  // This iterator is only at |end()| for a root.
  typename Children::const_iterator position_in_parent_children_;

  // This iterator is at |end()| if the parent's timeline is empty, or if this
  // object is a root.
  TimelineConstIterator position_in_parent_timeline_;

  Children children_;
};

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
