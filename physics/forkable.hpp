#pragma once

#include <map>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {

using geometry::Instant;

namespace physics {

//TODO(phl): Fix all the comments.
template<typename Tr4jectory>
class Forkable {

  not_null<Forkable*> NewFork(Instant const& time);

  void DeleteFork(not_null<Forkable**> const fork);

  bool is_root() const;

  not_null<Trajectory const*> root() const;
  not_null<Trajectory*> root();

  Instant const* fork_time() const;

  class Iterator {
   public:
    static Iterator New(not_null<Forkable*> const forkable,
                        Instant const& time);
    static Iterator New(not_null<Forkable*> const forkable,
                        not_null<Forkable*> const ancestor,
                        typename Tr4jectory::TimelineConstIterator const
                            position_in_ancestor_timeline);

    Iterator& operator++();

   private:
    // |ancestry_| has one more element than |forks_|.  The first element in
    // |ancestry_| is the root.  There is no element in |forks_| for the root.
    // It is therefore empty for a root trajectory.
    typename Timeline::const_iterator current_;
    std::list<not_null<Trajectory const*>> ancestry_;  // Pointers not owned.
    std::list<Fork> forks_;
  };

 private:
  // There may be several forks starting from the same time, hence the multimap.
  using Children = std::multimap<Instant, Forkable>;

  // The two iterators denote entries in the containers of the parent.
  // |timeline| is past the end if the fork happened at the fork point of the
  // grandparent.  Note that this implies that the containers should not be
  // swapped.
  struct Fork {
  };

  // Null for a root.
  Forkable* const parent_;

  // At |end()| for a root.
  typename Children::const_iterator position_in_parent_children_;
  typename Tr4jectory::TimelineConstIterator position_in_parent_timeline_;

  Children children_;
};

}  // namespace physics
}  // namespace principia
