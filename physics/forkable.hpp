#pragma once

#include <list>
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

  void DeleteFork(not_null<Forkable**> const forkable);

  bool is_root() const;

  not_null<Forkable const*> root() const;
  not_null<Forkable*> root();

  bool ContainsTime(Instant const& time) const;  //TODO(phl):Useful?

  Instant const* ForkTime() const;  // optional

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
    typename Tr4jectory::TimelineConstIterator current_;
    std::list<not_null<Forkable const*>> ancestry_;  // Pointers not owned.
    std::list<Fork> forks_;
  };

 protected:
  virtual typename Tr4jectory::TimelineConstIterator timeline_end() = 0;
  virtual typename Tr4jectory::TimelineConstIterator timeline_find(
      Instant const& time) = 0;
  virtual void timeline_insert(
      typename Tr4jectory::TimelineConstIterator begin,
      typename Tr4jectory::TimelineConstIterator end) = 0;
  virtual bool timeline_is_empty() const = 0;

 private:
  // There may be several forks starting from the same time, hence the multimap.
  using Children = std::multimap<Instant, Forkable>;

  // A constructor for creating a child object during forking.
  Forkable(not_null<Forkable*> const parent,
           typename Children::const_iterator position_in_parent_children,
           typename Tr4jectory::TimelineConstIterator
               position_in_parent_timeline);

  // Null for a root.
  Forkable* const parent_;

  // The two iterators denote entries in the containers of the parent.
  // |timeline| is past the end if the fork happened at the fork point of the
  // grandparent.  Note that this implies that the containers should not be
  // swapped.
  // At |end()| for a root.
  typename Children::const_iterator position_in_parent_children_;
  typename Tr4jectory::TimelineConstIterator position_in_parent_timeline_;  // May be at end.

  Children children_;
};

}  // namespace physics
}  // namespace principia

#include "physics/forkable_body.hpp"
