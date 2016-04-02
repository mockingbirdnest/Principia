#include "ksp_plugin/interface.hpp"

#include <vector>

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using ksp_plugin::LineSegment;

namespace interface {

bool principia__IteratorAtEnd(Iterator const* const iterator) {
  journal::Method<journal::IteratorAtEnd> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->AtEnd());
}

void principia__IteratorDelete(Iterator** const iterator) {
  journal::Method<journal::IteratorDelete> m({iterator}, {iterator});
  TakeOwnership(iterator);
  return m.Return();
}

void principia__IteratorIncrement(Iterator* const iterator) {
  journal::Method<journal::IteratorIncrement> m({iterator});
  CHECK_NOTNULL(iterator)->Increment();
  return m.Return();
}

XYZSegment principia__IteratorGetXYZSegment(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetXYZSegment> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<RenderedTrajectory<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<XYZSegment>(
      [](LineSegment<World> const& line_segment) -> XYZSegment {
    return {ToXYZ((line_segment.begin - World::origin).coordinates() / Metre),
            ToXYZ((line_segment.end - World::origin).coordinates() / Metre)};
  }));
}

int principia__IteratorSize(Iterator const* const iterator) {
  journal::Method<journal::IteratorSize> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->Size());
}


}  // namespace interface
}  // namespace principia
