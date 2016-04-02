#include "ksp_plugin/interface.hpp"

#include <vector>

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using ksp_plugin::Positions;
using ksp_plugin::World;

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

XYZ principia__IteratorGetXYZ(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetXYZ> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<Positions<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<XYZ>(
      [](Position<World> const& position) -> XYZ {
    return ToXYZ((position - World::origin).coordinates() / Metre);
  }));
}

int principia__IteratorSize(Iterator const* const iterator) {
  journal::Method<journal::IteratorSize> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->Size());
}


}  // namespace interface
}  // namespace principia
