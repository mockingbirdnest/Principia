#include "ksp_plugin/interface.hpp"

#include <vector>

#include "journal/method.hpp"
#include "journal/profiles.hpp"

namespace principia {
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

XYZSegment CDECL
principia__IteratorGetXYZSegment(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetXYZSegment> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const* typed_iterator =
      dynamic_cast<TypedIterator<XYZSegment, std::vector> const*>(iterator);
  return m.Return(typed_iterator->Get());
}

int principia__IteratorSize(Iterator* const iterator) {
  journal::Method<journal::IteratorSize> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->Size());
}


}  // namespace interface
}  // namespace principia
