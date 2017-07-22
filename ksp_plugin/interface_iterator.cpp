
#include "ksp_plugin/interface.hpp"

#include <vector>

#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace interface {

using base::check_not_null;
using ksp_plugin::TypedIterator;
using ksp_plugin::World;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;

bool principia__IteratorAtEnd(Iterator const* const iterator) {
  journal::Method<journal::IteratorAtEnd> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->AtEnd());
}

void principia__IteratorDelete(Iterator** const iterator) {
  journal::Method<journal::IteratorDelete> m({iterator}, {iterator});
  TakeOwnership(iterator);
  return m.Return();
}

QP principia__IteratorGetQP(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetQP> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<QP>(
      [](DiscreteTrajectory<World>::Iterator const& iterator) -> QP {
        return ToQP(iterator.degrees_of_freedom());
      }));
}

double principia__IteratorGetTime(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetTime> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  auto const plugin = typed_iterator->plugin();
  return m.Return(typed_iterator->Get<double>(
      [plugin](DiscreteTrajectory<World>::Iterator const& iterator) -> double {
        return ToGameTime(*plugin, iterator.time());
      }));
}

XYZ principia__IteratorGetXYZ(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetXYZ> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<XYZ>(
      [](DiscreteTrajectory<World>::Iterator const& iterator) -> XYZ {
        return ToXYZ(iterator.degrees_of_freedom().position());
      }));
}

void principia__IteratorIncrement(Iterator* const iterator) {
  journal::Method<journal::IteratorIncrement> m({iterator});
  CHECK_NOTNULL(iterator)->Increment();
  return m.Return();
}

int principia__IteratorSize(Iterator const* const iterator) {
  journal::Method<journal::IteratorSize> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->Size());
}

}  // namespace interface
}  // namespace principia
