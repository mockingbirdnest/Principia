
#include "ksp_plugin/interface.hpp"

#include <vector>

#include "geometry/rp2_point.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace interface {

using base::check_not_null;
using geometry::RP2Line;
using geometry::RP2Lines;
using geometry::RP2Point;
using ksp_plugin::Camera;
using ksp_plugin::TypedIterator;
using ksp_plugin::VesselSet;
using ksp_plugin::World;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Length;

bool __cdecl principia__IteratorAtEnd(Iterator const* const iterator) {
  journal::Method<journal::IteratorAtEnd> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->AtEnd());
}

void __cdecl principia__IteratorDelete(Iterator** const iterator) {
  journal::Method<journal::IteratorDelete> m({iterator}, {iterator});
  TakeOwnership(iterator);
  return m.Return();
}

QP __cdecl principia__IteratorGetDiscreteTrajectoryQP(
    Iterator const* const iterator) {
  journal::Method<journal::IteratorGetDiscreteTrajectoryQP> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<QP>(
      [](DiscreteTrajectory<World>::Iterator const& iterator) -> QP {
        return ToQP(iterator->degrees_of_freedom);
      }));
}

double __cdecl principia__IteratorGetDiscreteTrajectoryTime(
    Iterator const* const iterator) {
  journal::Method<journal::IteratorGetDiscreteTrajectoryTime> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  auto const plugin = typed_iterator->plugin();
  return m.Return(typed_iterator->Get<double>(
      [plugin](DiscreteTrajectory<World>::Iterator const& iterator) -> double {
        return ToGameTime(*plugin, iterator->time);
      }));
}

XYZ __cdecl principia__IteratorGetDiscreteTrajectoryXYZ(
    Iterator const* const iterator) {
  journal::Method<journal::IteratorGetDiscreteTrajectoryXYZ> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<DiscreteTrajectory<World>> const*>(iterator));
  return m.Return(typed_iterator->Get<XYZ>(
      [](DiscreteTrajectory<World>::Iterator const& iterator) -> XYZ {
        return ToXYZ(iterator->degrees_of_freedom.position());
      }));
}

Iterator* __cdecl principia__IteratorGetRP2LinesIterator(
    Iterator const* const iterator) {
  journal::Method<journal::IteratorGetRP2LinesIterator> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<RP2Lines<Length, Camera>> const*>(iterator));
  return m.Return(typed_iterator->Get<Iterator*>(
      [](RP2Line<Length, Camera> const& rp2_line) -> Iterator* {
        return new TypedIterator<RP2Line<Length, Camera>>(rp2_line);
      }));
}

XY __cdecl principia__IteratorGetRP2LineXY(Iterator const* const iterator) {
  journal::Method<journal::IteratorGetRP2LineXY> m({iterator});
  CHECK_NOTNULL(iterator);
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<RP2Line<Length, Camera>> const*>(iterator));
  return m.Return(typed_iterator->Get<XY>(
      [](RP2Point<Length, Camera> const& rp2_point) -> XY {
        return ToXY(rp2_point);
      }));
}

char const* __cdecl principia__IteratorGetVesselGuid(
    Iterator const* const iterator) {
  journal::Method<journal::IteratorGetVesselGuid> m({iterator});
  auto const typed_iterator = check_not_null(
      dynamic_cast<TypedIterator<VesselSet> const*>(iterator));
  return m.Return(typed_iterator->Get<char const*>(
      [](Vessel* const vessel) -> char const* {
        return vessel->guid().c_str();
      }));
}

void __cdecl principia__IteratorIncrement(Iterator* const iterator) {
  journal::Method<journal::IteratorIncrement> m({iterator});
  CHECK_NOTNULL(iterator)->Increment();
  return m.Return();
}

void __cdecl principia__IteratorReset(Iterator* const iterator) {
  journal::Method<journal::IteratorReset> m({iterator});
  CHECK_NOTNULL(iterator)->Reset();
  return m.Return();
}

int __cdecl principia__IteratorSize(Iterator const* const iterator) {
  journal::Method<journal::IteratorSize> m({iterator});
  return m.Return(CHECK_NOTNULL(iterator)->Size());
}

}  // namespace interface
}  // namespace principia
