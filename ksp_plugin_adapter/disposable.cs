using System;
using System.Collections.Generic;

namespace principia {
namespace ksp_plugin_adapter {

static class DisposableIteratorExtensions {
  public static IEnumerable<Node> Nodes(
      this DisposableIterator nodes_iterator) {
    for (;
         !nodes_iterator.IteratorAtEnd();
         nodes_iterator.IteratorIncrement()) {
      yield return nodes_iterator.IteratorGetNode();
    }
  }

  public static IEnumerable<TQP> DiscreteTrajectoryPoints(
      this DisposableIterator apsis_iterator) {
    for (;
         !apsis_iterator.IteratorAtEnd();
         apsis_iterator.IteratorIncrement()) {
      yield return new TQP {
          t = apsis_iterator.IteratorGetDiscreteTrajectoryTime(),
          qp = apsis_iterator.IteratorGetDiscreteTrajectoryQP()
      };
    }
  }
}

class DisposableIterator : IDisposable {
  public DisposableIterator(IntPtr iterator) {
    iterator_ = iterator;
  }

  ~DisposableIterator() {
    Dispose(false);
  }

  protected virtual void Dispose(bool disposing) {
    if (iterator_ != IntPtr.Zero) {
      Interface.IteratorDelete(ref iterator_);
    }
  }

  public void Dispose() {
    Dispose(true);
    GC.SuppressFinalize(this);
  }

  // Exclusively for use by the marshaller.
  public IntPtr IntPtr => iterator_;

  private IntPtr iterator_ = IntPtr.Zero;
}

class DisposablePlanetarium : IDisposable {
  public DisposablePlanetarium(IntPtr planetarium) {
    planetarium_ = planetarium;
  }

  ~DisposablePlanetarium() {
    Dispose(false);
  }

  protected virtual void Dispose(bool disposing) {
    if (planetarium_ != IntPtr.Zero) {
      Interface.PlanetariumDelete(ref planetarium_);
    }
  }

  public void Dispose() {
    Dispose(true);
    GC.SuppressFinalize(this);
  }

  // Exclusively for use by the marshaller.
  public IntPtr IntPtr => planetarium_;

  private IntPtr planetarium_ = IntPtr.Zero;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
