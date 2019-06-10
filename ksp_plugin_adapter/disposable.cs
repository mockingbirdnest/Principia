using System;

namespace principia {
namespace ksp_plugin_adapter {

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
