using System;
using System.Runtime.InteropServices;

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

  public IntPtr IntPtr {
    get {
      return iterator_;
    }
  }

  private IntPtr iterator_ = IntPtr.Zero;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
