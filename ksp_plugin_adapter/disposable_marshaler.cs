using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class DisposableIteratorMarshaler : MonoMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {}

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (managed_object == null) {
      return IntPtr.Zero;
    }
    var disposable_iterator = managed_object as DisposableIterator;
    return disposable_iterator.IntPtr;
  }

  public override object MarshalNativeToManaged(IntPtr iterator) {
    return new DisposableIterator(iterator);
  }

  private static readonly DisposableIteratorMarshaler instance_ =
      new DisposableIteratorMarshaler();
}

internal class DisposablePlanetariumMarshaler : MonoMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {}

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (managed_object == null) {
      return IntPtr.Zero;
    }
    var disposable_planetarium = managed_object as DisposablePlanetarium;
    return disposable_planetarium.IntPtr;
  }

  public override object MarshalNativeToManaged(IntPtr planetarium) {
    return new DisposablePlanetarium(planetarium);
  }

  private static readonly DisposablePlanetariumMarshaler instance_ =
      new DisposablePlanetariumMarshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
