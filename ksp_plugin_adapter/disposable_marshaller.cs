using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class DisposableIteratorMarshaller : MonoMarshaler {
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

  private static readonly DisposableIteratorMarshaller instance_ =
      new DisposableIteratorMarshaller();
}

internal class DisposablePlanetariumMarshaller : MonoMarshaler {
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

  private static readonly DisposablePlanetariumMarshaller instance_ =
      new DisposablePlanetariumMarshaller();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
