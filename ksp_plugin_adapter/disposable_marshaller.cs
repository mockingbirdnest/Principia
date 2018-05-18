using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class DisposableIteratorMarshaller : ICustomMarshaler {
  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {}

  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    if (managed_object == null) {
      return IntPtr.Zero;
    }
    var disposable_iterator = managed_object as DisposableIterator;
    return disposable_iterator.IntPtr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr iterator) {
    return new DisposableIterator(iterator);
  }

  private readonly static DisposableIteratorMarshaller instance_ =
      new DisposableIteratorMarshaller();
}

internal class DisposablePlanetariumMarshaller : ICustomMarshaler {
  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {}

  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    if (managed_object == null) {
      return IntPtr.Zero;
    }
    var disposable_planetarium = managed_object as DisposablePlanetarium;
    return disposable_planetarium.IntPtr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr planetarium) {
    return new DisposablePlanetarium(planetarium);
  }

  private readonly static DisposablePlanetariumMarshaller instance_ =
      new DisposablePlanetariumMarshaller();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
