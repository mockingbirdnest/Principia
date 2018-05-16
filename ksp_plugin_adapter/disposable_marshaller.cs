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
    throw Log.Fatal("No disposable marshaller for in parameters");
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr iterator) {
    return new DisposableIterator(iterator);
  }

  private readonly static DisposableIteratorMarshaller instance_ =
      new DisposableIteratorMarshaller();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
