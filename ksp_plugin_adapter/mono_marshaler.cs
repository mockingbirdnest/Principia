using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

// It seems that, for out parameters and returned value, Mono calls
// CleanUpNativeData on data that was not allocated by MarshalManagedToNative,
// see https://github.com/mono/mono/blob/941a335ea0f20c22a02a7947945f53787a56b2d3/mono/metadata/marshal-ilgen.c#L4555
// This appears to violate the contract documented by .Net, see
// https://docs.microsoft.com/en-us/dotnet/api/system.runtime.interopservices.icustommarshaler?view=netframework-4.8#implementing-the-icustommarshaler-interface
// although in reality it seems that the .Net tests themselves don't expect
// CleanUpManagedData and GetNativeDataSize to be called, see
// https://github.com/dotnet/runtime/blob/4f9ae42d861fcb4be2fcd5d3d55d5f227d30e723/src/coreclr/tests/src/Interop/ICustomMarshaler/Primitives/ICustomMarshaler.cs#L147.
// Amid all this confusion, this class attempts to maintain a sane behavior,
// i.e., one that is close enough to the .Net documentation.
internal abstract class MonoMarshaler : ICustomMarshaler {
  // Subclasses must override the following methods to implement the contract
  // defined by .Net.
  public abstract void CleanUpNativeDataImplementation(IntPtr native_data);
  public abstract IntPtr MarshalManagedToNativeImplementation(
      object managed_object);

  // We have no evidence that this method is ever called.
  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {
    IntPtr actual_native_data = IntPtr.Zero;
    lock (allocated_intptrs_) {
      if (allocated_intptrs_.Contains(native_data)) {
        actual_native_data = native_data;
        allocated_intptrs_.Remove(native_data);
      }
    }
    if (actual_native_data != IntPtr.Zero) {
      CleanUpNativeDataImplementation(actual_native_data);
    }
  }

  // We have no evidence that this method is ever called.
  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    IntPtr result = MarshalManagedToNativeImplementation(managed_object);
    lock (allocated_intptrs_) {
      if (result != IntPtr.Zero) {
        allocated_intptrs_.Add(result);
      }
    }
    return result;
  }

  public abstract object MarshalNativeToManaged(IntPtr native_data);

  private static readonly HashSet<IntPtr> allocated_intptrs_ =
      new HashSet<IntPtr>();
}

} // namespace ksp_plugin_adapter
}  // namespace principia
