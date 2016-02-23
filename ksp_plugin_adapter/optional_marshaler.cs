using System;
using System.Globalization;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

class OptionalMarshaler<T> : ICustomMarshaler where T : struct {
  void ICustomMarshaler.CleanUpManagedData(object ManagedObj) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr pNativeData) {
    if (pNativeData == IntPtr.Zero) {
      return;
    }
    Marshal.FreeHGlobal(pNativeData);
  }

  int ICustomMarshaler.GetNativeDataSize() {
    // Not a value type.
    return -1;
  }

  IntPtr ICustomMarshaler.MarshalManagedToNative(object ManagedObj) {
    if (ManagedObj == null) {
      return IntPtr.Zero;
    }
    var value_if_correct_type = ManagedObj as T?;
    if (value_if_correct_type == null) {
      throw new MarshalDirectiveException(
          String.Format(CultureInfo.InvariantCulture,
                        "{0} must be used on a {1} or a {2}.",
                        GetType().Name,
                        typeof(T?).Name,
                        typeof(T).Name));
    }
    T value = value_if_correct_type.Value;
    IntPtr ptr = Marshal.AllocHGlobal(Marshal.SizeOf(value));
    Marshal.StructureToPtr(value, ptr, fDeleteOld: false);
    return ptr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr pNativeData) {
    if (pNativeData == IntPtr.Zero) {
      return null;
    } else {
      return Marshal.PtrToStructure(pNativeData, typeof(T));
    }
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
