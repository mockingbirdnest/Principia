using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

// A marshaler that knows how to encode/decode UTF-16 strings whose ownership is
// not taken from C++.
internal class UnownedUTF16Marshaler : MonoMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    throw Log.Fatal("use |MarshalAs(UnmanagedType.LPWStr)| for in parameters");
    Interface.DeleteU16String(ref native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    throw Log.Fatal("use |MarshalAs(UnmanagedType.LPWStr)| for in parameters");
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    return Marshal.PtrToStringUni(native_data);
  }

  private static readonly UnownedUTF16Marshaler instance_ =
      new UnownedUTF16Marshaler();
}

// A marshaler for UTF-16 strings whose ownership is taken from C++.  Useful for
// out parameters and returned values.
internal class OwnedUTF16Marshaler : UnownedUTF16Marshaler {
  public new static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    var result = base.MarshalNativeToManaged(native_data);
    Interface.DeleteString(ref native_data);
    return result;
  }

  private static readonly OwnedUTF16Marshaler instance_ =
      new OwnedUTF16Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
