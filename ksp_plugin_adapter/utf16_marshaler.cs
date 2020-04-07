using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

// A marshaler that knows how to encode/decode UTF-16 strings.
internal abstract class UTF16Marshaler : MonoMarshaler {
  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    throw Log.Fatal("use |MarshalAs(UnmanagedType.LPWStr)| for in parameters");
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    throw Log.Fatal("use |MarshalAs(UnmanagedType.LPWStr)| for in parameters");
  }

  public object MarshalNativeToManagedImplementation(IntPtr native_data) {
    return Marshal.PtrToStringUni(native_data);
  }
}

// A marshaler for UTF-16 strings whose ownership is taken from C++.  Useful for
// out parameters and returned values.
internal class OwnershipTransferUTF16Marshaler : UTF16Marshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    var result = MarshalNativeToManagedImplementation(native_data);
    Interface.DeleteU16String(ref native_data);
    return result;
  }

  private static readonly OwnershipTransferUTF16Marshaler instance_ =
      new OwnershipTransferUTF16Marshaler();
}

// A marshaler for UTF-16 strings whose ownership is not taken from C++.
internal class NoOwnershipTransferUTF16Marshaler : UTF16Marshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    return Marshal.PtrToStringUni(native_data);
  }

  private static readonly NoOwnershipTransferUTF16Marshaler instance_ =
      new NoOwnershipTransferUTF16Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
