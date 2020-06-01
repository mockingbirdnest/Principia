using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// A marshaler that knows how to encode/decode UTF-8 strings.
internal abstract class UTF8Marshaler : MonoMarshaler {
  public object MarshalNativeToManagedImplementation(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    }
    int size;
    for (size = 0; Marshal.ReadByte(native_data, size) != 0; ++size) {}
    byte[] bytes = new byte[size];
    Marshal.Copy(native_data, bytes, 0, size);
    var result = utf8_.GetString(bytes, 0, size);
    return result;
  }

  protected static readonly Encoding utf8_ = new UTF8Encoding(
      encoderShouldEmitUTF8Identifier : false,
      throwOnInvalidBytes             : true);
}

// A marshaler for UTF-8 strings whose ownership is taken from C++.  Useful for
// out parameters and returned values.
internal class OwnershipTransferUTF8Marshaler : UTF8Marshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    throw Log.Fatal("use NoOwnershipTransferUTF8Marshaler for in parameters");
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    throw Log.Fatal("use NoOwnershipTransferUTF8Marshaler for in parameters");
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    var result = MarshalNativeToManagedImplementation(native_data);
    Interface.DeleteString(ref native_data);
    return result;
  }

  private static readonly OwnershipTransferUTF8Marshaler instance_ =
      new OwnershipTransferUTF8Marshaler();
}

// A marshaler for UTF-8 strings whose ownership is not taken from C++.
internal class NoOwnershipTransferUTF8Marshaler : UTF8Marshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (!(managed_object is string value)) {
      throw new NotSupportedException();
    }
    int size = utf8_.GetByteCount(value);
    IntPtr buffer = Marshal.AllocHGlobal(size + 1);
    byte[] bytes = new byte[size + 1];
    utf8_.GetBytes(value, 0, value.Length, bytes, 0);
    bytes[size] = 0;
    Marshal.Copy(bytes, 0, buffer, size + 1);
    return buffer;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    return MarshalNativeToManagedImplementation(native_data);
  }

  private static readonly NoOwnershipTransferUTF8Marshaler instance_ =
      new NoOwnershipTransferUTF8Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
