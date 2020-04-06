using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// A marshaler that knows how to encode/decode UTF-8 strings.
internal abstract class UTF8Marshaler : MonoMarshaler {
  protected static readonly Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

// A marshaler for UTF-8 strings whose ownership is not taken from C++.
internal class UnownedUTF8Marshaler : UTF8Marshaler {
  public static ICustomMarshaler GetInstance(string s) {
    UnityEngine.Debug.LogError("UnownedUTF8Marshaler.GetInstance");
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    UnityEngine.Debug.LogError("UnownedUTF8Marshaler.CleanUpNativeData " +
                               GetHashCode() + " " + native_data);
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (!(managed_object is string value)) {
      throw Log.Fatal(string.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(string).Name));
    }
    int size = utf8_.GetByteCount(value);
    IntPtr buffer = Marshal.AllocHGlobal(size + 1);
    byte[] bytes = new byte[size + 1];
    utf8_.GetBytes(value, 0, value.Length, bytes, 0);
    bytes[size] = 0;
    Marshal.Copy(bytes, 0, buffer, size + 1);
    UnityEngine.Debug.LogError("UnownedUTF8Marshaler.MarshalManagedToNative " +
                               GetHashCode() + " " + buffer + " " +
                               managed_object);
    return buffer;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    int size;
    for (size = 0; Marshal.ReadByte(native_data, size) != 0; ++size) {}
    byte[] bytes = new byte[size];
    Marshal.Copy(native_data, bytes, 0, size);
    var result = utf8_.GetString(bytes, 0, size);
    UnityEngine.Debug.LogError("UnownedUTF8Marshaler.MarshalNativeToManaged " +
                               GetHashCode() + " " + native_data + " " +
                               result);
    return result;
  }

  private static readonly UnownedUTF8Marshaler instance_ =
      new UnownedUTF8Marshaler();
}

// A marshaler for UTF-8 strings whose ownership is taken from C++.  Useful for
// out parameters and returned values.
internal class OwnedUTF8Marshaler : UnownedUTF8Marshaler {
  public new static ICustomMarshaler GetInstance(string s) {
    UnityEngine.Debug.LogError("OwnedUTF8Marshaler.GetInstance");
    return instance_;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    var result = base.MarshalNativeToManaged(native_data);
    Interface.DeleteString(ref native_data);
    UnityEngine.Debug.LogError("OwnedUTF8Marshaler.MarshalNativeToManaged " +
                               GetHashCode() + " " + native_data + " " +
                               result);
    return result;
  }

  private static readonly OwnedUTF8Marshaler instance_ =
      new OwnedUTF8Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
