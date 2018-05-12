using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class UTF8Marshaler : ICustomMarshaler {
  public abstract void CleanUpNativeData(IntPtr native_data);
  public abstract IntPtr MarshalManagedToNative(object managed_object);
  public abstract object MarshalNativeToManaged(IntPtr native_data);

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  protected readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

// A marshaler for in parameter UTF-8 strings whose ownership is not taken from
// the caller.
internal class InUTF8Marshaler : UTF8Marshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNative(object managed_object) {
    var value = managed_object as String;
    if (value == null) {
      throw Log.Fatal(String.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(String).Name));
    }
    int size = utf8_.GetByteCount(value);
    IntPtr buffer = Marshal.AllocHGlobal(size + 1);
    while (bytes_.Length < size + 1) {
      bytes_ = new byte[2 * bytes_.Length];
    }
    utf8_.GetBytes(value, 0, value.Length, bytes_, 0);
    bytes_[size] = 0;
    Marshal.Copy(bytes_, 0, buffer, size + 1);
    return buffer;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    throw Log.Fatal("use |OutUTF8Marshaler| for out parameters");
  }

  private readonly static InUTF8Marshaler instance_ = new InUTF8Marshaler();
  private byte[] bytes_ = new byte[1];
}

// A marshaler for out parameter or return value UTF-8 strings whose ownership
// is not taken by the caller.
internal class OutUTF8Marshaler : UTF8Marshaler {
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {}
  public override IntPtr MarshalManagedToNative(object managed_object) {
    throw Log.Fatal("use |InUTF8Marshaler| for in parameters");
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    int size;
    for (size = 0; Marshal.ReadByte(native_data, size) != 0; ++size) {}
    while (bytes_.Length < size) {
      bytes_ = new byte[2 * bytes_.Length];
    }
    Marshal.Copy(native_data, bytes_, 0, size);
    return utf8_.GetString(bytes_, 0, size);
  }

  private readonly static OutUTF8Marshaler instance_ = new OutUTF8Marshaler();
  private byte[] bytes_ = new byte[1];
}

// A marshaler for out parameter or return value UTF-8 strings whose ownership
// is taken by the caller.
internal class OutOwnedUTF8Marshaler : OutUTF8Marshaler {
  public static new ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {
    Interface.DeleteString(ref native_data);
  }

  private readonly static OutOwnedUTF8Marshaler instance_ =
      new OutOwnedUTF8Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
