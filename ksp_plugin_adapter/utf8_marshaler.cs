using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

// A marshaler for in parameter UTF-8 strings whose ownership is not taken from
// the caller.
internal class InUTF8Marshaler : ICustomMarshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
    Console.WriteLine("freeh " + Convert.ToString(native_data.ToInt64(), 16));
  }

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  unsafe IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    var value = managed_object as String;
    if (value == null) {
      throw Log.Fatal(String.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(String).Name));
    }
    int length = utf8_.GetByteCount(value);
    IntPtr ptr = Marshal.AllocHGlobal(length + 1);
    byte* bytes = (byte*)ptr;
    Console.WriteLine("alloch " + Convert.ToString(ptr.ToInt64(), 16));
    fixed (char* utf16_units = value) {
      utf8_.GetBytes(utf16_units, value.Length, bytes, length);
    }
    bytes[length] = 0;
    return ptr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr native_data) {
    throw Log.Fatal("use |OutUTF8Marshaler| for out parameters");
  }

  private readonly static InUTF8Marshaler instance_ = new InUTF8Marshaler();
  private readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

// A marshaler for out parameter or return value UTF-8 strings whose ownership
// is not taken by the caller.
internal class OutUTF8Marshaler : ICustomMarshaler {
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    throw Log.Fatal("use |InUTF8Marshaler| for in parameters");
  }

  unsafe object ICustomMarshaler.MarshalNativeToManaged(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    } else {
      sbyte* begin = (sbyte*)native_data;
      sbyte* end;
      for (end = begin; *end != 0; ++end) {}
      return new String(begin, 0, (int)(end - begin), utf8_);
    }
  }

  private readonly static OutUTF8Marshaler instance_ = new OutUTF8Marshaler();
  private readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
