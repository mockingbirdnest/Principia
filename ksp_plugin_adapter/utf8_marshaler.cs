using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class UTF8Marshaler : ICustomMarshaler {
  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return;
    }
    Marshal.FreeHGlobal(native_data);
    Console.WriteLine("freeh " + Convert.ToString(native_data.ToInt64(), 16));
  }

  int ICustomMarshaler.GetNativeDataSize() {
    // I think this is supposed to return -1, and I also think it doesn't
    // matter, but honestly I'm not sure...
    return -1;
  }

  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    if (managed_object == null) {
      // This is not our job.
      throw Log.Fatal("The runtime returns null for null objects");
    }
    var value = managed_object as String;
    if (value == null) {
      throw Log.Fatal(String.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(String).Name));
    }
    int length = utf8_.GetByteCount(value);
    IntPtr ptr = Marshal.AllocHGlobal(length + 1);
    Console.WriteLine("alloch " + Convert.ToString(ptr.ToInt64(), 16));
    Marshal.Copy(utf8_.GetBytes(value), 0, ptr, length);
    Marshal.WriteByte(ptr, length, 0);
    return ptr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    } else {
      int length = 0;
      while(Marshal.ReadByte(native_data, length) != 0) { ++length; }
      byte[] bytes = new byte[length];
      Marshal.Copy(native_data, bytes, 0, length);
      string s = utf8_.GetString(bytes);
      Console.WriteLine(s);
      return utf8_.GetString(bytes);
    }
  }

  private readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

internal class InUTF8Marshaler : UTF8Marshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  private readonly static InUTF8Marshaler instance_ =
      new InUTF8Marshaler();
}

internal class OutUTF8Marshaler : UTF8Marshaler, ICustomMarshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {}

  // Don't leak.
  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    throw Log.Fatal("use |InUTF8Marshaler| for in parameters");
  }

  private readonly static OutUTF8Marshaler instance_ =
      new OutUTF8Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
