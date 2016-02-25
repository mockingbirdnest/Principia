using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal class UTF8Marshaler : ICustomMarshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return;
    }
    Marshal.FreeHGlobal(native_data);
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
    IntPtr ptr = Marshal.AllocHGlobal(utf8_.GetByteCount);
    Marshal.StructureToPtr(value, ptr, fDeleteOld: false);
    return ptr;
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    } else {
      return Marshal.PtrToStructure(native_data, typeof(T));
    }
  }

  private readonly static UTF8Marshaler instance_ = new UTF8Marshaler();
  private readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
