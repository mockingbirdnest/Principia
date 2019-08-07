using System;
using System.Globalization;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class OptionalMarshaler<T> : ICustomMarshaler where T : struct {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(string s) {
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
    T value;
    if (managed_object is T value_if_boxed) {
      value = value_if_boxed;
    } else {
      if (managed_object is Boxed<T> value_if_strongly_boxed) {
        value = value_if_strongly_boxed.all;
      } else {
        throw Log.Fatal(
            string.Format(
                CultureInfo.InvariantCulture,
                "|{0}<{1}>| must be used on a boxed |{1}| or on a |{2}<{1}>|.",
                GetType().Name,
                typeof(T).Name,
                typeof(Boxed<>).Name));
      }
    }
    IntPtr ptr = Marshal.AllocHGlobal(Marshal.SizeOf(value));
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

  private static readonly OptionalMarshaler<T> instance_ =
      new OptionalMarshaler<T>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
