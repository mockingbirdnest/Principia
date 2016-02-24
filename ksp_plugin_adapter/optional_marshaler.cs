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
    // While we expect that the marshaling attribute will be used on a |T?|,
    // we get it boxed, and boxing has special behaviour on |Nullable|;
    // specifically, |T?| is boxed to either |T| or |null|, depending on whether
    // it has a value.  In the latter case, we lose the type information, so we
    // cannot test whether |object is T?|.  Instead we check whether
    // |object == null|, if it is not, we check that it's a |T|.
    if (managed_object == null) {
      return IntPtr.Zero;
    }
    var value_if_correct_type = managed_object as T?;
    if (value_if_correct_type == null) {
      throw Log.Fatal(
          String.Format(
              CultureInfo.InvariantCulture,
              "|{0}| must be used on a |{1}| or a (possibly boxed) |{2}|.",
              GetType().Name,
              typeof(T?).Name,
              typeof(T).Name));
    }
    T value = value_if_correct_type.Value;
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

  private readonly static OptionalMarshaler<T> instance_ =
      new OptionalMarshaler<T>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
