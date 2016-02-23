using System;
using System.Globalization;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

public class OptionalMarshaler<T> : ICustomMarshaler where T : struct {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public void CleanUpManagedData(object ManagedObj) {}

  public void CleanUpNativeData(IntPtr pNativeData) {
    if (pNativeData == IntPtr.Zero) {
      return;
    }
    Marshal.FreeHGlobal(pNativeData);
  }

  public int GetNativeDataSize() {
    // I think this is supposed to return -1, and I also think it doesn't
    // matter, but honestly I'm not sure...
    return -1;
  }

  public IntPtr MarshalManagedToNative(object ManagedObj) {
    // While we expect that the marshaling attribute will be used on a |T?|,
    // we get it boxed, and boxing has special behaviour on |Nullable|;
    // specifically, |T?| is boxed to either |T| or |null|, depending on whether
    // it has a value.  In the latter case, we lose the type information, so we
    // cannot test whether |object is T?|.  Instead we check whether
    // |object == null|, if it it's not, we check that it's a |T|.
    if (ManagedObj == null) {
      return IntPtr.Zero;
    }
    var value_if_correct_type = ManagedObj as T?;
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
    Log.Error(value.ToString());
    IntPtr ptr = Marshal.AllocHGlobal(Marshal.SizeOf(value));
    Marshal.StructureToPtr(value, ptr, fDeleteOld: false);
    return ptr;
  }

  public object MarshalNativeToManaged(IntPtr pNativeData) {
    if (pNativeData == IntPtr.Zero) {
      return null;
    } else {
      return Marshal.PtrToStructure(pNativeData, typeof(T));
    }
  }

  private readonly static OptionalMarshaler<T> instance_ =
      new OptionalMarshaler<T>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
