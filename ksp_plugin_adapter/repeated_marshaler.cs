using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class RepeatedMarshaler<T, TMarshaler> : MonoMarshaler
    where T : class where TMarshaler : ICustomMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    int i = 0;
    for (;;) {
      IntPtr native_t = Marshal.ReadIntPtr(native_data, i);
      if (native_t == IntPtr.Zero) {
        break;
      }
      t_marshaler_instance_.CleanUpNativeData(native_t);
    }
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (!(managed_object is T[] value)) {
      throw new NotSupportedException();
    }
    int sizeof_t = Marshal.SizeOf(typeof(T));
    IntPtr native_array = Marshal.AllocHGlobal(sizeof_t * (value.Length + 1));
    for (int i = 0; i < value.Length; ++i) {
      IntPtr native_t = t_marshaler_instance_.MarshalManagedToNative(value[i]);
      Marshal.WriteIntPtr(native_array, i, native_t);
    }
    Marshal.WriteIntPtr(native_array, value.Length, IntPtr.Zero);
    return native_array;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    throw new NotSupportedException();
  }

  // The stupid language won't let us access a static method of a generic
  // parameter type, so we go for reflection.
  private static readonly ICustomMarshaler t_marshaler_instance_ =
      (ICustomMarshaler)typeof(TMarshaler).GetMethod("GetInstance")?.
          Invoke(null, new object[]{null});
  private static readonly RepeatedMarshaler<T, TMarshaler> instance_ =
      new RepeatedMarshaler<T, TMarshaler>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
