using System;
using System.Globalization;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class OptionalMarshaler<T> : MonoMarshaler where T : struct {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
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

  public override object MarshalNativeToManaged(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    } else {
      return Marshal.PtrToStructure(native_data, typeof(T));
    }
  }

  private static readonly OptionalMarshaler<T> instance_ =
      new OptionalMarshaler<T>();
}

internal class OptionalMarshaler<T, TMarshaler> : MonoMarshaler
    where TMarshaler : MonoMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    if (managed_object == null) {
      // This is not our job.
      throw Log.Fatal("The runtime returns null for null objects");
    }
    T value;
    if (managed_object is T value_if_boxed) {
      value = value_if_boxed;
    } else {
      throw Log.Fatal(
          string.Format(
              CultureInfo.InvariantCulture,
              "|{0}<{1}>| must be used on a boxed |{1}| or on a |{2}<{1}>|.",
              GetType().Name,
              typeof(T).Name,
              typeof(Boxed<>).Name));
    }
    return t_marshaler_instance_.MarshalManagedToNativeImplementation(value);
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    if (native_data == IntPtr.Zero) {
      return null;
    } else {
      return t_marshaler_instance_.MarshalNativeToManaged(native_data);
    }
  }

  // The stupid language won't let us access a static method of a generic
  // parameter type, so we go for reflection.
  private static readonly MonoMarshaler t_marshaler_instance_ =
      (MonoMarshaler)typeof(TMarshaler).GetMethod("GetInstance")?.
          Invoke(null, new object[]{ null });
  private static readonly OptionalMarshaler<T, TMarshaler> instance_ =
      new OptionalMarshaler<T, TMarshaler>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
