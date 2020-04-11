using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

// An adapter that uses TMarshaler to marshal objects of type T, but
// additionally takes ownership of the native data passed in the native-to-
// managed direction.  T may be a struct or a class, but it *must* correspond to
// an interchange message defined in journal.proto.
internal class OwnershipTransferMarshaler<T, TMarshaler> : MonoMarshaler
    where TMarshaler : MonoMarshaler {
  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  public override void CleanUpNativeDataImplementation(IntPtr native_data) {
    t_marshaler_instance_.CleanUpNativeDataImplementation(native_data);
  }

  public override IntPtr MarshalManagedToNativeImplementation(
      object managed_object) {
    return t_marshaler_instance_.MarshalManagedToNativeImplementation(
        managed_object);
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    object managed_object =
        t_marshaler_instance_.MarshalNativeToManaged(native_data);
    Interface.DeleteInterchange(ref native_data);
    return managed_object;
  }

  // The stupid language won't let us access a static method of a generic
  // parameter type, so we go for reflection.
  private static readonly MonoMarshaler t_marshaler_instance_ =
      (MonoMarshaler)typeof(TMarshaler).GetMethod("GetInstance")?.
          Invoke(null, new object[]{null});
  private static readonly OwnershipTransferMarshaler<T, TMarshaler> instance_ =
      new OwnershipTransferMarshaler<T, TMarshaler>();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
