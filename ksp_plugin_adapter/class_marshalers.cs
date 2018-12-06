using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static partial class Interface {

  static IntPtr At(this IntPtr pointer, Int64 offset) {
    return new IntPtr(pointer.ToInt64() + offset);
  }

  internal class InBodyParametersMarshaler : ICustomMarshaler {

    [StructLayout(LayoutKind.Sequential)]
    internal partial class BodyParametersRepresentation {
      public String name;
      public String gravitational_parameter;
      public String reference_instant;
      public String mean_radius;
      public String axis_right_ascension;
      public String axis_declination;
      public String reference_angle;
      public String angular_frequency;
      public String reference_radius;
      public String j2;
      public IntPtr geopotential;
      public int geopotential_size;
    }

    public static ICustomMarshaler GetInstance(String s) {
      return instance_;
    }

    public void CleanUpNativeData(IntPtr native_data) {
      var representation = new BodyParametersRepresentation();
      Marshal.PtrToStructure(native_data, representation);
      for (int i = 0; i < representation.geopotential_size; ++i) {
        Marshal.DestroyStructure(
            representation.geopotential.At(
                i * Marshal.SizeOf(typeof(BodyGeopotentialElement))),
            typeof(BodyGeopotentialElement));
      }
      Marshal.FreeHGlobal(representation.geopotential);
      Marshal.FreeHGlobal(native_data);
    }

    public IntPtr MarshalManagedToNative(object managed_object) {
      var parameters = managed_object as BodyParameters;
      var representation = new BodyParametersRepresentation{
          angular_frequency       = parameters.angular_frequency,
          axis_declination        = parameters.axis_declination,
          axis_right_ascension    = parameters.axis_right_ascension,
          gravitational_parameter = parameters.gravitational_parameter,
          j2                      = parameters.j2,
          mean_radius             = parameters.mean_radius,
          name                    = parameters.name,
          reference_angle         = parameters.reference_angle,
          reference_instant       = parameters.reference_instant,
          reference_radius        = parameters.reference_radius};
      representation.geopotential_size = parameters.geopotential?.Length ?? 0;
      if (representation.geopotential_size == 0) {
        representation.geopotential = IntPtr.Zero;
      } else {
        int sizeof_element = Marshal.SizeOf(typeof(BodyGeopotentialElement));
        representation.geopotential = Marshal.AllocHGlobal(
            sizeof_element * parameters.geopotential.Length);
        for (int i = 0; i < parameters.geopotential.Length; ++i) {
          Marshal.StructureToPtr(
              parameters.geopotential[i],
              representation.geopotential.At(i * sizeof_element),
              fDeleteOld: false);
        }
      }
      IntPtr buffer = Marshal.AllocHGlobal(Marshal.SizeOf(representation));
      Marshal.StructureToPtr(representation, buffer, fDeleteOld: false);
      return buffer;
    }

    public object MarshalNativeToManaged(IntPtr native_data) {
      throw Log.Fatal("InBodyParametersMarshaler.MarshalNativeToManaged");
    }

    public void CleanUpManagedData(object managed_data) {
      throw Log.Fatal("InBodyParametersMarshaler.CleanUpManagedData");
    }

    int ICustomMarshaler.GetNativeDataSize() {
      return -1;
    }

    private readonly static InBodyParametersMarshaler instance_ =
        new InBodyParametersMarshaler();
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
