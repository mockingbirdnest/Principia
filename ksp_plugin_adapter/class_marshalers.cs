using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static partial class Interface {
  private static IntPtr At(this IntPtr pointer, long offset) {
    return new IntPtr(pointer.ToInt64() + offset);
  }

  internal class InBodyParametersMarshaler : ICustomMarshaler {

    [StructLayout(LayoutKind.Sequential)]
    internal class BodyParametersRepresentation {
      public string name;
      public string gravitational_parameter;
      public string reference_instant;
      public string axis_right_ascension;
      public string axis_declination;
      public string reference_angle;
      public string angular_frequency;
      public string reference_radius;
      public string j2;
      public double min_radius;
      public double mean_radius;
      public double max_radius;
      public IntPtr geopotential;
      public int geopotential_size;
    }

    public static ICustomMarshaler GetInstance(string s) {
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
          max_radius              = parameters.max_radius,
          mean_radius             = parameters.mean_radius,
          min_radius              = parameters.min_radius,
          name                    = parameters.name,
          reference_angle         = parameters.reference_angle,
          reference_instant       = parameters.reference_instant,
          reference_radius        = parameters.reference_radius,
          geopotential_size       = parameters.geopotential?.Length ?? 0
      };
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

    private static readonly InBodyParametersMarshaler instance_ =
        new InBodyParametersMarshaler();
  }

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
