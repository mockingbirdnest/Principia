using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal partial struct XYZ {
  public static explicit operator XYZ(Vector3d v) {
    return new XYZ{x = v.x, y = v.y, z = v.z};
  }
  public static explicit operator Vector3d(XYZ v) {
    return new Vector3d{x = v.x, y = v.y, z = v.z};
  }
}

internal partial struct WXYZ {
  public static explicit operator WXYZ(UnityEngine.QuaternionD q) {
    return new WXYZ{w = q.w, x = q.x, y = q.y, z = q.z};
  }
  public static explicit operator UnityEngine.QuaternionD(WXYZ q) {
    return new UnityEngine.QuaternionD{w = q.w, x = q.x, y = q.y, z = q.z};
  }
}

internal static partial class Interface {
  internal const string dll_path = "principia";

  internal static KeplerianElements Elements(this Orbit orbit) {
    double mean_motion = 2 * Math.PI / orbit.period;
    return new KeplerianElements{
        eccentricity                           = orbit.eccentricity,
        semimajor_axis                         = double.NaN,
        mean_motion                            = mean_motion,
        inclination_in_degrees                 = orbit.inclination,
        longitude_of_ascending_node_in_degrees = orbit.LAN,
        argument_of_periapsis_in_degrees       = orbit.argumentOfPeriapsis,
        mean_anomaly                           =
            orbit.meanAnomalyAtEpoch - orbit.epoch * mean_motion};
  }

  [DllImport(dllName           : dll_path,
             EntryPoint        = "principia__ActivateRecorder",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void ActivateRecorder(bool activate);

  [DllImport(dllName           : dll_path,
             EntryPoint        = "principia__InitGoogleLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void InitGoogleLogging();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
