using System;
using System.Collections.Generic;
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
#if __MonoCS__
  internal const string kDllPath = "GameData/Principia/principia.so";
#else
  internal const string kDllPath = "GameData/Principia/principia.dll";
#endif

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__ActivateRecorder",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void ActivateRecorder(bool activate,
                                               bool verbose);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
