using System;
using System.IO;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Loader {
  [DllImport("kernel32", SetLastError=true, CharSet = CharSet.Ansi)]
  static extern IntPtr LoadLibrary(
      [MarshalAs(UnmanagedType.LPStr)]string lpFileName);

  const int RTLD_NOW = 2;
  [DllImport("dl")]
  static extern IntPtr dlopen(
      [MarshalAs(UnmanagedType.LPTStr)] string filename,
      int flags = RTLD_NOW);


  // TODO(egg): return a string detailing the error.
  public static bool LoadPrincipiaDll() {
    if (loaded_principia_dll_) {
      return true;
    }
    bool is_32_bit = IntPtr.Size == 4;
    String dll = null;
    Func<String, IntPtr> load;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        load = LoadLibrary;
        if (is_32_bit) {
          dll = @"GameData\Principia\Win32\principia.dll";
        } else {
          dll = @"GameData\Principia\x64\principia.dll";
        }
        break;
      case PlatformID.Unix:
        load = (name) => dlopen(name);
        if (is_32_bit) {
          return false;
        } else {
          dll = @"GameData/Principia/Linux64/principia.so";
        }
        break;
      default:
        return false;
    }
    if (!File.Exists(dll)) {
      return false;
    }
    load(dll);
    loaded_principia_dll_ = true;
    return true;
  }

  private static bool loaded_principia_dll_ = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
