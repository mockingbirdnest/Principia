using Microsoft.Win32;
using System;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Loader {
  internal static string LoadPrincipiaDllAndInitGoogleLogging() {
    if (loaded_principia_dll_) {
      return null;
    }
    bool is_32_bit = IntPtr.Size == 4;
    if (is_32_bit) {
      return "32-bit platforms are no longer supported; " +
             "use the 64-bit KSP executable.";
    }
    String[] possible_dll_paths = null;
    bool? is_cxx_installed;
    string required_cxx_packages;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        is_cxx_installed = IsVCRedistInstalled();
        required_cxx_packages =
            "the Visual C++ Redistributable Packages for Visual Studio " +
            "2017 on x64";
        possible_dll_paths =
            new String[] {@"GameData\Principia\x64\principia.dll"};
        break;
      // Both Mac and Linux report |PlatformID.Unix|, so we treat them together
      // (we probably don't actually encounter |PlatformID.MacOSX|).
      case PlatformID.Unix:
      case PlatformID.MacOSX:
        possible_dll_paths = new String[] {
            @"GameData/Principia/Linux64/principia.so",
            @"GameData/Principia/MacOS64/principia.so"};
        is_cxx_installed = null;
        required_cxx_packages = "libc++ and libc++abi 6.0-2 (Linux) or " +
                                "El Capitan or later (MacOS)";
        break;
      default:
        return "The operating system " + Environment.OSVersion +
               " is not supported at this time.";
    }
    if (!possible_dll_paths.Any(File.Exists)) {
      return "The principia DLL was not found at '" +
             String.Join("', '", possible_dll_paths) + "'.";
    }
    try {
      loaded_principia_dll_ = true;
      // No kernel32 on *nix, so we throw an exception and immediately resume
      // after the try block.
      try {
        // We dynamically link glog, protobuf, the serialization DLL, and an
        // optimized subset of the physics library on Windows, so we need that
        // to be in the DLL search path for the main DLL to load.
        if (!SetDllDirectory(@"GameData\Principia\x64")) {
          return "Failed to set DLL directory (error code " +
                 Marshal.GetLastWin32Error() + ").";
        }
      } catch {}
      Log.InitGoogleLogging();
      return null;
    } catch (Exception e) {
      UnityEngine.Debug.LogException(e);
      if (is_cxx_installed == false) {
        return "Dependencies, namely " + required_cxx_packages +
               ", were not found.";
      } else {
        return "An unknown error occurred; detected OS " +
               Environment.OSVersion + " 64-bit; tried loading dll at '" +
               String.Join("', '", possible_dll_paths) + "'. Note that " +
               required_cxx_packages + " are required.";
      }
    }
  }

  private static bool IsVCRedistInstalled() {
    // NOTE(phl): This GUID is specific to:
    //   Microsoft Visual C++ 2017 Redistributable (x64) - 14.16.27012
    // It will need to be updated when new versions of Visual C++
    // Redistributable are released by Microsoft.
    RegistryKey key = Registry.LocalMachine.OpenSubKey(
         @"Software\Classes\Installer\Dependencies\" +
         @"VC,redist.x64,amd64,14.16,bundle",
         writable : false);
    if (key == null) {
      return false;
    } else {
      string version = (string)key.GetValue("Version");
      // NOTE(phl): This string needs to be updated when new versions of Visual
      // C++ Redistributable are released by Microsoft.
      return version != null && version == "14.16.27012.6";
    }
  }

  [DllImport("kernel32", CharSet = CharSet.Unicode, SetLastError = true)]
  [return: MarshalAs(UnmanagedType.Bool)]
  static extern bool SetDllDirectory(string lpPathName);

  [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
  private static extern IntPtr LoadLibrary(
      [MarshalAs(UnmanagedType.LPStr)]string lpFileName);

  private const int RTLD_NOW = 2;
  [DllImport("dl")]
  private static extern IntPtr dlopen(
      [MarshalAs(UnmanagedType.LPTStr)] string filename,
      int flags = RTLD_NOW);

  internal static bool loaded_principia_dll_ { get; private set; } = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
