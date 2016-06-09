using Microsoft.Win32;
using System;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Loader {
  public static string LoadPrincipiaDllAndInitGoogleLogging() {
    if (loaded_principia_dll_) {
      return null;
    }
    bool is_32_bit = IntPtr.Size == 4;
    if (is_32_bit) {
      return "32-bit platforms are no longer supported; " +
             "use the 64-bit KSP executable.";
    }
    String[] possible_dll_paths = null;
    bool can_determine_cxx_installed;
    Func<bool> is_cxx_installed;
    string required_cxx_packages;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        can_determine_cxx_installed = true;
        is_cxx_installed = () => IsVCRedistInstalled(is_32_bit);
        required_cxx_packages =
            "the Visual C++ Redistributable Packages for Visual Studio " +
            "2015 on " + (is_32_bit ? "x86" : "x64");
        possible_dll_paths =
            new String[] {@"GameData\Principia\" +
                          (is_32_bit ? "Win32" : "x64") + @"\principia.dll"};
        break;
      // Both Mac and Linux report |PlatformID.Unix|, so we treat them together
      // (we probably don't actually encounter |PlatformID.MacOSX|.
      case PlatformID.Unix:
      case PlatformID.MacOSX:
        possible_dll_paths = new String[] {
            @"GameData/Principia/Linux64/principia.so",
            @"GameData/Principia/MacOS64/principia.so"};
        can_determine_cxx_installed = false;
        is_cxx_installed = null;
        required_cxx_packages = "libc++ and libc++abi 3.5-2";
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
      Log.InitGoogleLogging();
      return null;
    } catch (Exception e) {
      UnityEngine.Debug.LogException(e);
      if (can_determine_cxx_installed && !is_cxx_installed()) {
        return "Dependencies, namely " + required_cxx_packages +
               ", were not found.";
      } else {
        return "An unknown error occurred; detected OS " +
               Environment.OSVersion + " " + (is_32_bit ? "32" : "64") +
               "-bit; tried loading dll at '" +
               String.Join("', '", possible_dll_paths) + "'. Note that " +
               required_cxx_packages + " are required.";
      }
    }
  }

  private static bool IsVCRedistInstalled(bool is_32_bit) {
    // NOTE(phl): These GUIDs are specific to:
    //   Microsoft Visual C++ 2015 Redistributable (x64) - 14.0.23918
    //   Microsoft Visual C++ 2015 Redistributable (x86) - 14.0.23918
    // They will need to be updated when new versions of Visual C++
    // Redistributable are released by Microsoft.
    RegistryKey key = Registry.LocalMachine.OpenSubKey(
         @"Software\Classes\Installer\Dependencies\" +
             (is_32_bit ? "{2e085fd2-a3e4-4b39-8e10-6b8d35f55244}"
                        : "{dab68466-3a7d-41a8-a5cf-415e3ff8ef71}"),
         writable : false);
    if (key == null) {
      return false;
    } else {
      string version = (string)key.GetValue("Version");
      // NOTE(phl): This string needs to be updated when new versions of Visual
      // C++ Redistributable are released by Microsoft.
      return version != null && version == "14.0.23918.0";
    }
  }

  [DllImport("kernel32", SetLastError=true, CharSet = CharSet.Ansi)]
  private static extern IntPtr LoadLibrary(
      [MarshalAs(UnmanagedType.LPStr)]string lpFileName);

  private const int RTLD_NOW = 2;
  [DllImport("dl")]
  private static extern IntPtr dlopen(
      [MarshalAs(UnmanagedType.LPTStr)] string filename,
      int flags = RTLD_NOW);

  private static bool loaded_principia_dll_ = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
