using Microsoft.Win32;
using System;
using System.IO;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Loader {
  public static string LoadPrincipiaDllAndInitGoogleLogging() {
    if (loaded_principia_dll_) {
      return null;
    }
    bool is_32_bit = IntPtr.Size == 4;
    String dll = null;
    Func<String, IntPtr> load;
    Func<bool> is_cxx_installed;
    string required_cxx_packages;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        load = LoadLibrary;
        is_cxx_installed = () => IsVCRedistInstalled(is_32_bit);
        if (is_32_bit) {
          required_cxx_packages =
              "the Visual C++ Redistributable Packages for Visual Studio " +
              "2015 on x86";
          dll = @"GameData\Principia\Win32\principia.dll";
        } else {
          required_cxx_packages =
              "the Visual C++ Redistributable Packages for Visual Studio " +
              "2015 on x64";
          dll = @"GameData\Principia\x64\principia.dll";
        }
        break;
      case PlatformID.Unix:
        load = (name) => dlopen(name);
        if (is_32_bit) {
          return "Linux 32-bit is not supported at this time.";
        } else {
          dll = @"GameData/Principia/Linux64/principia.so";
          // TODO(egg): figure out how to check whether the right versions of
          // libc++ and libc++abi are installed.
          is_cxx_installed = () => false;
          required_cxx_packages = "libc++ and libc++abi 3.5-2";
        }
        break;
      default:
        return "The operating system " + Environment.OSVersion +
               " is not supported at this time.";
    }
    if (!File.Exists(dll)) {
      return "The principia DLL was not found at '" + dll + "'.";
    }
    try {
      load(dll);
      loaded_principia_dll_ = true;
      Log.InitGoogleLogging();
      return null;
    } catch (Exception e) {
      UnityEngine.Debug.LogException(e);
      if (!is_cxx_installed()) {
        return "Dependencies, namely " + required_cxx_packages +
               ", were not found.";
      } else {
        return "An unknown error occurred; detected OS " +
               Environment.OSVersion + " " + (is_32_bit ? "32" : "64") +
               "-bit; tried loading dll at '" + dll + "'.";
      }
    }
  }

  private static bool IsVCRedistInstalled(bool is_32_bit) {
     RegistryKey key = Registry.LocalMachine.OpenSubKey(
         @"Software\Microsoft\VisualStudio\14.0\VC\Runtimes" +
             (is_32_bit ? "x86" : "x64"),
         writable : false);
    if (key == null) {
      return false;
    } else {
      int? installed = (int?)key.GetValue("Installed");
      return installed == 1;
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
