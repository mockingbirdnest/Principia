using System;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using Microsoft.Win32;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Loader {
  internal static string LoadPrincipiaDllAndInitGoogleLogging() {
    if (loaded_principia_dll) {
      return null;
    }
    bool is_32_bit = IntPtr.Size == 4;
    if (is_32_bit) {
      return "32-bit platforms are no longer supported; " +
             "use the 64-bit KSP executable.";
    }
    string[] possible_dll_paths;
    string dll_filename;
    bool? is_cxx_installed;
    string required_cxx_packages;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        is_cxx_installed = IsVCRedistInstalled();
        required_cxx_packages =
            "the Microsoft Visual C++ 2015-2022 Redistributable (x64) - " +
            "14.38.33130";
        dll_filename = "principia.dll";
        possible_dll_paths = new []
            { @"GameData\Principia\x64\" + dll_filename };
        break;
      // Both Mac and Linux report `PlatformID.Unix`, so we treat them together
      // (we probably don't actually encounter `PlatformID.MacOSX`).
      case PlatformID.Unix:
      case PlatformID.MacOSX:
        dll_filename = "principia.so";
        possible_dll_paths = new [] {
            @"GameData/Principia/Linux64/" + dll_filename,
            @"GameData/Principia/MacOS64/" + dll_filename,
        };
        is_cxx_installed = null;
        required_cxx_packages = "libc++abi1-17, libc++1-17, and libunwind-17 " +
                                "or later (Linux) or High Sierra or later " +
                                "(MacOS)";
        break;
      default:
        return "The operating system " +
               Environment.OSVersion +
               " is not supported at this time.";
    }
    if (!possible_dll_paths.Any(File.Exists)) {
      string[] where_did_they_put_the_dll = Directory.GetFiles(
          Directory.GetCurrentDirectory(),
          dll_filename,
          SearchOption.AllDirectories);
      string incorrectly_installed_in = "";
      if (where_did_they_put_the_dll.Any()) {
        incorrectly_installed_in = "  It was incorrectly installed in " +
                                   string.Join(", ",
                                               where_did_they_put_the_dll) +
                                   ".";
      }
      return "The principia DLL was not found at '" +
             string.Join("', '", possible_dll_paths) +
             "' in directory '" +
             Directory.GetCurrentDirectory() +
             "'." +
             incorrectly_installed_in;
    }
    string non_ascii_path_error = null;
    foreach (char c in Directory.GetCurrentDirectory()) {
      if (c >= 128) {
        non_ascii_path_error = Directory.GetCurrentDirectory() +
                               " contains the non-ASCII character " +
                               c +
                               "; this is known to confuse Mono.";
        break;
      }
    }
    try {
      loaded_principia_dll = true;
      Log.InitGoogleLogging();
      return null;
    } catch (Exception e) {
      UnityEngine.Debug.LogException(e);
      if (non_ascii_path_error != null) {
        return non_ascii_path_error;
      } else if (is_cxx_installed == false) {
        return "Dependencies, namely " +
               required_cxx_packages +
               ", were not found.";
      } else {
        return "An unknown error occurred; detected OS " +
               Environment.OSVersion +
               " 64-bit; tried loading dll at '" +
               string.Join("', '", possible_dll_paths) +
               "'. Note that " +
               required_cxx_packages +
               " are required.";
      }
    }
  }

  private static bool IsVCRedistInstalled() {
    // NOTE(phl): This GUID is specific to:
    //   Microsoft Visual C++ 2015-2022 Redistributable (x64) - 14.38.33130
    // It will need to be updated when new versions of Visual C++
    // Redistributable are released by Microsoft.
    RegistryKey key = Registry.LocalMachine.OpenSubKey(
        @"Software\Classes\Installer\Dependencies\" +
        @"VC,redist.x64,amd64,14.38,bundle",
        writable : false);
    if (key == null) {
      return false;
    } else {
      string version = (string)key.GetValue("Version");
      // NOTE(phl): This string needs to be updated when new versions of Visual
      // C++ Redistributable are released by Microsoft.
      return version != null && version == "14.38.33130.0";
    }
  }

  [DllImport("kernel32", CharSet = CharSet.Unicode, SetLastError = true)]
  [return: MarshalAs(UnmanagedType.Bool)]
  static extern bool SetDllDirectory(string lpPathName);

  [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
  private static extern IntPtr LoadLibrary(
      [MarshalAs(UnmanagedType.LPStr)] string lpFileName);

  internal static bool loaded_principia_dll { get; private set; } = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
