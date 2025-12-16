using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
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
    string[] possible_dll_path_patterns;
    string dll_filename;
    bool? is_cxx_installed;
    string required_cxx_packages;
    switch (Environment.OSVersion.Platform) {
      case PlatformID.Win32NT:
        is_cxx_installed = IsVCRedistInstalled();
        required_cxx_packages =
            "the Microsoft Visual C++ 2015-2022 Redistributable (x64) - " +
            vc_redist_bundle_ + "." + vc_redist_version_;
        dll_filename = "principia.dll";
        possible_dll_path_patterns = new []
            { @"GameData\Principia\Windows\{0}\" + dll_filename };
        break;
      // Both Mac and Linux report `PlatformID.Unix`, so we treat them together
      // (we probably don't actually encounter `PlatformID.MacOSX`).
      case PlatformID.Unix:
      case PlatformID.MacOSX:
        dll_filename = "principia.so";
        possible_dll_path_patterns = new [] {
            @"GameData/Principia/Linux/{0}/" + dll_filename,
            @"GameData/Principia/macOS/{0}/" + dll_filename,
        };
        is_cxx_installed = null;
        required_cxx_packages = "libc++abi1-20, libc++1-20, and libunwind-20 " +
                                "or later (Linux) or High Sierra or later " +
                                "(MacOS)";
        break;
      default:
        return "The operating system " +
               Environment.OSVersion +
               " is not supported at this time.";
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
      // First try to load the x64 DLL.  This should always work, and it makes
      // it possible to query the CPU flags.  Note that we can only call a
      // limited set of function before we load the final DLL.  In particular,
      // we cannot `InitGoogleLogging` as this seems to reference dangling
      // pointers on the second load.
      string error_x64 = LoadPrincipiaDllForPlatform(
          possible_dll_path_patterns,
          dll_filename: dll_filename,
          platform: "x64");
      if (error_x64 != null) {
        return error_x64;
      }
      Interface.GetCPUIDFeatureFlags(out bool has_avx, out bool has_fma);
      if (has_fma && !has_avx) {
        return "Principia does not run on processors with FMA support but no " +
               "AVX support.";
      } else if (has_fma) {
        // If it turns out that the machine has FMA (and therefore AVX), change
        // our mind and load the x64_AVX_FMA DLL.  It's necessary to unload the
        // DLL, otherwise the base address reported by `InitGoogleLogging` might
        // (incorrectly) be the one of the `x64` DLL.
        UnloadPrincipiaDll(platform: "x64");
        string error_x64_avx_fma = LoadPrincipiaDllForPlatform(
            possible_dll_path_patterns,
            dll_filename: dll_filename,
            platform: "x64_AVX_FMA");
        if (error_x64_avx_fma != null) {
          return error_x64_avx_fma;
        }
      }
      loaded_principia_dll = true;
      Log.InitGoogleLogging();
      Log.Info("Processor " +
               (has_avx ? "has" : "does not have") +
               " AVX support and " +
               (has_fma ? "has" : "does not have") +
               " FMA support.");
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
               string.Join("', '", possible_dll_path_patterns) +
               "'. Note that " +
               required_cxx_packages +
               " are required.";
      }
    }
  }

  private static string LoadPrincipiaDllForPlatform(
      string[] possible_dll_path_patterns,
      string dll_filename,
      string platform) {
    var possible_dll_paths = new List<string>();
    foreach (string pattern in possible_dll_path_patterns) {
      possible_dll_paths.Add(
          string.Format(CultureInfo.InvariantCulture, pattern, platform));
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
             string.Join("', '", possible_dll_paths.ToArray()) +
             "' in directory '" +
             Directory.GetCurrentDirectory() +
             "'." +
             incorrectly_installed_in;
    }

    UnityEngine.Debug.Log("Loading the " +
                          platform +
                          " Principia DLL " +
                          dll_filename +
                          " from paths: " +
                          string.Join(", ", possible_dll_paths));
    if (Environment.OSVersion.Platform == PlatformID.Win32NT) {
      principia_dll_ = LoadLibrary(possible_dll_paths[0]);
    } else {
      const int RTLD_NOW = 2;
      const int RTLD_GLOBAL = 8;
      foreach (string path in possible_dll_paths) {
        principia_dll_ = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
        if (principia_dll_ != IntPtr.Zero) {
          break;
        }
      }
    }
    UnityEngine.Debug.Log("Principia DLL loaded at address " +
                          principia_dll_.ToString("X16"));
    Interface.LoadSymbols();
    return null;
  }

  private static string UnloadPrincipiaDll(string platform) {
    if (Environment.OSVersion.Platform == PlatformID.Win32NT) {
      if (!FreeLibrary(principia_dll_)) {
        return "Unable to free the " + platform + "  Principia DLL.";
      }
    } else {
      if (dlclose(principia_dll_) != 0) {
        return "Unable to close the " + platform + " Principia DLL";
      }
    }
    return null;
  }

  private static bool IsVCRedistInstalled() {
    RegistryKey key = Registry.LocalMachine.OpenSubKey(
        @"Software\Classes\Installer\Dependencies\" +
        @"VC,redist.x64,amd64," +
        vc_redist_bundle_ +
        ",bundle",
        writable : false);
    if (key == null) {
      return false;
    } else {
      string version = (string)key.GetValue("Version");
      return version != null &&
             version == vc_redist_bundle_ + "." + vc_redist_version_ + ".0";
    }
  }

  private static IntPtr principia_dll_;

  internal static T LoadFunction<T>(string function) where T : Delegate {
    IntPtr function_pointer;
    if (Environment.OSVersion.Platform == PlatformID.Win32NT) {
      function_pointer = GetProcAddress(principia_dll_, function);
    } else {
      function_pointer = dlsym(principia_dll_, function);
    }
    T result = Marshal.GetDelegateForFunctionPointer(
        function_pointer, typeof(T)) as T;
    if (result == null) {
      throw new EntryPointNotFoundException(function);
    }
    return result;
  }

  [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
  private static extern IntPtr LoadLibrary(
      [MarshalAs(UnmanagedType.LPStr)] string lpFileName);

  [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
  private static extern bool FreeLibrary(IntPtr hLibModule);

  [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
  private static extern IntPtr GetProcAddress(
      IntPtr hModule, [MarshalAs(UnmanagedType.LPStr)] string lpProcName);

  [DllImport("dl")]
  private static extern int dlclose(IntPtr handle);

  [DllImport("dl")]
  private static extern IntPtr dlopen(string filename, int flags);

  [DllImport("dl")]
  private static extern IntPtr dlsym(IntPtr handle, string symbol);

  // NOTE(phl): These strings need to be updated when new versions of Visual
  // C++ Redistributable are released by Microsoft.
  private static readonly string vc_redist_bundle_ = "14.44";
  private static readonly string vc_redist_version_ = "35211";

  internal static bool loaded_principia_dll { get; private set; } = false;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
