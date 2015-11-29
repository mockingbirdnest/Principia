using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Log {
  internal static String[] kSeverityNames = {"INFO", "WARNING", "ERROR",
                                             "FATAL"};

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__InitGoogleLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void InitGoogleLogging();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__SetBufferedLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void SetBufferedLogging(int max_severity);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__GetBufferedLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern int GetBufferedLogging();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__SetBufferDuration",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void SetBufferDuration(int seconds);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__GetBufferDuration",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern int GetBufferDuration();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__SetSuppressedLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void SetSuppressedLogging(int min_severity);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__GetSuppressedLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern int GetSuppressedLogging();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__SetVerboseLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void SetVerboseLogging(int level);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__GetVerboseLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern int GetVerboseLogging();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__SetStderrLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void SetStderrLogging(int min_severity);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__GetStderrLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern int GetStderrLogging();

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__LogInfo",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Info(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__LogWarning",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Warning(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__LogError",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Error(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : Interface.kDllPath,
             EntryPoint        = "principia__LogFatal",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Fatal(
      [MarshalAs(UnmanagedType.LPStr)] String message);
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
