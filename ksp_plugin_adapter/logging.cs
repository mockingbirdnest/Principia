using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Log {

  [DllImport(dllName           : PluginAdapter.kDllPath,
             EntryPoint        = "principia__InitGoogleLogging",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void InitGoogleLogging();

  [DllImport(dllName           : PluginAdapter.kDllPath,
             EntryPoint        = "principia__LogInfo",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Info(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : PluginAdapter.kDllPath,
             EntryPoint        = "principia__LogWarning",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Warning(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : PluginAdapter.kDllPath,
             EntryPoint        = "principia__LogError",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Error(
      [MarshalAs(UnmanagedType.LPStr)] String message);

  [DllImport(dllName           : PluginAdapter.kDllPath,
             EntryPoint        = "principia__LogFatal",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void Fatal(
      [MarshalAs(UnmanagedType.LPStr)] String message);

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
