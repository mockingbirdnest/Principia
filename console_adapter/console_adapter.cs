using System;
using System.Runtime.InteropServices;

namespace principia {
namespace console_adapter {

internal class ConsoleAdapter {
  private const string kDllPath = "GameData/Principia/principia.dll";

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern int Say33();

  [DllImport(dllName           : kDllPath,
             CallingConvention = CallingConvention.Cdecl)]
  private static extern IntPtr SayHello();

  internal static void Main(string[] args) {
    Console.WriteLine(Say33());
    IntPtr hello_ptr = SayHello();
    string hello = Marshal.PtrToStringAnsi(hello_ptr);
    Console.WriteLine(hello);
    Console.ReadLine();
  }
}

}  // namespace console_adapter
}  // namespace principia
