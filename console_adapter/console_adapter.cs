using System;
using System.Runtime.InteropServices;

namespace princpia {
namespace console_adapter {

internal class ConsoleAdapter {

  [DllImport("test_plugin.dll", CallingConvention = CallingConvention.Cdecl)]
  private static extern int Say33();

  [DllImport("test_plugin.dll", CallingConvention = CallingConvention.Cdecl)]
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
}  // namespace princpia
