using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;


namespace principia {
namespace tools {

internal static class DbgHelp {
  internal const UInt32 SYMOPT_LOAD_LINES = 0x00000010;
  internal const Int32 MAX_SYM_NAME = 2000;

  [StructLayout(LayoutKind.Sequential,
                Size = 88 + 2 * (MAX_SYM_NAME - 1),
                CharSet = CharSet.Unicode)]
  internal class SYMBOL_INFOW {
    public Int32 SizeOfStruct = 88;
    public Int32 TypeIndex;
    public Int64 Reserved0;
    public Int64 Reserved1;
    public Int32 Index;
    public Int32 Size;
    public Int64 ModBase;
    public Int32 Flags;
    public Int64 Value;
    public Int64 Address;
    public Int32 Register;
    public Int32 Scope;
    public Int32 Tag;
    public Int32 NameLen;
    public Int32 MaxNameLen = MAX_SYM_NAME;
    // The entire name goes here, but in order to access it we would need
    // a fixed-size buffer, which would require making this an unsafe struct
    // instead of a class.  We don't need the name at this point.
    public char Name0;
  }

  [StructLayout(LayoutKind.Sequential)]
  internal class IMAGEHLP_LINEW64 {
    public Int32 SizeOfStruct = Marshal.SizeOf<IMAGEHLP_LINEW64>();
    public IntPtr Key;
    public Int32 LineNumber;
    private IntPtr FileName_;
    public Int64 Address;

    public string FileName {
      get {
        StringBuilder result = new StringBuilder();
        for (int i = 0;; i += 2) {
          char code_unit = (char)Marshal.ReadInt16(FileName_, i);
          if (code_unit == 0) {
            return result.ToString();
          }
          result.Append(code_unit);
        }
      }
    }
  }

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymInitializeW(
      IntPtr hProcess,
      [MarshalAs(UnmanagedType.LPWStr)] string UserSearchPath,
      [MarshalAs(UnmanagedType.Bool)] bool fInvadeProcess);

  [DllImport("dbghelp.dll", CharSet = CharSet.Unicode)]
  internal static extern Int32 SymAddrIncludeInlineTrace(
      IntPtr hProcess,
      Int64 Address);

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymQueryInlineTrace(
      IntPtr hProcess,
      Int64 StartAddress,
      Int32 StartContext,
      Int64 StartRetAddress,
      Int64 CurAddress,
      out Int32 CurContext,
      out Int32 CurFrameIndex);

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymGetLineFromInlineContextW(
      IntPtr hProcess,
      Int64 dwAddr,
      Int32 InlineContext,
      Int64 qwModuleBaseAddress,
      out Int32 pdwDisplacement,
      [MarshalAs(UnmanagedType.LPStruct)] IMAGEHLP_LINEW64 Line);


  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymGetLineFromAddrW64(
      IntPtr hProcess,
      Int64 dwAddr,
      out Int32 pdwDisplacement,
      [MarshalAs(UnmanagedType.LPStruct)] IMAGEHLP_LINEW64 Line);

  [DllImport("dbghelp.dll", CharSet = CharSet.Unicode)]
  internal static extern UInt32 SymSetOptions(
      UInt32 SymOptions);

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  internal static extern Int64 SymLoadModuleExW(
      IntPtr hProcess,
      IntPtr hFile,
      [MarshalAs(UnmanagedType.LPWStr)] string ImageName,
      [MarshalAs(UnmanagedType.LPWStr)] string ModuleName,
      Int64 BaseOfDll,
      Int32 DllSize,
      IntPtr Data,
      UInt32 Flags);

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymFromAddrW(
      IntPtr hProcess,
      Int64 Address,
      out Int64 Displacement,
      [MarshalAs(UnmanagedType.LPStruct)] SYMBOL_INFOW Symbol);

  [DllImport("dbghelp.dll", SetLastError = true, CharSet = CharSet.Unicode)]
  [return: MarshalAs(UnmanagedType.Bool)]
  internal static extern bool SymFromInlineContextW(
      IntPtr hProcess,
      Int64 Address,
      Int32 InlineContext,
      out Int64 Displacement,
      [MarshalAs(UnmanagedType.LPStruct)] SYMBOL_INFOW Symbol);
}

} // namespace tools
} // namespace principia
