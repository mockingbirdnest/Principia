using System;
using System.IO;
using System.Net;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;
using static principia.tools.DbgHelp;

namespace principia {
namespace tools {

class StackTraceDecoder {
  // Returns the base address for the given DLL.
  private static Int64 GetBaseAddress(bool unity_crash,
                                      string unity_regex,
                                      string principia_cpp_regex,
                                      StreamReader stream) {
    var base_address_regex =
        new Regex(unity_crash ? unity_regex
                              : @"^I.*" + principia_cpp_regex +
                                    @":.*\] Base address is ([0-9A-F]+)$");
    Match base_address_match;
    do {
      base_address_match = base_address_regex.Match(stream.ReadLine());
    } while (!base_address_match.Success);
    string base_address_string = base_address_match.Groups[1].ToString();
    return Convert.ToInt64(base_address_string, 16);
  }

  // If the IMAGEHLP_LINEW64 represents a line of Principia code, writes a
  // GitHub link to the console and returns true.  Otherwise, returns false.
  private static bool ParseLine(IMAGEHLP_LINEW64 line, string commit) {
    var file_regex = new Regex(@".*\\principia\\([a-z_]+)\\(\S+)");
    Match file_match = file_regex.Match(line.FileName);
    if (file_match.Success) {
      string file = $"{file_match.Groups[1]}/{file_match.Groups[2]}";
      string url = "https://github.com/mockingbirdnest/Principia/blob/" +
                   $"{commit}/{file}#L{line.LineNumber}";
      Console.WriteLine($"[`{file}:{line.LineNumber}`]({url})");
      return true;
    } else {
      return false;
    }
  }

  private static void Win32Check(bool success,
                                 [CallerMemberName] string member = "",
                                 [CallerFilePath] string file = "",
                                 [CallerLineNumber] int line = 0) {
    if (!success) {
      Console.WriteLine($"Error {Marshal.GetLastWin32Error()}");
      Console.WriteLine($"{file}:{line} ({member})");
      Environment.Exit(1);
    }
  }

  private static void LogComment(string comment) {
    Console.WriteLine($"<!--- {comment} -->");
  }

  private static void Main(string[] args) {
    bool unity_crash;
    string commit = null;
    if (args.Length == 3) {
      unity_crash = false;
    } else if (args.Length == 4) {
      var match = Regex.Match(args[3],
                              "--unity-crash-at-commit=([0-9a-f]{40})");
      if (match.Success) {
        unity_crash = true;
        commit = match.Groups[1].ToString();
      } else {
        PrintUsage();
        return;
      }
    } else {
      PrintUsage();
      return;
    }
    string info_file_uri = args[0];
    string principia_pdb_file = args[1];
    string physics_pdb_file = args[2];
    var web_client = new WebClient();
    var stream = new StreamReader(web_client.OpenRead(info_file_uri),
                                  Encoding.UTF8);
    if (!unity_crash) {
      var version_regex = new Regex(
          @"^I.*\] Principia version " + 
          @"([0-9]{10}-\w+)-[0-9]+-g([0-9a-f]{40}) built");
      Match version_match;
      do {
        version_match = version_regex.Match(stream.ReadLine());
      } while (!version_match.Success);
      string tag = version_match.Groups[1].ToString();
      commit = version_match.Groups[2].ToString();
    }
    Int64 principia_base_address =
        GetBaseAddress(unity_crash,
                       @"GameData\\Principia\\principia.dll:principia.dll " +
                           @"\(([0-9A-F]+)\)",
                       "interface\\.cpp",
                       stream);
    Int64 physics_base_address =
        GetBaseAddress(unity_crash,
                       @"GameData\\Principia\\principia.dll:principia.dll " +
                           @"\(([0-9A-F]+)\)",
                       "ksp_physics_lib\\.cpp",
                       stream);
    LogComment($"Using Principia base address {principia_base_address:X}");
    LogComment($"Using Physics base address {physics_base_address:X}");
    var stack_regex = new Regex(
        unity_crash ? @"\(0x([0-9A-F]+)\) .*"
                    : @"@\s+[0-9A-F]+\s+.* \[0x([0-9A-F]+)(\+[0-9]+)?\]");
    Match stack_match;
    if (unity_crash) {
      Match stack_start_match;
      do {
        stack_start_match =
            Regex.Match(stream.ReadLine(),
                        @"========== OUTPUTING STACK TRACE ==================");
      } while (!stack_start_match.Success);
    }
    do {
      stack_match = stack_regex.Match(stream.ReadLine());
    } while (!stack_match.Success);
    IntPtr handle = new IntPtr(1729);
    SymSetOptions(SYMOPT_LOAD_LINES);
    Win32Check(SymInitializeW(handle, null, fInvadeProcess: false));
    Win32Check(SymLoadModuleExW(handle,
                                IntPtr.Zero,
                                principia_pdb_file.Replace(".pdb", ".dll"),
                                null,
                                principia_base_address,
                                0,
                                IntPtr.Zero,
                                0) != 0);
    Win32Check(SymLoadModuleExW(handle,
                                IntPtr.Zero,
                                physics_pdb_file.Replace(".pdb", ".dll"),
                                null,
                                physics_base_address,
                                0,
                                IntPtr.Zero,
                                0) != 0);

    for (;
         stack_match.Success;
         stack_match = stack_regex.Match(stream.ReadLine())) {
      Int64 address = Convert.ToInt64(stack_match.Groups[1].ToString(), 16);
      IMAGEHLP_LINEW64 line = new IMAGEHLP_LINEW64();
      if (SymGetLineFromAddrW64(handle,
                                address,
                                out Int32 displacement,
                                line)) {
        if (!ParseLine(line, commit)) {
           LogComment($"Not in Principia code: {stack_match.Groups[0]}");
        }
      } else if (Marshal.GetLastWin32Error() == 126) {
        LogComment($"Not in loaded modules: {stack_match.Groups[0]}");
      } else {
        Win32Check(false);
      }
      Int32 inline_trace = SymAddrIncludeInlineTrace(handle, address);
      if (inline_trace != 0) {
        LogComment($"{inline_trace} inline frames");
        Win32Check(SymQueryInlineTrace(handle,
                                       address,
                                       0,
                                       address,
                                       address,
                                       out Int32 current_context,
                                       out Int32 current_frame_index));
        for (int i = 0; i < inline_trace; ++i) {
          Win32Check(SymGetLineFromInlineContextW(
              handle, address, current_context + i, 0, out Int32 dsp, line));
          if (!ParseLine(line, commit)) {
            LogComment("Inline frame not in Principia code");
          }
        }
      }
    }
  }

  private static void PrintUsage() {
    Console.WriteLine("Usage: stacktrace_decoder " +
                      "<info_file_uri> <principia_pdb_file> " +
                      "<physics_pdb_file> [--unity-crash-at-commit=<sha1>]");
  }
}

}  // namespace tools
}  // namespace principia
