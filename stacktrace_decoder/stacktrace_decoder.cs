using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace principia {
namespace tools {

class StackTraceDecoder {
  const string kDBH =
      @"\Program Files (x86)\Windows Kits\8.1\Debuggers\x86\dbh.exe";

  private static void Main(string[] args) {
    if (args.Length != 2) {
      Console.WriteLine("Usage: stacktrace_decoder <info_file_uri> <pdb_file>");
      return;
    }
    string info_file_uri = args[0];
    string pdb_file = args[1];
    var web_client = new WebClient();
    var stream = new StreamReader(web_client.OpenRead(info_file_uri));
    var version_regex = new Regex(
        @"^I.*\] Principia version " + 
        @"([0-9]{10}-[A-Za-z]+)-[0-9]+-g([0-9a-f]{40}) built");
    Match version_match;
    do {
      version_match = version_regex.Match(stream.ReadLine());
    } while (!version_match.Success);
    string tag = version_match.Groups[1].ToString();
    string commit = version_match.Groups[2].ToString();
    var base_address_regex = new Regex(@"^I.*\] Base address is ([0-9A-F]+)$");
    string base_address_string =
        base_address_regex.Match(stream.ReadLine()).Groups[1].ToString();
    Int64 base_address = Convert.ToInt64(base_address_string, 16);
    var stack_regex = new Regex(
        @"^\s+@\s+[0-9A-F]+\s+\(No symbol\) \[0x([0-9A-F]+)\]");
    Match stack_match;
    do {
      stack_match = stack_regex.Match(stream.ReadLine());
    } while (!stack_match.Success);
    var file_regex = new Regex(
        @"file\s+:\s+.*\\principia\\([a-z_]+)\\([a-z_]+\.[ch]pp)");
    var line_regex = new Regex(@"line\s+:\s+([0-9]+)");
    for (;
         stack_match.Success;
         stack_match = stack_regex.Match(stream.ReadLine())) {
      Int64 address = Convert.ToInt64(stack_match.Groups[1].ToString(), 16);
      Int64 dbh_base_address = 0x1000000;
      string rebased_address =
          Convert.ToString(address - base_address + dbh_base_address, 16);
      var p = new Process();
      p.StartInfo.UseShellExecute = false;
      p.StartInfo.RedirectStandardOutput = true;
      p.StartInfo.FileName = kDBH;
      p.StartInfo.Arguments =
          '"' + pdb_file + "\" laddr \"" + rebased_address + '"';
      p.Start();
      string output = p.StandardOutput.ReadToEnd();
      p.WaitForExit();
      Match file_match = file_regex.Match(output);
      if (file_match.Success) {
        string file = file_match.Groups[1].ToString() + '/' +
                      file_match.Groups[2].ToString();
        string line = line_regex.Match(output).Groups[1].ToString();
        string url = "https://github.com/mockingbirdnest/Principia/blob/" +
                     commit + '/' + file + "#L" + line;
        Console.WriteLine("[`" + file + ":" + line + "`](" + url + ")");
      }
    }
  }
}

}  // namespace tools
}  // namespace principia
