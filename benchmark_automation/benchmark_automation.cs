using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace benchmark_automation {

class BenchmarkAutomation {
  private const String benchmark_executable = @".\Release\benchmarks.exe";

  private static void Main(string[] args) {
    DirectoryInfo directory = new DirectoryInfo(args[0]);
    DateTime date = DateTime.UtcNow;
    String mathematica_date = date.ToString("{yyyy, M, d, H, m, s.fffffff}");
    String mathematica_output_file =
        Path.Combine(args[1],
                     "principia_benchmark_results_" +
                         date.ToString("yyyy_M_d__H_m_s_fffffff.wl"));
    Console.WriteLine("DateList[TimeZone->0]: " + mathematica_date);
    Console.WriteLine("Results will be written to " + mathematica_output_file);
    StreamWriter mathematica_stream =
        new StreamWriter(
            File.Open(mathematica_output_file, FileMode.CreateNew));
    mathematica_stream.WriteLine("{");
    mathematica_stream.WriteLine(mathematica_date);
    FileInfo[] files = directory.GetFiles("*.cpp");
    foreach (FileInfo file in files) {
      StreamReader stream = file.OpenText();
      String command_line;
      do {
        command_line = stream.ReadLine();
      } while (command_line == "");
      if (command_line.StartsWith(@"// " + benchmark_executable)) {
        // Get rid of the // NOLINT comments and of the actual command, leaving
        // only the arguments.
        command_line =
            command_line.Split(
                separator : new String[]{benchmark_executable, "//"},
                options   : StringSplitOptions.None)[2].Trim();
        Console.WriteLine(command_line);
        Process process = new Process {
          StartInfo = new ProcessStartInfo {
            FileName = benchmark_executable,
            Arguments = command_line,
            UseShellExecute = false,
            RedirectStandardOutput = true,
            CreateNoWindow = true
          }
        };
        process.Start();
        while (!process.StandardOutput.EndOfStream) {
          String line = process.StandardOutput.ReadLine();
          Console.WriteLine(line);
          String[] words =
              line.Split(separator : new Char[]{' '},
                         options   : StringSplitOptions.RemoveEmptyEntries);
          if (words[0].EndsWith("_mean")) {
            String benchmark_name = words[0].Substring(0, words[0].Length - 5);
            Int64 μ = Int64.Parse(words[1]);
            Console.WriteLine(benchmark_name + ": μ = " + μ + " ns");
            mathematica_stream.WriteLine(",");
            mathematica_stream.WriteLine("{");
            mathematica_stream.WriteLine("\"" + benchmark_name + "\",");
            mathematica_stream.WriteLine(μ + ",");
          } else if (words[0].EndsWith("_stddev")) {
            String benchmark_name = words[0].Substring(0, words[0].Length - 7);
            Int64 σ = Int64.Parse(words[1]);
            Console.WriteLine(benchmark_name + ": σ = " + σ + " ns");
            mathematica_stream.WriteLine(σ);
            mathematica_stream.WriteLine("}");
            mathematica_stream.Flush();
          }
        }
      } else {
        Console.WriteLine("No benchmark command found at the beginning of " +
                          file.Name);
      }
    }
    mathematica_stream.WriteLine("}");
    mathematica_stream.Close();
  }
}

}  // namespace benchmark_automation
}  // namespace principia
