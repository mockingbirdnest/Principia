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
    DirectoryInfo benchmark_directory = new DirectoryInfo(args[0]);
    DirectoryInfo mathematica_directory = new DirectoryInfo(args[1]);
    DirectoryInfo jenkins_directory = new DirectoryInfo(args[2]);
    DateTime date = DateTime.UtcNow;
    String mathematica_date = date.ToString("{yyyy, M, d, H, m, s.fffffff}");
    String mathematica_output_file =
        Path.Combine(mathematica_directory.FullName,
                     "principia_benchmark_results_" +
                         date.ToString("yyyy_M_d__H_m_s_fffffff.wl"));
    Console.WriteLine("DateList[TimeZone->0]: " + mathematica_date);
    Console.WriteLine("Results will be written to " + mathematica_output_file);
    StreamWriter mathematica_stream =
        new StreamWriter(
            File.Open(mathematica_output_file, FileMode.CreateNew));
    String csv_benchmark_names = "";
    String csv_means = "";
    mathematica_stream.WriteLine("{");
    mathematica_stream.WriteLine(mathematica_date);
    FileInfo[] files = benchmark_directory.GetFiles("*.cpp");
    foreach (FileInfo file in files) {
      StreamReader stream = file.OpenText();
      String command_line;
      do {
        command_line = stream.ReadLine();
      } while (command_line == "");
      if (command_line.StartsWith("// " + benchmark_executable)) {
        // Get rid of the // NOLINT comments and of the actual command, leaving
        // only the arguments.
        command_line =
            command_line.Split(
                separator : new String[]{benchmark_executable, "//"},
                options   : StringSplitOptions.None)[2].Trim();
        bool has_repetitions = command_line.Contains("--benchmark_repetitions");
        Console.WriteLine(
            "Running benchmarks with arguments from " + file.Name);
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
          String[] words =
              line.Split(separator : new Char[]{' '},
                         options   : StringSplitOptions.RemoveEmptyEntries);
          if (has_repetitions) {
            const String mean_postfix = "_mean";
            const String stddev_postfix = "_stddev";
            if (words[0].EndsWith(mean_postfix)) {
              String benchmark_name =
                  words[0].Substring(
                      startIndex : 0,
                      length     : words[0].Length - mean_postfix.Length);
              Int64 μ = Int64.Parse(words[1]);
              Console.WriteLine(benchmark_name + ": μ = " + μ + " ns");
              mathematica_stream.WriteLine(",");
              mathematica_stream.WriteLine("{");
              mathematica_stream.WriteLine("\"" + benchmark_name + "\",");
              mathematica_stream.WriteLine(μ + ",");
              CommaSeparatedAppend(
                  ref csv_benchmark_names,
                  "\"" + benchmark_name.Replace("\"", "\"\"") + "\"");
              CommaSeparatedAppend(ref csv_means, μ.ToString());
            } else if (words[0].EndsWith(stddev_postfix)) {
              String benchmark_name =
                  words[0].Substring(
                      startIndex : 0,
                      length     : words[0].Length - stddev_postfix.Length);
              Int64 σ = Int64.Parse(words[1]);
              Console.WriteLine(benchmark_name + ": σ = " + σ + " ns");
              mathematica_stream.WriteLine(σ);
              mathematica_stream.WriteLine("}");
              mathematica_stream.Flush();
            }
          } else if (words[0].StartsWith("BM")) {
            String benchmark_name = words[0];
            Int64 μ = Int64.Parse(words[1]);
            Int64 σ = 0;
            CommaSeparatedAppend(
                ref csv_benchmark_names,
                "\"" + benchmark_name.Replace("\"", "\"\"") + "\"");
            CommaSeparatedAppend(ref csv_means, μ.ToString());
            Console.WriteLine(benchmark_name + ": μ = " + μ + " ns");
            Console.WriteLine(benchmark_name + ": σ = " + σ + " ns");
            mathematica_stream.WriteLine(",");
            mathematica_stream.WriteLine("{");
            mathematica_stream.WriteLine("\"" + words[0] + "\",");
            mathematica_stream.WriteLine(μ + ",");
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
    File.WriteAllText(
        Path.Combine(jenkins_directory.FullName,
                     "jenkins_benchmark_results.csv"),
        csv_benchmark_names + "\n" + csv_means);
  }

  private static void CommaSeparatedAppend(ref string csv, string value) {
    if (csv == "") {
      csv = value;
    } else {
      csv = csv + ", " + value;
    }
  }
}

}  // namespace benchmark_automation
}  // namespace principia
