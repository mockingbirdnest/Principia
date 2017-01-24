using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading.Tasks;

namespace principia {
namespace parallel_test_runner {

class ParallelTestRunner {
  static void Main(string[] args) {
    var processes = new List<Process>();
    int test_process_counter = 0;
    foreach (string directory in args) {
      string[] test_binaries = Directory.GetFiles(directory, "*_tests.exe");
      foreach (string test_binary in test_binaries) {
        var info = new ProcessStartInfo(test_binary, "--gtest_list_tests");
        info.UseShellExecute = false;
        info.RedirectStandardOutput = true;
        var list_tests = Process.Start(info);
        var output = list_tests.StandardOutput;
        string test_case = null;
        while(!output.EndOfStream) {
          string line = output.ReadLine();
          if (line[0] != ' ') {
            test_case = line;
          } else {
            var process = new Process();
            process.StartInfo.UseShellExecute = false;
            process.StartInfo.RedirectStandardOutput = true;
            process.StartInfo.RedirectStandardError = true;
            process.StartInfo.FileName = test_binary;
            process.StartInfo.Arguments =
                "--gtest_filter=" + test_case + line.Split(' ')[2];
            if (process.StartInfo.Arguments ==
                "--gtest_filter=PlayerTest.Benchmarks") {
              continue;
            }
            process.StartInfo.Arguments +=
                " --gtest_output=xml:TestResults\\gtest_results_" +
                test_process_counter++ + ".xml";
            processes.Add(process);
          }
        }
      }
    }
    var stopwatch = new Stopwatch();
    stopwatch.Start();
    Task[] tasks = new Task[processes.Count];
    for (int i = 0; i < processes.Count; ++i) {
      var process = processes[i];
      tasks[i] = new Task(async () => {
        process.Start();
        while (!process.StandardOutput.EndOfStream) {
          string output = await process.StandardOutput.ReadLineAsync();
          if (output.StartsWith("[ ") && !output.StartsWith("[  PASSED  ]")) {
            Console.WriteLine(output);
          }
        }
        if (process.ExitCode != 0) {
          Console.WriteLine(process.StandardOutput.ReadToEnd());
          Console.WriteLine("Exit code " + process.ExitCode + " from " +
                            process.StartInfo.FileName + " " +
                            process.StartInfo.Arguments);
          Console.WriteLine("stderr follows:");
          Console.WriteLine(process.StandardError.ReadToEnd());
          Environment.Exit(process.ExitCode);
        }
      });
      tasks[i].Start();
    }

    Console.WriteLine(tasks.Length);
    Task.WaitAll(tasks);
    Console.WriteLine("Done (" + stopwatch.ElapsedMilliseconds + " ms)");
  }
}

}  // namespace parallel_test_runner
}  // namespace principia
