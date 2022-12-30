using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {
class Renamespacer {
  static void Main(string[] args) {
    // Parse the arguments.
    DirectoryInfo directory = null;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split = arg.Split(new []{"--", ":"}, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "directory") {
          directory = new DirectoryInfo(value);
        }
      }
    }

    // Find the files to process.
    FileInfo[] hpp_files = directory.GetFiles("*.hpp");
    FileInfo[] cpp_files = directory.GetFiles("*.cpp");
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();

    foreach (FileInfo input_file in all_files) {
      string basename = Regex.Replace(input_file.Name,  "\\.[hc]pp$", "");
      StreamWriter writer = File.CreateText(input_file.DirectoryName + "\\" +
                                            input_file.Name + ".new");
      var included_files = new List<string>();
      bool has_emitted_usings = false;

      StreamReader reader = input_file.OpenText();
      while (!reader.EndOfStream) {
        string line = reader.ReadLine();
        if (line.StartsWith("#include")) {
          included_files.Add(
              Regex.Replace(
                  Regex.Replace(line, "#include.*/", ""),
                  "\\.[hc]pp$", ""));
          writer.WriteLine(line);
        } else if (line.StartsWith("namespace internal_")) {
          writer.WriteLine("namespace " + basename + " {");
          writer.WriteLine("namespace internal {");
        } else if (line.StartsWith("using internal_")) {
          writer.WriteLine(Regex.Replace(line, "_" + basename, ""));
        } else if (line.StartsWith("using ")) {
          if (!has_emitted_usings) {
            foreach (string included_file in included_files) {
              writer.WriteLine("using namespace principia::" + basename + ";");
            }
            has_emitted_usings = true;
          }
        } else if (line.StartsWith("}  // namespace internal_")) {
          writer.WriteLine("}  // namespace internal");
        } else {
          writer.WriteLine(line);
        }
      }
    }
  }
}

}  // namespace renamespacer
}  // namespace principia
