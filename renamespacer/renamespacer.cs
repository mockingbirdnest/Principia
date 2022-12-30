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
    bool move_files = false;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split = arg.Split(new []{"--", ":"}, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "directory") {
          directory = new DirectoryInfo(value);
        }
      } else if (arg == "--move") {
        move_files = true;
      }
    }

    // Find the files to process.
    FileInfo[] hpp_files = directory.GetFiles("*.hpp");
    FileInfo[] cpp_files = directory.GetFiles("*.cpp");
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();

    foreach (FileInfo input_file in all_files) {
      string input_filename = input_file.FullName;
      string output_filename = input_file.DirectoryName + "\\" +
                               input_file.Name + ".new";
      string basename = Regex.Replace(input_file.Name,  "\\.[hc]pp$", "");
      string file_namespace = Regex.Replace(basename, "_body", "");
      StreamWriter writer = File.CreateText(output_filename);
      var included_files = new List<string>();
      bool has_closed_file_namespace = false;
      bool has_closed_internal_namespace = false;
      bool has_emitted_new_style_usings = false;
      bool has_opened_file_namespace = false;
      bool has_seen_namespaces = false;
      bool has_seen_usings = false;

      StreamReader reader = input_file.OpenText();
      while (!reader.EndOfStream) {
        string line = reader.ReadLine();
        if (line.StartsWith("#include \"") &&
            (line.EndsWith(".hpp\"") || line.EndsWith(".cpp\""))) {
          included_files.Add(
              Regex.Replace(
                  Regex.Replace(line, "#include \"", ""),
                  "\\.[hc]pp\"$", ""));
          writer.WriteLine(line);
        } else if (line.StartsWith("namespace internal_")) {
          writer.WriteLine("namespace " + file_namespace + " {");
          writer.WriteLine("namespace internal {");
          has_opened_file_namespace = true;
        } else if (line.StartsWith("namespace ")) {
          writer.WriteLine(line);
          has_seen_namespaces = true;
        } else if (has_closed_internal_namespace &&
                   line.StartsWith("using internal_")) {
          writer.WriteLine(Regex.Replace(line, "_" + file_namespace, ""));
        } else if (line.StartsWith("using ")) {
          if (!line.StartsWith("using " + file_namespace + "::")) {
            writer.WriteLine(line);
          }
          has_seen_usings = true;
        } else if (line.StartsWith("}  // namespace internal_")) {
          writer.WriteLine("}  // namespace internal");
          has_closed_internal_namespace = true;
        } else if (line.StartsWith("}  // namespace") &&
                   !has_closed_file_namespace) {
          writer.WriteLine("}  // namespace " + file_namespace);
          writer.WriteLine(line);
          has_closed_file_namespace = true;
        } else if (has_seen_usings && !has_emitted_new_style_usings) {
          foreach (string included_file in included_files) {
            writer.WriteLine("using namespace principia::" +
                              included_file.Replace("/", "::") + ";");
          }
          writer.WriteLine(line);
          has_emitted_new_style_usings = true;
        } else if (has_seen_namespaces && !has_opened_file_namespace) {
          writer.WriteLine("namespace " + file_namespace + " {");
          writer.WriteLine(line);
          has_opened_file_namespace = true;
        } else {
          writer.WriteLine(line);
        }
      }
      writer.Close();
      reader.Close();
      if (move_files) {
        File.Move(output_filename, input_filename);
      }
    }
  }
}

}  // namespace renamespacer
}  // namespace principia
