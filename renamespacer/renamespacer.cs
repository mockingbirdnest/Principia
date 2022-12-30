using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {
class Renamespacer {
  static void ProcessClientFile(FileInfo input_file,
                                string project_name,
                                bool move_files) {
    string input_filename = input_file.FullName;
    string output_filename = input_file.DirectoryName + "\\" +
                              input_file.Name + ".new";

    string project_namespace = project_name;

    var using_directives = new SortedSet<string>();
    bool has_emitted_new_style_usings = false;
    bool has_seen_usings = false;

    StreamReader reader = input_file.OpenText();
    StreamWriter writer = File.CreateText(output_filename);
    while (!reader.EndOfStream) {
      string line = reader.ReadLine();
      if (line.StartsWith("#include \"" + project_name) &&
          (line.EndsWith(".hpp\"") || line.EndsWith(".cpp\""))) {
        // Collect the includes for files in our project, we'll need them to
        // generate the using-directives.
        using_directives.Add(
            "principia::" + Regex.Replace(
                Regex.Replace(line, "#include \"", ""),
                "\\.[hc]pp\"$", "").Replace("/", "::"));
        writer.WriteLine(line);
      } else if (line.StartsWith("using ")) {
        // We are in the block of usings.  Preserve them unless they are for the
        // project that we are processing, or for an internal namespace.
        if (!line.StartsWith("using " + project_namespace) &&
            !line.StartsWith("using internal_")) {
          writer.WriteLine(line);
        } else if (line.StartsWith("using namespace")) {
          // Collect (and don't emit) existing using-directives as we'll want to
          // reemit them in order.
          using_directives.Add(
              line.Replace("using namespace ", "").Replace(";", ""));
          writer.WriteLine(line);
        }
        has_seen_usings = true;
      } else if (has_seen_usings && !has_emitted_new_style_usings) {
        // The first line that follows the using-declarations.  Emit the new
        // style using-directives here.
        foreach (string using_directive in using_directives) {
          writer.WriteLine("using namespace " + using_directive + ";");
        }
        writer.WriteLine(line);
        has_emitted_new_style_usings = true;
      } else {
        writer.WriteLine(line);
      }
    }
    writer.Close();
    reader.Close();
    if (move_files) {
      File.Move(output_filename, input_filename, true);
    }
  }

  static void ProcessProjectFile(FileInfo input_file,
                                 string project_name,
                                 bool move_files) {
    string input_filename = input_file.FullName;
    string output_filename = input_file.DirectoryName + "\\" +
                              input_file.Name + ".new";
    string file_basename = Regex.Replace(input_file.Name,  "\\.[hc]pp$", "");

    string project_namespace = project_name;
    string file_namespace = Regex.Replace(file_basename, "_body|_test", "");

    var using_directives = new SortedSet<string>();
    bool has_closed_file_namespace = false;
    bool has_closed_internal_namespace = false;
    bool has_emitted_new_style_usings = false;
    bool has_opened_file_namespace = false;
    bool has_seen_namespaces = false;
    bool has_seen_usings = false;

    StreamReader reader = input_file.OpenText();
    StreamWriter writer = File.CreateText(output_filename);
    while (!reader.EndOfStream) {
      string line = reader.ReadLine();
      if (line.StartsWith("#include \"" + project_name) &&
          (line.EndsWith(".hpp\"") || line.EndsWith(".cpp\""))) {
        // Collect the includes for files in our project, we'll need them to
        // generate the using-directives.
        using_directives.Add(
            "principia::" + Regex.Replace(
                Regex.Replace(line, "#include \"", ""),
                "\\.[hc]pp\"$", "").Replace("/", "::"));
        writer.WriteLine(line);
      } else if (line.StartsWith("namespace internal_")) {
        // The internal namespace gets wrapped into the file namespace and is
        // named "internal".
        writer.WriteLine("namespace " + file_namespace + " {");
        writer.WriteLine("namespace internal {");
        has_opened_file_namespace = true;
      } else if (line.StartsWith("namespace ")) {
        // Record that we have seen the first opening namespace.  Note that
        // this code assumes that all the opening namespaces are in sequence.
        writer.WriteLine(line);
        has_seen_namespaces = true;
      } else if (has_closed_internal_namespace &&
                  line.StartsWith("using internal_")) {
        // A using after the internal namespace has been closed, that exposes
        // the declaration to the outside.
        writer.WriteLine(Regex.Replace(line, "_" + file_namespace, ""));
      } else if (line.StartsWith("using ")) {
        // We are in the block of usings.  Preserve them unless they are for the
        // project that we are processing, or for an internal namespace.
        if (!line.StartsWith("using " + project_namespace) &&
            !line.StartsWith("using internal_")) {
          writer.WriteLine(line);
        } else if (line.StartsWith("using namespace")) {
          // Collect (and don't emit) existing using-directives as we'll want to
          // reemit them in order.
          using_directives.Add(
              line.Replace("using namespace ", "").Replace(";", ""));
          writer.WriteLine(line);
        }
        has_seen_usings = true;
      } else if (line.StartsWith("}  // namespace internal_")) {
        // Change the style of the line that closes the internal namespace.
        writer.WriteLine("}  // namespace internal");
        has_closed_internal_namespace = true;
      } else if (line.StartsWith("}  // namespace") &&
                  !has_closed_file_namespace) {
        // The close of a namespace, and we have not closed the file namespace
        // yet.  Do so now.
        writer.WriteLine("}  // namespace " + file_namespace);
        writer.WriteLine(line);
        has_closed_file_namespace = true;
      } else if ((has_seen_usings ||
                 (line != "" && has_opened_file_namespace)) &&
                 !has_emitted_new_style_usings) {
        // The first line that follows the using-declarations.  Emit the new
        // style using-directives here.  Note that this covers the first line
        // after the namespace that's not empty, in case there are no usings.
        foreach (string using_directive in using_directives) {
          writer.WriteLine("using namespace " + using_directive + ";");
        }
        if (!has_seen_usings) {
          writer.WriteLine("");
        }
        writer.WriteLine(line);
        has_emitted_new_style_usings = true;
      } else if (has_seen_namespaces && !has_opened_file_namespace) {
        // The first line that follows the opening namespaces.  Open the file
        // namespace if we haven't done so yet.
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
      File.Move(output_filename, input_filename, true);
    }
  }

  static void Main(string[] args) {
    // Parse the arguments.
    DirectoryInfo project = null;
    var clients = new List<DirectoryInfo>();
    bool move_files = false;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split = arg.Split(new []{"--", ":"}, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          project = new DirectoryInfo(value);
        } else if (option == "client") {
          clients.Add(new DirectoryInfo(value));
        }
      } else if (arg == "--move") {
        move_files = true;
      }
    }

    // Process the files in the project.
    string project_name = project.Name;
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
    foreach (FileInfo input_file in all_files) {
      ProcessProjectFile(input_file, project_name, move_files);
    }

    foreach (DirectoryInfo client in clients) {
      FileInfo[] client_hpp_files = client.GetFiles("*.hpp");
      FileInfo[] client_cpp_files = client.GetFiles("*.cpp");
      FileInfo[] all_client_files =
          client_hpp_files.Union(client_cpp_files).ToArray();
      foreach (FileInfo input_file in all_client_files) {
        ProcessClientFile(input_file, project_name, move_files);
      }
    }
  }
}

}  // namespace renamespacer
}  // namespace principia
