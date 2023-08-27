using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.IO;
using System.Linq;
using static principia.sourcerer.Filenames;

namespace principia {
namespace sourcerer {

class MakeEditorconfig {
  public static void Run(string[] args) {
    var projects = new List<DirectoryInfo>();
    var extras = new Dictionary<string, string>();
    string? solution = null;
    bool dry_run = true;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split =
            arg.Split(new []{ "--", ":" }, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          projects.Add(new DirectoryInfo(value));
        } else if (option == "solution") {
          solution = Path.GetFullPath(value);
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        } else if (option == "extra") {
          // --extra:gtest/gtest.h=gtest/gtest-matchers.h+gtest/internal/gtest-internal.h
          string[] facade_and_implementations = value.Split('=');
          string facade = facade_and_implementations[0];
          string[] implementations = facade_and_implementations[1].Split('+');
          foreach (string implementation in implementations) {
            extras.Add(implementation, facade);
          }
        } else {
          throw new ArgumentException("Unknown option " + option);
        }
      }
    }

    var corresponding_header = new Dictionary<FileInfo, FileInfo>();
    foreach (DirectoryInfo project in projects) {
      FileInfo[] hpp_files = project.GetFiles("*.hpp");
      FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);

      var file_name_to_file_info = new Dictionary<string, FileInfo>();
      foreach (FileInfo input_file in hpp_files) {
        file_name_to_file_info.Add(Path.GetFullPath(input_file.Name,
                                                    input_file.DirectoryName!),
                                   input_file);
      }

      foreach (FileInfo input_file in body_hpp_files) {
        string corresponding_header_file_name =
            Filenames.CorrespondingHeader(input_file);
        // Benchmarks and tests don't have a corresponding header.
        if (file_name_to_file_info.ContainsKey(
                corresponding_header_file_name)) {
          var corresponding_header_file =
              file_name_to_file_info[corresponding_header_file_name];
          if (corresponding_header_file != input_file) {
            corresponding_header.Add(input_file, corresponding_header_file);
          }
        }
      }
    }

    string cpp_include_cleanup_alternate_files =
        "cpp_include_cleanup_alternate_files";

    // Generate the alternate files map.
    string input_filename = solution! + "/.editorconfig";
    string output_filename = input_filename + ".new";
    var writer = new StreamWriter(output_filename);
    using (var reader = new StreamReader(input_filename)) {
      while (!reader.EndOfStream) {
        string line = reader.ReadLine()!;
        if (line.StartsWith(cpp_include_cleanup_alternate_files)) {
          // The alternate files map is overwritten.
          bool first = true;
          foreach (var correspondence in corresponding_header) {
            var body = Filenames.ToSlash(
                Filenames.SolutionRelativePath(solution!, correspondence.Key));
            var header = Filenames.ToSlash(
                Filenames.SolutionRelativePath(solution!,
                                               correspondence.Value));
            if (first) {
              writer.Write("cpp_include_cleanup_alternate_files = ");
              first = false;
            } else {
              writer.Write(", ");
            }
            writer.Write(header + ":" + body);
          }
          foreach (var extra in extras) {
            var body = extra.Key;
            var header = extra.Value;
            writer.Write(", " + header + ":" + body);
          }
          writer.Write("\r\n");
        } else {
          // Uninteresting lines are copied verbatim to the new file.
          writer.WriteLine(line);
        }
      }
      writer.Close();
    }

    if (!dry_run) {
      System.IO.File.Move(output_filename, input_filename, overwrite: true);
    }
  }
}

} // namespace sourcerer
} // namespace principia
