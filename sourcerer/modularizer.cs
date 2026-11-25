using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;
using static principia.sourcerer.Rewriter;

namespace principia {
namespace sourcerer {

class Modularizer {

  public static void Run(string[] args) {
    // Parse the arguments.
    DirectoryInfo? project = null;
    var clients = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
    var extra_headers = new HashSet<string>();
    bool dry_run = true;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split =
            arg.Split(new []{ "--", ":" }, StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          project = new DirectoryInfo(value);
        } else if (option == "client") {
          clients.Add(new DirectoryInfo(value));
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        } else if (option == "exclude") {
          excluded.Add(value);
        } else if (option == "extra_header") {
          extra_headers.Add(value);
        } else {
          throw new ArgumentException("Unknown option " + option);
        }
      }
    }
    if (project == null) {
      throw new NullReferenceException();
    }

    // Parse all the files in our project.
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);
    FileInfo[] all_body_files = body_hpp_files.Union(cpp_files).ToArray();
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
    var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file =
          Parser.ParseFile(input_file, IsBody(input_file, extra_headers));
      file_info_to_file.Add(input_file, parser_file);
    }

    // Fix the using declarations in our project.
    foreach (var kv in file_info_to_file) {
      FileInfo input_file = kv.Key;
      Parser.File parser_file = kv.Value;
      FixPragmaOnce(parser_file);
      bool has_system_header = FixSystemHeaders(parser_file);
      RewriteFile(input_file, parser_file, dry_run);
    }
  }

  // Removes `#pragma once` directives from the file and insert an
  // `export module` directive in its place.
  private static void FixPragmaOnce(Parser.File file) {
    foreach (Node n in file.children) {
      if (n is PreprocessorDirective{ is_pragma_once: true } po) {
        var export_module =
            new ExportModule(parent: null, file.file_module_full_name);
        file.InsertChild(po.position_in_parent, export_module);
        file.RemoveChild(po);
      }
    }
  }

  // Removes include directives for system headers from the file.  Returns true
  // iff there was such an include.
  private static bool FixSystemHeaders(Parser.File file) {
    bool has_system_header = false;
    foreach (Node n in file.children) {
      if (n is Include{ is_system: true } inc) {
        file.RemoveChild(inc);
        has_system_header = true;
      }
    }
    return has_system_header;
  }
}

} // namespace sourcerer
} // namespace principia
