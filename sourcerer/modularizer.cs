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
            arg.Split([ "--", ":" ], StringSplitOptions.None);
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
    var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    var file_name_to_file_info = new Dictionary<string, FileInfo>();
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);
    FileInfo[] all_body_files = [.. body_hpp_files.Union(cpp_files)];
    FileInfo[] all_files = [.. hpp_files.Union(cpp_files)];
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file =
          Parser.ParseFile(input_file, IsBody(input_file, extra_headers));
      file_info_to_file.Add(input_file, parser_file);
      file_name_to_file_info.Add(Path.GetFullPath(input_file.Name,
                                                  input_file.DirectoryName!),
                                 input_file);
    }

    // Map the bodies to their headers.
    var corresponding_header = new Dictionary<FileInfo, FileInfo>();
    foreach (FileInfo input_file in all_body_files) {
      string corresponding_header_file_name =
          Filenames.CorrespondingHeader(input_file);
      // Benchmarks and tests don't have a corresponding header.
      if (file_name_to_file_info.TryGetValue(corresponding_header_file_name, out FileInfo? corresponding_header_file)) {
        if (corresponding_header_file != input_file) {
          corresponding_header.Add(input_file, corresponding_header_file);
        }
      }
    }

    // Fix our body files.
    var needs_import_std = new HashSet<FileInfo>();
    foreach (FileInfo body_file in all_body_files) {
      if (excluded.Contains(body_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[body_file];
      FixPragmaOnce(parser_file, needs_export_module: false);
      bool has_system_header = FixSystemHeaders(parser_file);
      if (has_system_header) {
        needs_import_std.Add(
            corresponding_header.GetValueOrDefault(body_file, body_file));
      }
      RewriteFile(body_file, parser_file, dry_run);
    }

    // Fix our non-body header files.  This must come last because we need the
    // information that we collected from the body files.
    foreach (FileInfo hpp_file in hpp_files) {
      if (excluded.Contains(hpp_file.Name) || IsBodyHpp(hpp_file)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[hpp_file];
      ExportModule? export_module =
          FixPragmaOnce(parser_file, needs_export_module: true);
      bool has_system_header = FixSystemHeaders(parser_file);
      if (has_system_header || needs_import_std.Contains(hpp_file)) {
        FixImportStd(parser_file, export_module!);
      }
      RewriteFile(hpp_file, parser_file, dry_run);
    }
  }

  private static void FixImportStd(Parser.File file,
                                   ExportModule export_module) {
    var import_std = new Import(parent: null, "std");
    int position_in_parent = export_module.position_in_parent + 1;
    file.InsertChild(position_in_parent, import_std);

    // Surround the import with blank lines if needed.
    var blank = new Text("", parent: null);
    if (!file.children[position_in_parent + 1].is_blank_line) {
      file.InsertChild(position_in_parent + 1, blank);
    }
    if (!file.children[position_in_parent - 1].is_blank_line) {
      file.InsertChild(position_in_parent, blank);
    }
  }

  // Removes `#pragma once` directives from the file and insert an
  // `export module` directive in its place if needed.
  private static ExportModule? FixPragmaOnce(Parser.File file,
                                             bool needs_export_module) {
    // Collect the pragmas to fix.
    var pragma_onces = new List<PreprocessorDirective>();
    var blank_lines_after_pragma_once = new List<Text>();
    bool right_after_pragma_once = false;
    foreach (Node n in file.children) {
      if (n is PreprocessorDirective{ is_pragma_once: true } po) {
        pragma_onces.Add(po);
        right_after_pragma_once = true;
      } else if (right_after_pragma_once && n.is_blank_line) {
        blank_lines_after_pragma_once.Add((n as Text)!);
      } else {
        right_after_pragma_once = false;
      }
    }

    // And fix them.
    ExportModule? export_module = null;
    foreach (PreprocessorDirective po in pragma_onces) {
      if (needs_export_module) {
        int position_in_parent = po.position_in_parent;
        export_module =
            new ExportModule(parent: null, name: file.file_module_full_name);
        file.InsertChild(position_in_parent, export_module);
      } else {
        // If the `#pragma once` is not replaced by an `export module`, remove
        // the blank lines after it.
        foreach (Text blank in blank_lines_after_pragma_once) {
          file.RemoveChild(blank);
        }
      }
      file.RemoveChild(po);
    }
    return export_module;
  }

  // Removes include directives for system headers from the file.  Returns true
  // iff there was such an include.
  private static bool FixSystemHeaders(Parser.File file) {
    bool has_system_header = false;
    var system_headers = new List<Include>();
    foreach (Node n in file.children) {
      if (n is Include{ is_system: true } inc) {
        system_headers.Add(inc);
        has_system_header = true;
      }
    }
    foreach (Include inc in system_headers) {
      file.RemoveChild(inc);
    }
    return has_system_header;
  }
}

} // namespace sourcerer
} // namespace principia
