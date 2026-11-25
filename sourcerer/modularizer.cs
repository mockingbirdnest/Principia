using System;
using System.Collections.Generic;
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
      if (file_name_to_file_info.TryGetValue(corresponding_header_file_name,
                                             out FileInfo?
                                                     corresponding_header_file)) {
        if (corresponding_header_file != input_file) {
          corresponding_header.Add(input_file, corresponding_header_file);
        }
      }
    }

    // Fix our body files.
    var needs_import_std = new HashSet<FileInfo>();
    var body_principia_headers = new Dictionary<FileInfo, List<Include>>();
    var body_third_party_headers = new Dictionary<FileInfo, List<Include>>();
    foreach (FileInfo body_file in all_body_files) {
      if (excluded.Contains(body_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[body_file];
      FixPragmaOnce(parser_file, needs_export_module: false);

      // The `hpp_file` here might be this file if it's a test or a benchmark.
      FileInfo hpp_file =
          corresponding_header.GetValueOrDefault(body_file, body_file);

      // Remove our headers and collect them to fix the header.
      var includes = FindIncludes(parser_file);
      SplitIncludes(includes,
                    out List<Include> principia_headers,
                    out List<Include> system_headers,
                    out List<Include> third_party_headers);

      if (hpp_file == body_file) {
        Node? last_thread_party_import =
            FixImportThirdPartyHeaders(parser_file,
                                       after: null,
                                       third_party_headers);
        Node? last_std_import =
            FixImportStd(parser_file, last_thread_party_import);
        FixImportPrincipiaHeaders(parser_file,
                                  after: last_std_import,
                                  principia_headers);
      } else {
        body_principia_headers.Add(hpp_file, principia_headers);
        FixHeaders(parser_file, principia_headers);
        if (system_headers.Count > 0) {
          FixHeaders(parser_file, system_headers);
          needs_import_std.Add(hpp_file);
        }
        body_third_party_headers.Add(hpp_file, third_party_headers);
        FixHeaders(parser_file, third_party_headers);
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
      var includes = FindIncludes(parser_file);
      SplitIncludes(includes,
                    out List<Include> principia_headers,
                    out List<Include> system_headers,
                    out List<Include> third_party_headers);
      FixHeaders(parser_file, principia_headers);
      FixPragmaOnce(parser_file, needs_export_module: false);
      ExportModule? export_module =
          FixPragmaOnce(parser_file, needs_export_module: true);
      if (system_headers.Count > 0) {
        FixHeaders(parser_file, system_headers);
        needs_import_std.Add(hpp_file);
      }
      if (needs_import_std.Contains(hpp_file)) {
        FixImportStd(parser_file, export_module!);
      }
      FixHeaders(parser_file, third_party_headers);
      RewriteFile(hpp_file, parser_file, dry_run);
    }
  }

  // Removes include directives for the given headers from the file.
  private static void FixHeaders(Parser.File file, List<Include> headers) {
    foreach (Include inc in headers) {
      file.RemoveChild(inc);
    }
  }

  private static Node? FixImportPrincipiaHeaders(Parser.File file,
                                                 Node? after,
                                                 List<Include>
                                                     principia_headers) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Import? last = null;
    foreach (Include inc in principia_headers) {
      var import =
          new Import(parent: null,
                     name: "principia." + string.Join('.', inc.path));
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return last;
  }

  private static Node? FixImportStd(Parser.File file,
                                    Node? after) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = null;
    var import_std = new Import(parent: null, "std");
    file.InsertChild(position_in_parent, import_std);
    last = import_std;

    // Surround the import with blank lines if needed.
    if (!file.children[position_in_parent + 1].is_blank_line) {
      var blank = new Text("", parent: null);
      file.InsertChild(position_in_parent + 1, blank);
      last = blank;
    }
    if (position_in_parent > 0 &&
        !file.children[position_in_parent - 1].is_blank_line) {
      var blank = new Text("", parent: null);
      file.InsertChild(position_in_parent, blank);
      last = blank;
    }
    return last;
  }

  private static Node? FixImportThirdPartyHeaders(Parser.File file,
                                                  Node? after,
                                                  List<Include>
                                                      third_party_headers) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Import? last = null;
    foreach (Include inc in third_party_headers) {
      var import =
          new Import(parent: null,
                     name: '"' + string.Join('/', inc.path) + '"');
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return last;
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

  private static void SplitIncludes(List<Include> includes,
                                    out List<Include> principia_headers,
                                    out List<Include> system_headers,
                                    out List<Include> third_party_headers) {
    principia_headers = includes.FindAll(inc => inc.is_principia);
    system_headers = includes.FindAll(inc => inc.is_system);
    third_party_headers =
        includes.FindAll(inc => !inc.is_principia && !inc.is_system);
  }
}

} // namespace sourcerer
} // namespace principia
