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
    //TODO(phl):Sorted set!
    var body_principia_includes = new Dictionary<FileInfo, List<Include>>();
    var body_third_party_includes = new Dictionary<FileInfo, List<Include>>();
    foreach (FileInfo body_file in all_body_files) {
      if (excluded.Contains(body_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[body_file];

      // The `hpp_file` here might be this file if it's a test or a benchmark.
      FileInfo hpp_file =
          corresponding_header.GetValueOrDefault(body_file, body_file);

      // Only process includes up to the first outermost namespace.  We don't
      // want to touch the final include (for the body) or includes of generated
      // code.
      var outermost_namespaces = FindOutermostNamespaces(parser_file);
      Namespace? first_namespace =
          outermost_namespaces.Count > 0 ? outermost_namespaces[0] : null;

      Node? last_export_module =
          FixPragmaOnce(parser_file,
                        hpp_file != body_file
                            ? file_info_to_file[hpp_file]
                            : null,
                        needs_export_module: false);

      // Find all the includes and categorize them.
      SplitIncludes(FindIncludes(parser_file, stop: first_namespace),
                    out List<Include> principia_includes,
                    out List<Include> system_includes,
                    out List<Include> third_party_includes);

      if (hpp_file == body_file) {
        Node? last_std_import;
        if (system_includes.Count > 0) {
          last_std_import = FixImportStd(parser_file, last_export_module);
        } else {
          last_std_import = last_export_module;
        }
        // If this is a test or benchmark, its imports go right here.
        Node? last_third_party_import =
            FixImportThirdPartyHeaders(parser_file,
                                       after: last_std_import,
                                       third_party_includes);
        FixImportPrincipiaHeaders(parser_file,
                                  after: last_third_party_import,
                                  principia_includes);
      } else {
        // If this is a body that has a header, its imports will go in the
        // header.
        body_principia_includes.Add(hpp_file, principia_includes);
        if (system_includes.Count > 0) {
          needs_import_std.Add(hpp_file);
        }
        body_third_party_includes.Add(hpp_file, third_party_includes);
      }

      // In all cases, the includes must go away.
      FixIncludes(parser_file, principia_includes);
      FixIncludes(parser_file, system_includes);
      FixIncludes(parser_file, third_party_includes);

      RewriteFile(body_file, parser_file, dry_run);
    }

    // Fix our non-body header files.  This must come last because we need the
    // information that we collected from the body files.
    foreach (FileInfo hpp_file in hpp_files) {
      if (excluded.Contains(hpp_file.Name) || IsBodyHpp(hpp_file)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[hpp_file];

      // Only process includes up to the first outermost namespace.  We don't
      // want to touch the final include (for the body) or includes of generated
      // code.
      var outermost_namespaces = FindOutermostNamespaces(parser_file);
      Namespace? first_namespace =
          outermost_namespaces.Count > 0 ? outermost_namespaces[0] : null;

      Node? last_export_module =
          FixPragmaOnce(parser_file, hpp_file: null, needs_export_module: true);

      SplitIncludes(FindIncludes(parser_file, stop: first_namespace),
                    out List<Include> principia_includes,
                    out List<Include> system_includes,
                    out List<Include> third_party_includes);
      var all_principia_includes = new List<Include>(principia_includes);
      all_principia_includes.AddRange(
          body_principia_includes.GetValueOrDefault(
              hpp_file,
              new List<Include>()));
      var all_third_party_includes = new List<Include>(third_party_includes);
      third_party_includes.AddRange(
          body_third_party_includes.GetValueOrDefault(
              hpp_file,
              new List<Include>()));

      // Add the imports collected from this file and its body.
      Node? last_std_import;
      if (system_includes.Count > 0 || needs_import_std.Contains(hpp_file)) {
        last_std_import = FixImportStd(parser_file, last_export_module);
      } else {
        last_std_import = last_export_module;
      }
      Node? last_third_party_import =
          FixImportThirdPartyHeaders(parser_file,
                                     after: last_std_import,
                                     all_third_party_includes);
      FixImportPrincipiaHeaders(parser_file,
                                after: last_third_party_import,
                                all_principia_includes);

      // Drop the includes.
      FixIncludes(parser_file, principia_includes);
      FixIncludes(parser_file, system_includes);
      FixIncludes(parser_file, third_party_includes);

      // Export the non-internal using declarations.
      FixExportedUsingDeclarations(parser_file);

      // Export everything from the internal namespaces except the using
      // directives.
      FixInternalNamespacesExportSequence(parser_file);

      RewriteFile(hpp_file, parser_file, dry_run);
    }
  }

  // Ensures that `after` is followed by a blank line, unless it is `null` or
  // the last child of `file`.  Returns the blank line, or `null` or the last
  // child.
  private static Node? EnsureBlankLineAfter(Node? after) {
    if (after != null) {
      Node parent = after.parent!;
      if (after.is_blank_line) {
        // If `after` is already a blank line we have nothing to do.
        return after;
      }
      if (after.position_in_parent + 1 < parent.children.Count) {
        after = parent.children[after.position_in_parent + 1];
        if (!after.is_blank_line) {
          // If `after` is not followed by a blank line, insert one and move
          // `after` to that blank line.
          var blank = Text.Blank();
          parent.InsertChild(after.position_in_parent, blank);
          after = blank;
        }
      }
    }
    return after;
  }

  private static void FixExportedUsingDeclarations(Parser.File file) {
    var using_declarations = FindUsingDeclarations(file, internal_only: false);
    foreach (UsingDeclaration using_declaration in using_declarations) {
      Node parent = using_declaration.parent!;
      int position_in_parent = using_declaration.position_in_parent;
      var export = new Export(parent: null);
      parent.RemoveChild(using_declaration);
      export.AddChild(using_declaration);
      parent.InsertChild(position_in_parent, export);
    }
  }

  // Removes include directives for the given headers from the file.
  private static void FixIncludes(Parser.File file, List<Include> headers) {
    foreach (Include inc in headers) {
      int position_in_parent = inc.position_in_parent;
      file.RemoveChild(inc);
      // Remove blank lines that follow the include.
      while (position_in_parent < file.children.Count &&
             file.children[position_in_parent].is_blank_line) {
        file.RemoveChild(file.children[position_in_parent]);
      }
    }
  }

  private static void FixInternalNamespacesExportSequence(Parser.File file) {
    var internal_using_directives =
        FindUsingDirectives(file, internal_only: true);
    Node? after = null;
    if (internal_using_directives.Count == 0) {
      var innermost_namespaces = FindInnermostNamespaces(file);
      if (innermost_namespaces.Count == 0) {
        // No using directives and no namespaces, this is a weird file, leave it
        // alone.
        return;
      } else {
        // Insert the export sequence after the first innermost namespace.
        after = innermost_namespaces[0];
      }
    } else {
      // Insert the export sequence after the last internal using directive.
      after = internal_using_directives[^1];
    }
    Node parent = after.parent!;
    after = EnsureBlankLineAfter(after);
    var es = new ExportSequence(parent: null);
    parent.InsertChild(after!.position_in_parent + 1, es);
    EnsureBlankLineAfter(es);
  }

  private static Node? FixImportPrincipiaHeaders(Parser.File file,
                                                 Node? after,
                                                 List<Include>
                                                     principia_includes) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = after;
    foreach (Include inc in principia_includes) {
      var import =
          new Import(parent: null,
                     name: "principia." + string.Join('.', inc.path));
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return EnsureBlankLineAfter(last);
  }

  private static Node? FixImportStd(Parser.File file,
                                    Node? after) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = null;
    var import_std = new Import(parent: null, "std");
    file.InsertChild(position_in_parent, import_std);
    last = import_std;

    // Surround the import with blank lines if needed.
    bool needs_blank_line_before =
        position_in_parent > 0 &&
        !file.children[position_in_parent - 1].is_blank_line;
    bool needs_blank_line_after =
        position_in_parent + 1 < file.children.Count &&
        !file.children[position_in_parent + 1].is_blank_line;

    if (needs_blank_line_before) {
      var blank = Text.Blank();
      file.InsertChild(import_std.position_in_parent, blank);
      last = blank;
    }
    if (needs_blank_line_after) {
      var blank = Text.Blank();
      file.InsertChild(import_std.position_in_parent + 1, blank);
      last = blank;
    }
    return last;
  }

  private static Node? FixImportThirdPartyHeaders(Parser.File file,
                                                  Node? after,
                                                  List<Include>
                                                      third_party_includes) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = after;
    foreach (Include inc in third_party_includes) {
      var import =
          new Import(parent: null,
                     name: '"' +
                           string.Join('/', inc.path) +
                           "."  +
                           inc.extension +
                           '"');
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return EnsureBlankLineAfter(last);
  }

  // Removes `#pragma once` directives from the file and insert an
  // `export module` directive in its place if needed.
  private static Node? FixPragmaOnce(Parser.File file,
                                     Parser.File? hpp_file,
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
    Node? node = null;
    foreach (PreprocessorDirective po in pragma_onces) {
      if (needs_export_module) {
        int position_in_parent = po.position_in_parent;
        node = new Export(parent: null);
        var module =
            new Module(parent: node, name: file.file_module_full_name);
        file.InsertChild(position_in_parent, node);
      } else if (hpp_file != null) {
        int position_in_parent = po.position_in_parent;
        node =
            new Text(text: "// Body of link:" + hpp_file.file_link_name,
                     parent: null);
        file.InsertChild(position_in_parent, node);
      } else {
        // If the `#pragma once` is not replaced by an `export module` or a
        // link, remove the blank lines after it.
        foreach (Text blank in blank_lines_after_pragma_once) {
          file.RemoveChild(blank);
        }
      }
      file.RemoveChild(po);
    }
    return EnsureBlankLineAfter(node);
  }

  private static void SplitIncludes(List<Include> includes,
                                    out List<Include> principia_includes,
                                    out List<Include> system_includes,
                                    out List<Include> third_party_includes) {
    principia_includes = includes.FindAll(inc => inc.is_principia);
    system_includes = includes.FindAll(inc => inc.is_system);
    third_party_includes =
        includes.FindAll(inc => !inc.is_principia && !inc.is_system);
  }
}

} // namespace sourcerer
} // namespace principia
