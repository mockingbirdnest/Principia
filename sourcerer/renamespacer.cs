using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;

namespace principia {
namespace sourcerer {

class Renamespacer {

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

    // First collect all the declarations.  We'll them to rewrite the clients.
    // This must happen before we restructure any namespace, as the namespace
    // names will change and would confuse our super-fancy name resolution.
    var declaration_to_file = new Dictionary<Parser.Declaration, Parser.File>();
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name) ||
          IsBody(input_file, extra_headers)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      var exported_declarations =
          Analyser.CollectExportedDeclarations(parser_file);
      foreach (var exported_declaration in exported_declarations) {
        declaration_to_file.Add(exported_declaration, parser_file);
      }
    }

    // Process the files in client projects to fix the using declarations.
    foreach (DirectoryInfo client in clients) {
      FileInfo[] client_hpp_files = client.GetFiles("*.hpp");
      FileInfo[] client_cpp_files = client.GetFiles("*.cpp");
      FileInfo[] all_client_files =
          client_hpp_files.Union(client_cpp_files).ToArray();
      foreach (FileInfo input_file in all_client_files) {
        if (excluded.Contains(input_file.Name)) {
          continue;
        }
        bool is_body = IsBody(input_file, extra_headers);
        // This file is not in our project, so we didn't parse it yet.
        Parser.File parser_file = Parser.ParseFile(input_file, is_body);
        FixUsingDeclarations(parser_file,
                                      internal_only: !is_body,
                                      declaration_to_file);
        RewriteFile(input_file, parser_file, dry_run);
      }
    }

    // Fix the using declarations in our project.
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      FixUsingDeclarations(parser_file,
                                    internal_only: true,
                                    declaration_to_file);
    }

    // Rewrite the namespaces in our project's header files.
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name) ||
          IsBody(input_file, extra_headers)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      FixLegacyInternalNamespaces(parser_file);
      FixMissingInternalNamespaces(parser_file,
                                            insert_using_declarations: true);
      FixCompatibilityNamespace(parser_file);
      RewriteFile(input_file, parser_file, dry_run);
    }

    // Fix the body files in our project.
    foreach (FileInfo input_file in all_body_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      if (IsTest(input_file)) {
        FixUselessInternalNamespace(parser_file);
        FixFileUsingDirective(parser_file);
      } else {
        FixLegacyInternalNamespaces(parser_file);
        FixMissingInternalNamespaces(parser_file,
                                              insert_using_declarations: false);
      }
      RewriteFile(input_file, parser_file, dry_run);
    }
  }

  // Insert a using directive for |file_namespace| at the correct place in
  // |using_directives|, but only if it's not there already.  If
  // |using_directives| is empty, the insertion takes place after |after|.
  // Return true iff a directive was added.
  private static bool InsertUsingDirectiveIfNeeded(
      string file_namespace,
      Node after,
      List<UsingDirective> using_directives) {
    // Check if the using directive for the file is already present, and if
    // not determine where it should be inserted.  This assumes that the using
    // directives are sorted.
    bool file_namespace_already_exists = false;
    Node file_namespace_insertion_point = after;
    int file_namespace_insertion_index = 0;
    foreach (UsingDirective ud in using_directives) {
      if (ud.ns == file_namespace) {
        file_namespace_already_exists = true;
        break;
      } else if (string.CompareOrdinal(file_namespace, ud.ns) < 0) {
        break;
      } else {
        file_namespace_insertion_point = ud;
      }
      ++file_namespace_insertion_index;
    }
    if (file_namespace_already_exists) {
      return false;
    }

    // Insert the using directive.  Note that we must update |using_directives|.
    var parent = file_namespace_insertion_point.parent;
    Debug.Assert(parent is Namespace, "Insertion point not within a namespace");
    int insertion_point_position_in_parent =
        file_namespace_insertion_point.position_in_parent;
    var preceding_nodes_in_parent = parent.children.
        Take(insertion_point_position_in_parent + 1).ToList();
    var following_nodes_in_parent = parent.children.
        Skip(insertion_point_position_in_parent + 1).ToList();
    parent.children = preceding_nodes_in_parent;
    var using_directive = new UsingDirective(parent, file_namespace);
    parent.AddChildren(following_nodes_in_parent);
    using_directives.Insert(file_namespace_insertion_index, using_directive);
    return true;
  }

  private static void FixCompatibilityNamespace(Parser.File file) {
    var last_namespace = FindLastOutermostNamespace(file);
    if (last_namespace == null) {
      // Strange file with no namespace at all, inserting into ::.
      return;
    }
    if (last_namespace.is_compatibility_namespace) {
      // The namespace is already there.
      return;
    }
    var parent = last_namespace.parent;
    Debug.Assert(parent is Parser.File, "Last namespace not within a file");
    int last_namespace_position_in_parent = last_namespace.position_in_parent;
    var preceding_nodes_in_parent = parent.children.
        Take(last_namespace_position_in_parent + 1).ToList();
    var following_nodes_in_parent = parent.children.
        Skip(last_namespace_position_in_parent + 1).ToList();
    parent.children = preceding_nodes_in_parent;
    var blank_line_before = new Text("", parent);
    var project_namespace = new Namespace(parent,
                                          file.project_namespace_full_name);
    // No blank line after because either we have a blank line and an include,
    // or the end of the file.
    var using_directive = new UsingDirective(project_namespace,
                                             file.file_namespace_full_name);
    parent.AddChildren(following_nodes_in_parent);
  }

  private static void FixFileUsingDirective(Parser.File file) {
    var internal_using_directives =
        FindUsingDirectives(file, internal_only: false);
    var internal_using_declarations =
        FindUsingDeclarations(file, internal_only: false);
    var innermost_namespaces = FindInnermostNamespaces(file);
    var innermost_namespace = innermost_namespaces[0];
    Node after;
    if (internal_using_declarations.Count == 0) {
      // If there is no using declaration to hook our new using directive,
      // insert it after the first child of the namespace.  This assumes/ that
      // the namespace starts with a blank line.
      after = innermost_namespace.children[0];
    } else {
      after = internal_using_declarations[^1];
    }
    bool inserted_using_directive = InsertUsingDirectiveIfNeeded(
        file_namespace: file.file_namespace_full_name,
        after: after,
        using_directives: internal_using_directives);
    if (inserted_using_directive && internal_using_declarations.Count == 0) {
      // insert a blank line after the newly-inserted, self-standing using
      // directive.
      var preceding_nodes_in_innermost_namespace =
          innermost_namespace.children.Take(2).ToList();
      var following_nodes_in_innermost_namespace =
          innermost_namespace.children.Skip(2).ToList();
      innermost_namespace.children = preceding_nodes_in_innermost_namespace;
      var blank_line_before = new Text("", innermost_namespace);
      innermost_namespace.AddChildren(following_nodes_in_innermost_namespace);
    }
  }

  // Replaces the legacy "internal_foo" namespaces with a namespace "foo"
  // containing a namespace "internal".
  private static void FixLegacyInternalNamespaces(Parser.File file) {
    var legacy_internal_namespaces = FindLegacyInternalNamespaces(file);
    foreach (Namespace internal_namespace in legacy_internal_namespaces) {
      var parent = internal_namespace.parent;
      Debug.Assert(parent is Namespace,
                   "Internal namespace not within a namespace");
      int internal_position_in_parent = internal_namespace.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(internal_position_in_parent).ToList();
      var following_nodes_in_parent = parent.children.
          Skip(internal_position_in_parent + 1).ToList();
      parent.children = preceding_nodes_in_parent;
      var file_namespace = new Namespace(parent,
                                         file.file_namespace_simple_name);
      file_namespace.AddChild(internal_namespace);
      file_namespace.AddChildren(following_nodes_in_parent);
      internal_namespace.name = "internal";
    }
  }

  // Replaces the using declaration appearing in internal namespaces with
  // using directives of the appropriate file namespaces.  Ensures that the
  // using directives are deduplicated and sorted.
  private static void FixUsingDeclarations(
      Parser.File file,
      bool internal_only,
      Dictionary<Declaration, Parser.File> declaration_to_file) {
    var internal_using_declarations =
        FindUsingDeclarations(file, internal_only);
    var internal_using_directives = FindUsingDirectives(file, internal_only);
    foreach (UsingDeclaration using_declaration in
             internal_using_declarations) {
      var referenced_file =
          FindFileReferencedByUsingDeclaration(using_declaration,
                                               declaration_to_file);
      if (referenced_file == null) {
        // Not a reference to an entity that is in the project being processed.
        continue;
      }

      // Insert the using directive if not already present.
      InsertUsingDirectiveIfNeeded(
          file_namespace: referenced_file.file_namespace_full_name,
          after: internal_using_declarations[^1],
          using_directives: internal_using_directives);

      // Erase the using declaration.
      var parent = using_declaration.parent;
      Debug.Assert(parent is Namespace,
                   "Using declaration not within a namespace");
      int using_position_in_parent = using_declaration.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(using_position_in_parent).ToList();
      var following_nodes_in_parent = parent.children.
          Skip(using_position_in_parent + 1).ToList();
      parent.children = preceding_nodes_in_parent;
      parent.AddChildren(following_nodes_in_parent);
    }
  }

  // Finds the innermost namespaces that are not internal, and moves all of
  // their declaration into a nested "internal" namespace.  Adds using
  // declarations to make all the now-internal declarations accessible to the
  // outside world.
  private static void FixMissingInternalNamespaces(
      Parser.File file,
      bool insert_using_declarations) {
    var innermost_namespaces = FindInnermostNamespaces(file);
    foreach (var ns in innermost_namespaces) {
      if (!ns.is_internal && !ns.is_compatibility_namespace) {
        var nodes_in_ns = ns.children.ToList();
        // Check if there are names to export from this namespace.  If not,
        // don't touch it, unless it's in a body.  It may be some local
        // reopening.  This produces a deduped and sorted set of symbols.
        var names = new SortedSet<string>();
        foreach (Node n in nodes_in_ns) {
          if (n is Declaration decl) {
            names.Add(decl.name);
          }
        }
        if (names.Count == 0 && insert_using_declarations) {
          continue;
        }
        ns.children.Clear();
        var file_namespace = new Namespace(ns,
                                           file.file_namespace_simple_name);
        var internal_namespace = new Namespace(file_namespace, "internal");
        internal_namespace.children.AddRange(nodes_in_ns);

        // Insert the using declarations.
        if (insert_using_declarations) {
          var blank_line_before = new Text("", file_namespace);
          foreach (string name in names) {
            var using_declaration =
                new UsingDeclaration(file_namespace, "internal::" + name);
          }
          var blank_line_after = new Text("", file_namespace);
        }
      }
    }
  }

  private static void FixUselessInternalNamespace(Parser.File file) {
    var legacy_internal_namespaces = FindLegacyInternalNamespaces(file);
    foreach (Namespace internal_namespace in legacy_internal_namespaces) {
      var parent = internal_namespace.parent;
      Debug.Assert(parent is Namespace,
                   "Internal namespace not within a namespace");
      int internal_position_in_parent = internal_namespace.position_in_parent;
      var preceding_nodes_in_parent =
          parent.children.Take(internal_position_in_parent).ToList();
      var nodes_in_internal_namespace = internal_namespace.children.ToList();
      var following_nodes_in_parent = parent.children.
          Skip(internal_position_in_parent + 1).ToList();
      parent.children = preceding_nodes_in_parent;
      parent.AddChildren(nodes_in_internal_namespace);
      parent.AddChildren(following_nodes_in_parent);
    }
  }

  private static void RewriteFile(FileInfo input_file,
                                  Parser.File file,
                                  bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";

    using (StreamWriter writer = System.IO.File.CreateText(output_filename)) {
      RewriteNode(writer, file);
    }
    if (!dry_run) {
      System.IO.File.Move(output_filename, input_filename, overwrite: true);
    }
  }

  private static void RewriteNode(StreamWriter writer, Parser.Node node) {
    foreach (Parser.Node child in node.children) {
      if (child.must_rewrite) {
        writer.WriteLine(child.Cxx());
      } else {
        writer.WriteLine(child.text);
      }
      if (child is Parser.Namespace ns) {
        RewriteNode(writer, child);
        if (ns.must_rewrite) {
          writer.WriteLine(ns.ClosingCxx());
        } else {
          writer.WriteLine(ns.closing_text);
        }
      }
    }
  }
}

} // namespace sourcerer
} // namespace principia
