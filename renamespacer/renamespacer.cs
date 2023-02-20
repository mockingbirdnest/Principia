using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {

class Parser {
  public abstract class Node {
    // This links the new node as the last child of its parent.
    protected Node(int line_number, Node parent) {
      this.line_number = line_number;
      this.parent = parent;
      this.children = new List<Node>();
      if (parent != null) {
        position_in_parent = parent.children.Count;
        parent.children.Add(this);
      }
    }

    // This function must return a string ending with a new line.
    public virtual string Cxx(bool is_at_exit) {
      return "--FATAL--" + Environment.NewLine;
    }

    // Writes a single node.
    public abstract void WriteNode(string indent = "");

    public void WriteTree(string indent = "") {
      WriteNode(indent);
      foreach (Node child in children) {
        child.WriteTree(indent + "  ");
      }
    }

    public virtual bool must_rewrite {
      get => must_rewrite_;
      set { must_rewrite_ = value; }
    }

    public int line_number = 0;
    public int? last_line_number = null;
    public Node parent = null;
    public int position_in_parent = -1;
    public List<Node> children = null;
    protected bool must_rewrite_ = false;
  }

  public abstract class Declaration : Node {
    protected Declaration(int line_number, Node parent, string name) : base(
        line_number,
        parent) {
      this.name = name;
    }

    public string name = null;
  }

  public class Class : Declaration {
    public Class(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Class " + name);
    }
  }

  public class Constant : Declaration {
    public Constant(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) { }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Constant " + name);
    }
  }

  // A placeholder for a deleted node.
  public class Empty : Node {
    public Empty(int line_number, Node parent) : base(line_number, parent) {}

    public override string Cxx(bool is_at_exit) {
      return "";
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Empty");
    }

    public override bool must_rewrite => true;
  }

  public class File : Node {
    public File(int line_number, FileInfo file_info) : base(
        line_number,
        parent: null) {
      this.file_info = file_info;
      file_namespace = Regex.Replace(
          file_info.Name,
          "(_body|_test)?\\.[hc]pp",
          "");
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "File " + file_info.FullName);
    }

    public FileInfo file_info = null;
    public string file_namespace = null;
  }

  public class Function : Declaration {
    public Function(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) { }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Function " + name);
    }
  }

  public class Include : Node {
    public Include(int line_number, Node parent, string[] path) : base(
        line_number,
        parent) {
      this.path = path;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Include " + string.Join(", ", path));
    }

    public string[] path = null;
  }

  public class Namespace : Declaration {
    public Namespace(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) {
      if (parent is Namespace{ is_internal: true } ) {
        is_internal = true;
      } else {
        is_internal = name.StartsWith("internal");
      }
    }

    public override string Cxx(bool is_at_exit) {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      if (is_at_exit) {
        return "}  // namespace " + name + Environment.NewLine;
      } else {
        return "namespace " + name + " {" + Environment.NewLine;
      }
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent +
                        "Namespace " +
                        name +
                        (is_internal ? "Internal" : ""));
    }

    public bool is_internal = false;
  }

  public class UsingDeclaration : Declaration {
    public UsingDeclaration(int line_number, Node parent, string full_name) :
        base(line_number, parent, full_name.Split("::")[^1]) {
      this.full_name = full_name;
      string[] segments = full_name.Split("::");

      // Try to figure out if the name was declared in a preceding namespace.
      // This is useful to fix up the internal namespaces.
      if (segments.Length == 2) {
        if (parent is Namespace ns) {
          foreach (Node sibling in ns.children) {
            if (sibling is Namespace nested_ns &&
                nested_ns.name == segments[0]) {
              declared_in_namespace = nested_ns;
              break;
            }
          }
        }
      }
    }

    public override string Cxx(bool is_at_exit) {
      Debug.Assert(!is_at_exit, "no exit for using");
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      if (declared_in_namespace == null) {
        return "using " + full_name + ";" + Environment.NewLine;
      } else {
        return "using " +
               declared_in_namespace.name +
               "::" +
               name +
               ";" +
               Environment.NewLine;
      }
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent +
                        "UsingDeclaration " +
                        name +
                        (declared_in_namespace == null
                             ? ""
                             : " from " + declared_in_namespace.name));
    }

    public override bool must_rewrite =>
        must_rewrite_ || declared_in_namespace is { must_rewrite: true };

    public string full_name = null;
    public Namespace declared_in_namespace = null;
  }

  public class UsingDirective : Node {
    public UsingDirective(int line_number, Node parent, string ns) : base(
        line_number,
        parent) {
      this.ns = ns;
    }

    public override string Cxx(bool is_at_exit) {
      Debug.Assert(!is_at_exit, "no exit for using");
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      return "using namespace " + ns + ";" + Environment.NewLine;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "UsingDirective " + ns);
    }

    public string ns = null;
  }

  public class Struct : Declaration {
    public Struct(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Struct " + name);
    }
  }

  public class TypeAlias : Declaration {
    public TypeAlias(int line_number, Node parent, string name) : base(
        line_number,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "TypeAlias " + name);
    }
  }

  private static bool IsClass(string line) {
    return Regex.IsMatch(line, @"^class [\w]+[ ;].*$");
  }

  private static bool IsClosingNamespace(string line) {
    return line != "}  // namespace" && line.StartsWith("}  // namespace");
  }

  private static bool IsConstant(string line) {
    return Regex.IsMatch(line, @"^constexpr .* = .*$");
  }

  private static bool IsFunction(string line) {
    return Regex.IsMatch(line, @"^[\w].+ [^: ]+\(.*$");
  }

  private static bool IsOpeningNamespace(string line) {
    return line != "namespace {" &&
           line.StartsWith("namespace ") &&
           !Regex.IsMatch(line, @"^namespace [\w]+ = .*$");
  }

  private static bool IsOwnHeaderInclude(string line, FileInfo input_file) {
    string own_header = Regex.Replace(
        input_file.Name,
        @"(_body|_test)?\.[hc]pp",
        ".hpp");
    return line == "#include \"" + own_header + "\"";
  }

  private static bool IsPrincipiaInclude(string line) {
    // Principia header files end in .hpp.
    return line.StartsWith("#include \"") && line.EndsWith(".hpp\"");
  }

  private static bool IsStruct(string line) {
    return Regex.IsMatch(line, @"^struct [\w]+[ ;].*$");
  }

  private static bool IsTypeAlias(string line) {
    return Regex.IsMatch(line, @"^using [\w]+ =;$");
  }

  private static bool IsUsingDeclaration(string line) {
    return !IsUsingDirective(line) && Regex.IsMatch(line, @"^using [\w:]+;$");
  }

  private static bool IsUsingDirective(string line) {
    return line.StartsWith("using namespace ");
  }

  private static string ParseClass(string line) {
    return Regex.Replace(line.Replace("class ", ""), @"[; ].*$", "");
  }

  private static string ParseClosingNamespace(string line) {
    return line.Replace("}  // namespace ", "");
  }

  private static string ParseConstant(string line) {
    return Regex.Replace(Regex.Replace(line, @"^constexpr [^ ]+ ", ""),
                         @" = .*$",
                         "");
  }

  private static string ParseFunction(string line) {
    return Regex.Replace(Regex.Replace(line, @"\(.*$", ""), @"^.+ ", "");
  }

      private static string[] ParseIncludedPath(string line) {
    string path = line.Replace("#include \"", "").Replace(".hpp\"", "");
    return path.Split('/');
  }

  private static string ParseOpeningNamespace(string line) {
    return line.Replace("namespace ", "").Replace(" {", "");
  }

  private static string ParseStruct(string line) {
    return Regex.Replace(line.Replace("struct ", ""), @"[; ].*$", "");
  }

  private static string ParseTypeAlias(string line) {
    return Regex.Replace(line.Replace("using ", ""), @" =.*$", "");
  }

  private static string ParseUsingDeclaration(string line) {
    return line.Replace("using ", "").Replace(";", "");
  }

  private static string ParseUsingDirective(string line) {
    return line.Replace("using namespace ", "").Replace(";", "");
  }

  public static File ParseFile(FileInfo file_info) {
    var file = new File(line_number: 0, file_info);
    Node current = file;

    using (StreamReader reader = file_info.OpenText()) {
      int line_number = 1;
      while (!reader.EndOfStream) {
        string line = reader.ReadLine();
        if (IsPrincipiaInclude(line) && !IsOwnHeaderInclude(line, file_info)) {
          var include = new Include(line_number,
                                    parent: current,
                                    ParseIncludedPath(line));
        } else if (IsOpeningNamespace(line)) {
          current = new Namespace(line_number,
                                  parent: current,
                                  ParseOpeningNamespace(line));
        } else if (IsClosingNamespace(line)) {
          var name = ParseClosingNamespace(line);
          if (current is Namespace ns) {
            Debug.Assert(ns.name == name);
            ns.last_line_number = line_number;
          } else {
            Debug.Assert(false);
          }
          current = current.parent;
        } else if (IsClass(line)) {
          var klass = new Class(line_number, parent: current, ParseClass(line));
        } else if (IsConstant(line)) {
          var constant =
              new Constant(line_number, parent: current, ParseConstant(line));
        } else if (IsFunction(line)) {
          var function =
              new Function(line_number, parent: current, ParseFunction(line));
        } else if (IsStruct(line)) {
          var strukt = new Struct(
              line_number,
              parent: current,
              ParseStruct(line));
        } else if (IsTypeAlias(line)) {
          var type_alias = new TypeAlias(
              line_number,
              parent: current,
              ParseTypeAlias(line));
        } else if (IsUsingDirective(line)) {
          var using_directive = new UsingDirective(
              line_number,
              parent: current,
              ParseUsingDirective(line));
        } else if (IsUsingDeclaration(line)) {
          var using_declaration = new UsingDeclaration(
              line_number,
              parent: current,
              ParseUsingDeclaration(line));
        }
        ++line_number;
      }
      file.last_line_number = line_number;
    }
    return file;
  }

  private static string NamespaceForFile(FileInfo file_info) {
    return "principia::" +
           file_info.Directory.Name +
           "::" +
           Regex.Replace(file_info.Name, @"\.hpp|\.cpp", "");
  }

  public static List<Declaration> CollectExportedDeclarations(Node node) {
    var exported_declarations = new List<Declaration>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (!ns.is_internal) {
          exported_declarations.AddRange(CollectExportedDeclarations(child));
        }
      } else if (child is Declaration decl) {
        exported_declarations.Add(decl);
      }
    }
    return exported_declarations;
  }

  private static FileInfo FindFileReferencedByUsingDeclaration(
      UsingDeclaration using_declaration,
      Dictionary<Declaration, FileInfo> declaration_to_file) {
    var referenced_name = using_declaration.name;
    var referenced_innermost_namespace =
        using_declaration.full_name.Split("::")[^2];
    foreach (var pair in declaration_to_file) {
      var declaration = pair.Key;
      // Poor man's name resolution.
      if (declaration.name == referenced_name &&
          declaration.parent is Namespace ns &&
          ns.name == referenced_innermost_namespace) {
        return pair.Value;
      }
    }
    return null;
  }

  private static List<Namespace> FindInnermostNamespaces(Node node) {
    var innermost_namespaces = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        var nested_namespaces = FindInnermostNamespaces(child);
        if (nested_namespaces.Count == 0) {
          innermost_namespaces.Add(ns);
        } else {
          innermost_namespaces.AddRange(nested_namespaces);
        }
      }
    }
    return innermost_namespaces;
  }

  private static List<UsingDeclaration>
      FindInternalUsingDeclarations(Node node) {
    var internal_using_declarations = new List<UsingDeclaration>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDeclaration ud) {
              internal_using_declarations.Add(ud);
            }
          }
        }
        internal_using_declarations.AddRange(
            FindInternalUsingDeclarations(child));
      }
    }
    return internal_using_declarations;
  }

  private static List<UsingDirective>
      FindInternalUsingDirectives(Node node) {
    var internal_using_declarations = new List<UsingDirective>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDirective ud) {
              internal_using_declarations.Add(ud);
            }
          }
        }
        internal_using_declarations.AddRange(
            FindInternalUsingDirectives(child));
      }
    }
    return internal_using_declarations;
  }

  private static List<Namespace> FindLegacyInternalNamespaces(
      File file,
      Node node) {
    var legacy_internal_namespaces = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace internal_namespace &&
          internal_namespace.name.StartsWith("internal_")) {
        legacy_internal_namespaces.Add(internal_namespace);
      } else if (child is Namespace ns) {
        legacy_internal_namespaces.AddRange(
            FindLegacyInternalNamespaces(file, child));
      }
    }
    return legacy_internal_namespaces;
  }

  // Replaces the legacy "internal_foo" namespaces with a namespace "foo"
  // containing a namespace "internal".
  public static void FixLegacyInternalNamespaces(File file) {
    var legacy_internal_namespaces = FindLegacyInternalNamespaces(file, file);
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
      var file_namespace = new Namespace(internal_namespace.line_number,
                                         parent,
                                         file.file_namespace);
      file_namespace.children.Add(internal_namespace);
      file_namespace.children.AddRange(following_nodes_in_parent);
      file_namespace.last_line_number = parent.last_line_number;
      file_namespace.must_rewrite = true;
      internal_namespace.parent = file_namespace;
      internal_namespace.position_in_parent = 0;
      internal_namespace.name = "internal";
      internal_namespace.must_rewrite = true;
    }
  }

  // Replaces the using declaration appearing in internal namespaces with
  // using directives of the appropriate file namespaces.  Ensures that the
  // using directives are deduplicated and sorted.
  public static void FixInternalUsingDeclarations(
      File file,
      Dictionary<Declaration, FileInfo> declaration_to_file) {
    var internal_using_declarations = FindInternalUsingDeclarations(file);
    var internal_using_directives = FindInternalUsingDirectives(file);
    foreach (UsingDeclaration using_declaration in
             internal_using_declarations) {
      var file_info =
          FindFileReferencedByUsingDeclaration(using_declaration,
                                               declaration_to_file);
      if (file_info == null) {
        // Not a reference to an entity that is in the project being processed.
        continue;
      }

      // Check if the using directive for the file is already present, and if
      // not determine where it should be inserted.  This assumes that the using
      // directives are sorted.
      string file_namespace = NamespaceForFile(file_info);
      bool file_namespace_already_exists = false;
      Node file_namespace_insertion_point = internal_using_declarations[0];
      int file_namespace_insertion_index = 0;
      foreach (UsingDirective ud in internal_using_directives) {
        if (ud.ns == file_namespace) {
          file_namespace_already_exists = true;
          break;
        } else if (string.CompareOrdinal(file_namespace, ud.ns) < 0) {
          file_namespace_insertion_point = ud;
          break;
        }
        ++file_namespace_insertion_index;
      }
      if (file_namespace_already_exists) {
        continue;
      }

      // Insert the using directive.  Note that we must update
      // |internal_using_directives|.
      {
        var parent = file_namespace_insertion_point.parent;
        Debug.Assert(parent is Namespace,
                     "Insertion point not within a namespace");
        int insertion_point_position_in_parent =
            file_namespace_insertion_point.position_in_parent;
        var preceding_nodes_in_parent = parent.children.
            Take(insertion_point_position_in_parent).ToList();
        var following_nodes_in_parent = parent.children.
            Skip(insertion_point_position_in_parent).ToList();
        parent.children = preceding_nodes_in_parent;
        var using_directive = new UsingDirective(
            file_namespace_insertion_point.line_number,
            parent,
            file_namespace);
        using_directive.must_rewrite = true;
        internal_using_directives.Insert(file_namespace_insertion_index,
                                         using_directive);
        foreach (Node n in following_nodes_in_parent) {
          ++n.position_in_parent;
        }
        parent.children.AddRange(following_nodes_in_parent);
      }

      // Erase the using declaration.
      {
        var parent = using_declaration.parent;
        Debug.Assert(parent is Namespace,
                     "Using declaration not within a namespace");
        int using_position_in_parent = using_declaration.position_in_parent;
        var preceding_nodes_in_parent =
            parent.children.Take(using_position_in_parent).ToList();
        var following_nodes_in_parent = parent.children.
            Skip(using_position_in_parent + 1).ToList();
        parent.children = preceding_nodes_in_parent;
        var empty = new Empty(using_declaration.line_number, parent);
        parent.children.AddRange(following_nodes_in_parent);
      }
    }
  }

  // Finds the innermost namespaces that are not internal, and moves all of
  // their declaration into a nested "internal" namespace.  Adds using
  // declarations to make all the now-internal declarations accessible to the
  // outside world.
  public static void FixMissingInternalNamespaces(File file) {
    var innermost_namespaces = FindInnermostNamespaces(file);
    foreach (var ns in innermost_namespaces) {
      if (!ns.is_internal) {
        var nodes_in_ns = ns.children.ToList();
        ns.children.Clear();
        ns.must_rewrite = true;
        var internal_namespace = new Namespace(ns.line_number, ns, "internal");
        internal_namespace.children.AddRange(nodes_in_ns);
        internal_namespace.last_line_number = ns.last_line_number;
        internal_namespace.must_rewrite = true;

        // Insert the using declarations.  First dedupe and sort the symbols.
        var names = new SortedSet<string>();
        foreach (Node n in nodes_in_ns) {
          if (n is Declaration decl) {
            names.Add(decl.name);
          }
        }
        foreach (string name in names) {
          var using_declaration =
              new UsingDeclaration(ns.last_line_number.Value,
                                   ns,
                                   "internal::" + name);
          using_declaration.must_rewrite = true;
        }
      }
    }
  }
}

class Renamespacer {
  static void RewriteFile(FileInfo input_file, Parser.File file, bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";

    Parser.Node parent = file;
    int child_position = 0;
    Parser.Node node = parent.children[child_position];
    int node_line_number = node.line_number;
    using (StreamReader reader = input_file.OpenText()) {
      using (StreamWriter writer = File.CreateText(output_filename)) {
        int line_number = 1;
        bool is_at_exit = false;
        while (!reader.EndOfStream) {
          string line = reader.ReadLine();
          Debug.Assert(node_line_number >= line_number);

          bool has_rewritten = false;
          while (node_line_number == line_number) {
            if (node.must_rewrite) {
              has_rewritten = true;
              writer.Write(node.Cxx(is_at_exit));
            } else {
              has_rewritten = false;
            }

            if (node is Parser.Namespace ns && ns.line_number == line_number) {
              // Entering a namespace.
              parent = node;
              child_position = -1;
            }
            ++child_position;
            if (child_position == parent.children.Count) {
              // Exiting a namespace.
              node = parent;
              node_line_number = node.last_line_number.Value;
              child_position = node.position_in_parent;
              parent = node.parent;
              is_at_exit = true;
            } else {
              // Moving to the next sibling.
              node = parent.children[child_position];
              node_line_number = node.line_number;
              is_at_exit = false;
            }
          }

          if (!has_rewritten) {
            writer.WriteLine(line);
          }
          ++line_number;
        }
      }
    }
    if (!dry_run) {
      File.Move(output_filename, input_filename, overwrite: true);
    }
  }

  // Usage:
  //   renamespacer --project:quantities \
  //                --client:base --client:physics \
  //                --exclude:macros.hpp --dry_run:false
  // This will renamespace quantities and fix the references in the client
  // projects.  The files will be overwritten.
  static void Main(string[] args) {
    // Parse the arguments.
    DirectoryInfo project = null;
    var clients = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
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
        } else if (option == "exclude") {
          excluded.Add(value);
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        }
      }
    }

    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    var hpp_parsed_files = new Dictionary<FileInfo, Parser.File>();
    var declaration_to_file = new Dictionary<Parser.Declaration, FileInfo>();
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = Parser.ParseFile(input_file);
      Parser.FixLegacyInternalNamespaces(parser_file);
      Parser.FixMissingInternalNamespaces(parser_file);
      hpp_parsed_files.Add(input_file, parser_file);
      var exported_declarations =
          Parser.CollectExportedDeclarations(parser_file);
      foreach (var exported_declaration in exported_declarations) {
        declaration_to_file.Add(exported_declaration, input_file);
      }
      RewriteFile(input_file, parser_file, dry_run);
    }

    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    var cpp_parsed_files = new Dictionary<FileInfo, Parser.File>();
    foreach (FileInfo input_file in cpp_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = Parser.ParseFile(input_file);
      Parser.FixLegacyInternalNamespaces(parser_file);
      Parser.FixMissingInternalNamespaces(parser_file);
      cpp_parsed_files.Add(input_file, parser_file);
      RewriteFile(input_file, parser_file, dry_run);
    }

    // Process the files in client projects.
    foreach (DirectoryInfo client in clients) {
      FileInfo[] client_hpp_files = client.GetFiles("*.hpp");
      FileInfo[] client_cpp_files = client.GetFiles("*.cpp");
      FileInfo[] all_client_files =
          client_hpp_files.Union(client_cpp_files).ToArray();
      foreach (FileInfo input_file in all_client_files) {
        if (excluded.Contains(input_file.Name)) {
          continue;
        }
        Parser.File parser_file = Parser.ParseFile(input_file);
        Parser.FixInternalUsingDeclarations(parser_file, declaration_to_file);
        RewriteFile(input_file, parser_file, dry_run);
      }
    }
  }
}

} // namespace renamespacer
} // namespace principia
