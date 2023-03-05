using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {

class Parser {
  public abstract class Node {
    // This links the new node as the last child of its parent.
    protected Node(Node? parent) : this(text: null, parent) {}

    protected Node(string? text, Node? parent) {
      this.parent = parent;
      this.text = text;
      children = new List<Node>();
      if (parent != null) {
        position_in_parent = parent.children.Count;
        parent.children.Add(this);
      }
    }

    public void AddChild(Node child) {
      child.parent = this;
      child.position_in_parent = children.Count;
      children.Add(child);
    }

    public void AddChildren(List<Node> children) {
      foreach (Node child in children) {
        AddChild(child);
      }
    }

    public virtual string Cxx() {
      throw new InvalidOperationException();
    }

    // Writes a single node.
    public abstract void WriteNode(string indent = "");

    public void WriteTree(string indent = "") {
      WriteNode(indent);
      foreach (Node child in children) {
        child.WriteTree(indent + "  ");
      }
    }

    public int position_in_parent { get ; private set; } = -1;

    public string? text { get ; protected set; }

    public virtual bool must_rewrite => text == null;

    public Node? parent = null;
    public List<Node> children;
  }

  public abstract class Declaration : Node {
    protected Declaration(Node parent, string name) : this(
        text: null,
        parent,
        name) {}

    protected Declaration(string? text, Node parent, string name) : base(
        text,
        parent) {
      name_ = name;
    }

    public string name {
      get => name_;
      set {
        name_ = value;
        text = null;
      }
    }

    private string name_;
  }

  public class Class : Declaration {
    public Class(string text, Node parent, string name) : base(
        text,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Class " + name);
    }
  }

  public class Constant : Declaration {
    public Constant(string text, Node parent, string name) : base(
        text,
        parent,
        name) { }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Constant " + name);
    }
  }

  // A placeholder for a deleted node.
  public class File : Node {
    public File(FileInfo file_info) : base(
        parent: null) {
      this.file_info = file_info;
      file_namespace_simple_name = "_" +
                              Regex.Replace(file_info.Name,
                                            @"(_body|_test)?\.[hc]pp",
                                            "");
      project_namespace_full_name = "principia::" + file_info.Directory!.Name;
      file_namespace_full_name = project_namespace_full_name +
                                 "::" +
                                 file_namespace_simple_name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "File " + file_info.FullName);
    }

    public FileInfo file_info { get; }
    public string file_namespace_full_name { get; }
    public string file_namespace_simple_name { get; }
    public string project_namespace_full_name { get; }
  }

  public class Function : Declaration {
    public Function(string text, Node parent, string name) : base(
        text,
        parent,
        name) { }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Function " + name);
    }
  }

  public class Include : Node {
    public Include(string text, Node parent, string[] path) : base(
        text,
        parent) {
      this.path = path;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Include " + string.Join(", ", path));
    }

    public string[] path;
  }

  public class Namespace : Declaration {
    public Namespace(Node parent, string name) :
        this(text: null, parent, name) { }

    public Namespace(string? text, Node parent, string name) : base(
        text,
        parent,
        name) {
      if (parent is Namespace{ is_internal: true } ) {
        is_internal = true;
      } else {
        is_internal = name.StartsWith("internal") ||
                      name == "interface";
      }
    }

    public override string Cxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      return "namespace " + name + " {";
    }

    public string ClosingCxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      return "}  // namespace " + name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent +
                        "Namespace " +
                        name +
                        (is_internal ? "Internal" : ""));
    }

    // The compatibility namespace is identified by its ::.
    public bool is_compatibility_namespace => name.Contains("::");
    public bool is_internal = false;
    public string? closing_text;
  }

  public class Struct : Declaration {
    public Struct(string text, Node parent, string name) : base(
        text,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Struct " + name);
    }
  }

  public class Text : Node {
    public Text(string text, Node parent) : base(text, parent) {
      if (text == null) {
        throw new ArgumentNullException();
      }
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Text (" + text + ")");
    }
  }

  public class TypeAlias : Declaration {
    public TypeAlias(string text, Node parent, string name) : base(
        text,
        parent,
        name) {}

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "TypeAlias " + name);
    }
  }

  public class UsingDeclaration : Declaration {
    public UsingDeclaration(Node parent, string full_name) : this(
        text: null,
        parent,
        full_name) {}

    public UsingDeclaration(string? text, Node parent, string full_name) :
        base(text, parent, full_name.Split("::")[^1]) {
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

    public override string Cxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      if (declared_in_namespace == null) {
        return "using " + full_name + ";";
      } else {
        return "using " + declared_in_namespace.name + "::" + name + ";";
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
        declared_in_namespace is { must_rewrite: true };

    public string full_name;
    public Namespace? declared_in_namespace;
  }

  public class UsingDirective : Node {
    public UsingDirective(Node parent, string ns) :
        this(text: null, parent, ns) {}

    public UsingDirective(string? text, Node parent, string ns) : base(
        text,
        parent) {
      this.ns = ns;
    }

    public override string Cxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      return "using namespace " + ns + ";";
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "UsingDirective " + ns);
    }

    public string ns;
  }

  private static bool IsClass(string line) {
    return Regex.IsMatch(line, @"^class \w+[ ;].*$");
  }

  private static bool IsClosingNamespace(string line) {
    return line != "}  // namespace" && line.StartsWith("}  // namespace");
  }

  private static bool IsConstant(string line) {
    return Regex.IsMatch(line, @"^constexpr +[^ ]+ +\w+ +=.*$") ||
           Regex.IsMatch(line, @"^constexpr +[^ ]+ +\w+;") ||
           Regex.IsMatch(line, @"^inline +constexpr +[^ ]+ +\w+ +=.*$") ||
           Regex.IsMatch(line, @"^inline +constexpr +[^ ]+ +\w+;");
  }

  private static bool IsFunction(string line) {
    // There are type aliases in named_quantities.hpp that have a ( in them.
    // Note that if the return type is very long the function name will be on
    // its own line, in which case the parsing is even more cheesy than usual.
    return !Regex.IsMatch(line, @"^(using)") &&
           (Regex.IsMatch(line, @"^\w.+ [^: ]+\(.*$") ||
            Regex.IsMatch(line, @"^[A-Z][a-z]\w+\(.*$"));
  }

  private static bool IsOpeningNamespace(string line) {
    return line != "namespace {" &&
           line.StartsWith("namespace ") &&
           !Regex.IsMatch(line, @"^namespace \w+ = .*$");
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
    return Regex.IsMatch(line, @"^struct \w+[ ;].*$");
  }

  private static bool IsTypeAlias(string line) {
    return Regex.IsMatch(line, @"^using +\w+ +=.*$");
  }

  private static bool IsUsingDeclaration(string line) {
    return !IsUsingDirective(line) &&
           (Regex.IsMatch(line, @"^using +[\w:]+;$") ||
            Regex.IsMatch(line, @"^using +[\w:]+operator.*$"));
  }

  private static bool IsUsingDirective(string line) {
    return line.StartsWith("using namespace ");
  }

  private static string ParseClass(string line) {
    return Regex.Replace(line.Replace("class ", ""), @"[; ].*$", "");
  }

  private static string ParseClosingNamespace(string line) {
    // This gets the line *with comments* to be able to identify the namespace,
    // but beware the NOLINT my son!
    return Regex.Replace(line.Replace("}  // namespace ", ""),
                         @" +// NOLINT.*$",
                         "");
  }

  private static string ParseConstant(string line) {
    return Regex.Replace(
        Regex.Replace(Regex.Replace(line, @"^inline +", ""),
                      @"^constexpr +[^ ]+ +",
                      ""),
        @" +=.*$|;$",
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
    return Regex.Replace(Regex.Replace(line, @"using +", ""), @" +=.*$", "");
  }

  private static string ParseUsingDeclaration(string line) {
    return Regex.Replace(line, @"using +", "").Replace(";", "");
  }

  private static string ParseUsingDirective(string line) {
    return line.Replace("using namespace ", "").Replace(";", "");
  }

  private static string Uncomment(string line, ref bool in_comment) {
    string uncommented_line = line;
    if (in_comment) {
      int end_of_comment = uncommented_line.IndexOf("*/", 0);
      if (end_of_comment >= 0) {
        uncommented_line = uncommented_line.Substring(end_of_comment + 2);
        in_comment = false;
      } else {
        return "";
      }
    }
    Debug.Assert(!in_comment);
    uncommented_line = Regex.Replace(uncommented_line, @"//.*$", "");
    for (;;) {
      int start_of_comment = uncommented_line.IndexOf("/*", 0);
      if (start_of_comment < 0) {
        break;
      }
      int end_of_comment = uncommented_line.IndexOf("*/", start_of_comment + 2);
      if (end_of_comment > 0) {
        uncommented_line = uncommented_line.Substring(0, start_of_comment) +
                           uncommented_line.Substring(end_of_comment + 2);
      } else {
        uncommented_line = uncommented_line.Substring(0, start_of_comment);
        in_comment = true;
        break;
      }
    }
    // We normally don't have spaces at end of line, but of course we have
    // spaces before comments.
    return uncommented_line.TrimEnd(' ');
  }

  public static File ParseFile(FileInfo file_info) {
    var file = new File(file_info);
    Node current = file;

    using StreamReader reader = file_info.OpenText();
    bool in_comment = false;
    while (!reader.EndOfStream) {
      string line = reader.ReadLine()!;
      string uncommented_line = Uncomment(line, ref in_comment);
      if (IsPrincipiaInclude(uncommented_line) &&
          !IsOwnHeaderInclude(uncommented_line, file_info)) {
        var include = new Include(line,
                                  parent: current,
                                  ParseIncludedPath(uncommented_line));
      } else if (IsOpeningNamespace(uncommented_line)) {
        current = new Namespace(line,
                                parent: current,
                                ParseOpeningNamespace(uncommented_line));
      } else if (IsClosingNamespace(line)) {
        // Must use the raw line here as we use the comments to identify the
        // line.  Don't put funky stuff on the closing namespace line, please.
        var name = ParseClosingNamespace(line);
        if (current is Namespace ns) {
          Debug.Assert(ns.name == name, ns.name, name);
          ns.closing_text = line;
        } else {
          Debug.Assert(false);
        }
        current = current.parent!;
      } else if (IsClass(uncommented_line)) {
        var klass = new Class(line,
                              parent: current,
                              ParseClass(uncommented_line));
      } else if (IsConstant(uncommented_line)) {
        var constant = new Constant(line,
                                    parent: current,
                                    ParseConstant(uncommented_line));
      } else if (IsFunction(uncommented_line)) {
        var function = new Function(line,
                                    parent: current,
                                    ParseFunction(uncommented_line));
      } else if (IsStruct(uncommented_line)) {
        var strukt = new Struct(line,
                                parent: current,
                                ParseStruct(uncommented_line));
      } else if (IsTypeAlias(uncommented_line)) {
        var type_alias = new TypeAlias(line,
                                       parent: current,
                                       ParseTypeAlias(uncommented_line));
      } else if (IsUsingDirective(uncommented_line)) {
        var using_directive = new UsingDirective(
            line,
            parent: current,
            ParseUsingDirective(uncommented_line));
      } else if (IsUsingDeclaration(uncommented_line)) {
        var using_declaration = new UsingDeclaration(
            line,
            parent: current,
            ParseUsingDeclaration(uncommented_line));
      } else {
        var text = new Text(line, parent: current);
      }
    }
    return file;
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

  private static File? FindFileReferencedByUsingDeclaration(
      UsingDeclaration using_declaration,
      Dictionary<Declaration, File> declaration_to_file) {
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
      FindUsingDeclarations(Node node, bool internal_only) {
    var internal_using_declarations = new List<UsingDeclaration>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal || !internal_only) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDeclaration ud) {
              internal_using_declarations.Add(ud);
            }
          }
        }
        internal_using_declarations.AddRange(
            FindUsingDeclarations(child, internal_only));
      }
    }
    return internal_using_declarations;
  }

  private static List<UsingDirective>
      FindUsingDirectives(Node node, bool internal_only) {
    var internal_using_directives = new List<UsingDirective>();
    foreach (Node child in node.children) {
      if (child is Namespace ns) {
        if (ns.is_internal || !internal_only) {
          foreach (Node grandchild in child.children) {
            if (grandchild is UsingDirective ud) {
              internal_using_directives.Add(ud);
            }
          }
        }
        internal_using_directives.AddRange(
            FindUsingDirectives(child, internal_only));
      }
    }
    return internal_using_directives;
  }

  private static Namespace? FindLastOutermostNamespace(File file) {
    Namespace? last_outermost_namespace = null;
    foreach (Node child in file.children) {
      if (child is Namespace ns) {
        last_outermost_namespace = ns;
      }
    }
    return last_outermost_namespace;
  }

  private static List<Namespace> FindLegacyInternalNamespaces(
      Node node) {
    var legacy_internal_namespaces = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace internal_namespace &&
          internal_namespace.name.StartsWith("internal_")) {
        legacy_internal_namespaces.Add(internal_namespace);
      } else if (child is Namespace ns) {
        legacy_internal_namespaces.AddRange(
            FindLegacyInternalNamespaces(child));
      }
    }
    return legacy_internal_namespaces;
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

  public static void FixCompatibilityNamespace(File file) {
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
    Debug.Assert(parent is File, "Last namespace not within a file");
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

  public static void FixFileUsingDirective(File file) {
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
  public static void FixLegacyInternalNamespaces(File file) {
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
  public static void FixUsingDeclarations(
      File file,
      bool internal_only,
      Dictionary<Declaration, File> declaration_to_file) {
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
  public static void FixMissingInternalNamespaces(
      File file,
      bool insert_using_declarations) {
    var innermost_namespaces = FindInnermostNamespaces(file);
    foreach (var ns in innermost_namespaces) {
      if (!ns.is_internal && !ns.is_compatibility_namespace) {
        var nodes_in_ns = ns.children.ToList();
        // Check if there are names to export from this namespace.  If not,
        // don't touch it.  It may be some local reopening.  This produces a
        // deduped and sorted set of symbols.
        var names = new SortedSet<string>();
        foreach (Node n in nodes_in_ns) {
          if (n is Declaration decl) {
            names.Add(decl.name);
          }
        }
        if (names.Count == 0) {
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

  public static void FixUselessInternalNamespace(File file) {
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
}

class Renamespacer {
  static void RewriteFile(FileInfo input_file, Parser.File file, bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";

    using (StreamWriter writer = File.CreateText(output_filename)) {
      RewriteNode(writer, file);
    }
    if (!dry_run) {
      File.Move(output_filename, input_filename, overwrite: true);
    }
  }

  static void RewriteNode(StreamWriter writer, Parser.Node node) {
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

  static bool IsBody(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"_body\.hpp$|\.cpp$");
  }

  static bool IsTest(FileInfo file_info) {
    return Regex.IsMatch(file_info.Name, @"_test\.cpp$");
  }

  // Usage:
  //   renamespacer --project:quantities \
  //                --client:base --client:physics \
  //                --exclude:macros.hpp --dry_run:false
  // This will renamespace quantities and fix the references in the client
  // projects.  The files will be overwritten.
  static void Main(string[] args) {
    // Parse the arguments.
    DirectoryInfo? project = null;
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
    if (project == null) {
      throw new NullReferenceException();
    }

    // Parse all the files in our project.
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] body_hpp_files =
        Array.FindAll(hpp_files, f => Regex.IsMatch(f.Name, @"^.*_body\.hpp$"));
    FileInfo[] all_body_files = body_hpp_files.Union(cpp_files).ToArray();
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
    var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = Parser.ParseFile(input_file);
      file_info_to_file.Add(input_file, parser_file);
    }

    // First collect all the declarations.  We'll them to rewrite the clients.
    // This must happen before we restructure any namespace, as the namespace
    // names will change and would confuse our super-fancy name resolution.
    var declaration_to_file = new Dictionary<Parser.Declaration, Parser.File>();
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name) || IsBody(input_file)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      var exported_declarations =
          Parser.CollectExportedDeclarations(parser_file);
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
        // This file is not in our project, so we didn't parse it yet.
        Parser.File parser_file = Parser.ParseFile(input_file);
        Parser.FixUsingDeclarations(parser_file,
                                    internal_only: !IsBody(input_file),
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
      Parser.FixUsingDeclarations(parser_file,
                                  internal_only: true,
                                  declaration_to_file);
    }

    // Rewrite the namespaces in our project's header files.
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name) || IsBody(input_file)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      Parser.FixLegacyInternalNamespaces(parser_file);
      Parser.FixMissingInternalNamespaces(parser_file,
                                          insert_using_declarations: true);
      Parser.FixCompatibilityNamespace(parser_file);
      RewriteFile(input_file, parser_file, dry_run);
    }

    // Fix the body files in our project.
    foreach (FileInfo input_file in all_body_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[input_file];
      if (IsTest(input_file)) {
        Parser.FixUselessInternalNamespace(parser_file);
        Parser.FixFileUsingDirective(parser_file);
      } else {
        Parser.FixLegacyInternalNamespaces(parser_file);
        Parser.FixMissingInternalNamespaces(parser_file,
                                            insert_using_declarations: false);
      }
      RewriteFile(input_file, parser_file, dry_run);
    }
  }
}

} // namespace renamespacer
} // namespace principia
