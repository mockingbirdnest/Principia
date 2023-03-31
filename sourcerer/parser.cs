using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {

// Parses a file and produces a syntactic tree for it.
public class Parser {

  public static File ParseFile(FileInfo file_info, bool is_body) {
    var file = new File(file_info, is_body);
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
    public File(FileInfo file_info, bool is_body) : base(parent: null) {
      this.file_info = file_info;
      if (is_body) {
        file_namespace_simple_name = "_" +
                                     Regex.Replace(file_info.Name,
                                                   @"(_body|_test)?\.[hc]pp",
                                                   "");
      } else {
        file_namespace_simple_name = "_" +
                                     Regex.Replace(file_info.Name,
                                                   @"\.hpp",
                                                   "");
      }
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
    public Include(string text, Node parent, string[] path) :
        base(text, parent) {
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
        is_internal = name.StartsWith("internal") || name == "interface";
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

    public UsingDeclaration(string? text, Node parent, string full_name) : base(
        text,
        parent,
        full_name.Split("::")[^1]) {
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
    return Regex.Replace(Regex.Replace(Regex.Replace(line, @"^inline +", ""),
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
}

} // namespace renamespacer
} // namespace principia
