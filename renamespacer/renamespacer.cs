using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace principia {
namespace renamespacer {

class Parser {
  public abstract class Node {
    public Node(int line_number, Node parent) {
      this.line_number = line_number;
      this.parent = parent;
      this.children = new List<Node>();
      if (parent != null) {
        position_in_parent = parent.children.Count;
        parent.children.Add(this);
      }
    }

    public virtual string Cxx() {
      return "--FATAL--";
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
      set {
        must_rewrite_ = value;
      }
    }

    public int line_number = 0;
    public int? last_line_number = null;
    public Node parent = null;
    public int position_in_parent = -1;
    public List<Node> children = null;
    protected bool must_rewrite_ = false;
  }

  public class Class : Node {
    public Class(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Class " + name);
    }

    public string name = null;
  }

  public class File : Node {
    public File(int line_number, FileInfo file_info)
        : base(line_number, parent: null) {
      this.file_info = file_info;
      file_namespace = Regex.Replace(
          file_info.Name, "(_body|_test)?\\.[hc]pp", "");
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "File " + file_info.FullName);
    }

    public FileInfo file_info = null;
    public string file_namespace = null;
  }

  public class Include : Node {
    public Include(int line_number, Node parent, string[] path)
        : base(line_number, parent) {
      this.path = path;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Include " + String.Join(", ", path));
    }

    public string[] path = null;
  }

  public class Namespace : Node {
    public Namespace(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
      if (parent != null &&
          parent.GetType() == typeof(Namespace) &&
          ((Namespace)parent).is_internal) {
        is_internal = true;
      } else {
        is_internal = name.StartsWith("internal");
      }
    }

    public override string Cxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      return "namespace " + name + "{";
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Namespace " + name + (is_internal ? "Internal" : ""));
    }

    public string name = null;
    public bool is_internal = false;
  }

  public class UsingDeclaration : Node {
    public UsingDeclaration(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
      string[] segments = name.Split("::");
      declared_name = segments[segments.Length - 1];

      // Try to figure out if the name was declared in a preceding namespace.
      // This is useful to fixup the internal namespaces.
      if (segments.Length == 2) {
        if (parent.GetType() == typeof(Namespace)) {
          foreach (Node sibling in ((Namespace)parent).children) {
            if (sibling.GetType() == typeof(Namespace) &&
                ((Namespace)sibling).name == segments[0]) {
              declared_in_namespace = (Namespace)sibling;
              break;
            }
          }
        }
      }
    }

    public override string Cxx() {
      Debug.Assert(must_rewrite, "inconsistent rewrite");
      if (declared_in_namespace == null) {
        return "using " + name + ";";
      } else {
        return "using " + declared_in_namespace.name + "::" +
               declared_name + ";";
      }
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "UsingDeclaration " + declared_name +
                        (declared_in_namespace == null ?
                            "" : " from " + declared_in_namespace.name));
    }

    public override bool must_rewrite {
      get => must_rewrite_ ||
             (declared_in_namespace != null &&
              declared_in_namespace.must_rewrite);
    }

    public string name = null;
    public string declared_name = null;
    public Namespace declared_in_namespace = null;
  }

  public class UsingDirective : Node {
    public UsingDirective(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "UsingDirective " + name);
    }

    public string name = null;
  }

  public class Struct : Node {
    public Struct(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "Struct " + name);
    }

    public string name = null;
  }

  public class TypeAlias : Node {
    public TypeAlias(int line_number, Node parent, string name)
        : base(line_number, parent) {
      this.name = name;
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "TypeAlias " + name);
    }

    public string name = null;
  }

  private static bool IsClass(string line) {
    return line.StartsWith("class ");
  }

  private static bool IsClosingNamespace(string line) {
    return line != "}  // namespace" && line.StartsWith("}  // namespace");
  }

  private static bool IsOpeningNamespace(string line) {
    return line != "namespace {" && line.StartsWith("namespace ");
  }

  private static bool IsOwnHeaderInclude(string line,
                                         FileInfo input_file) {
    string own_header = Regex.Replace(
        input_file.Name, "(_body|_test)?\\.[hc]pp", ".hpp");
    return line == "#include \"" + own_header + "\"";
  }

  private static bool IsStruct(string line) {
    return line.StartsWith("struct ");
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

  private static bool IsPrincipiaInclude(string line) {
    // Principia header files end in .hpp.
    return line.StartsWith("#include \"") && line.EndsWith(".hpp\"");
  }

  private static string ParseClass(string line) {
    return Regex.Replace(line.Replace("class ", ""), " .*$", "");
  }

  private static string ParseClosingNamespace(string line) {
    return line.Replace("}  // namespace ", "");
  }

  private static string[] ParseIncludedPath(string line) {
    string path = line.Replace("#include \"", "").Replace(".hpp\"", "");
    return path.Split('/');
  }

  private static string ParseOpeningNamespace(string line) {
    return line.Replace("namespace ", "").Replace(" {", "");
  }

  private static string ParseStruct(string line) {
    return Regex.Replace(line.Replace("struct ", ""), " .*$", "");
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

    StreamReader reader = file_info.OpenText();
    int line_number = 1;
    while (!reader.EndOfStream) {
      string line = reader.ReadLine();
      if (IsPrincipiaInclude(line) && !IsOwnHeaderInclude(line, file_info)) {
        var include = new Include(
            line_number,
            parent: current,
            ParseIncludedPath(line));
      } else if (IsOpeningNamespace(line)) {
        current = new Namespace(
            line_number,
            parent: current,
            ParseOpeningNamespace(line));
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
      } else if (IsClass(line)) {
        var klass = new Class(
            line_number,
            parent: current,
            ParseClass(line));
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
      } else if (IsClosingNamespace(line)) {
        var name = ParseClosingNamespace(line);
        if (current is Namespace ns) {
          Debug.Assert(ns.name == name);
          ns.last_line_number = line_number;
        } else {
          Debug.Assert(false);
        }
        current = current.parent;
      }
      ++line_number;
    }
    file.last_line_number = line_number;
    return file;
  }

  public static List<Node> CollectExportedDeclarations(Node node) {
    var exported_declarations = new List<Node>();
    foreach (Node child in node.children) {
      if (child.GetType() == typeof(Namespace) &&
          !((Namespace)child).is_internal) {
        exported_declarations.AddRange(CollectExportedDeclarations(child));
      } else if (child.GetType() == typeof(Class) ||
                 child.GetType() == typeof(Struct) ||
                 child.GetType() == typeof(TypeAlias) ||
                 child.GetType() == typeof(UsingDeclaration)) {
        exported_declarations.Add(child);
      }
    }
    return exported_declarations;
  }

  public static List<Namespace> FindNamespacesToNormalize(File file,
                                                          Node node) {
    var namespaces_to_normalize = new List<Namespace>();
    foreach (Node child in node.children) {
      if (child is Namespace internal_namespace &&
          internal_namespace.name.StartsWith("internal_")) {
        namespaces_to_normalize.Add(internal_namespace);
      } else if (child is Namespace ns) {
        namespaces_to_normalize.AddRange(
            FindNamespacesToNormalize(file, child));
      }
    }
    return namespaces_to_normalize;
  }

  public static void NormalizeNamespaces(File file) {
    var namespaces_to_normalize = FindNamespacesToNormalize(file, file);
    foreach (Namespace internal_namespace in namespaces_to_normalize) {
      var parent = internal_namespace.parent;
      Debug.Assert(parent is Namespace,
                   "internal namespace not within a namespace");
      var file_namespace = new Namespace(internal_namespace.line_number,
                                         parent: null,
                                         file.file_namespace);
      parent.children[internal_namespace.position_in_parent] = file_namespace;
      file_namespace.parent = parent;
      file_namespace.position_in_parent = internal_namespace.position_in_parent;
      file_namespace.children.Add(internal_namespace);
      file_namespace.last_line_number = internal_namespace.last_line_number;
      file_namespace.must_rewrite = true;
      internal_namespace.parent = file_namespace;
      internal_namespace.position_in_parent = 0;
      internal_namespace.name = "internal";
      internal_namespace.must_rewrite = true;
    }
  }
}

class Renamespacer {
  static void RewriteFile(FileInfo input_file,
                          Parser.File file,
                          bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename = input_file.DirectoryName + "\\" +
                              input_file.Name + ".new";

    Parser.Node parent = file;
    int child_position = 0;
    Parser.Node node = parent.children[child_position];
    int node_line_number = node.line_number;
    StreamReader reader = input_file.OpenText();
    StreamWriter writer = File.CreateText(output_filename);
    int line_number = 1;
    while (!reader.EndOfStream) {
      string line = reader.ReadLine();
      Debug.Assert(node_line_number >= line_number);

      bool has_rewritten = false;
      while (node_line_number == line_number) {
        if (node.must_rewrite) {
          has_rewritten = true;
          writer.WriteLine(node.Cxx());
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
        } else {
          // Moving to the next sibling.
          node = parent.children[child_position];
          node_line_number = node.line_number;
        }
      }

      if (!has_rewritten) {
        writer.WriteLine(line);
      }
      ++line_number;
    }
    writer.Close();
    reader.Close();
    if (!dry_run) {
      File.Move(output_filename, input_filename, true);
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
        string[] split = arg.Split(new []{"--", ":"}, StringSplitOptions.None);
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
    var declaration_to_file = new Dictionary<Parser.Node, FileInfo>();
    foreach (FileInfo input_file in hpp_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file = Parser.ParseFile(input_file);
      Parser.NormalizeNamespaces(parser_file);
      hpp_parsed_files.Add(input_file, parser_file);
      var exported_declarations =
          Parser.CollectExportedDeclarations(parser_file);
      foreach (Parser.Node exported_declaration in exported_declarations) {
        declaration_to_file.Add(exported_declaration, input_file);
      }
      RewriteFile(input_file, parser_file, dry_run);
    }

    //FileInfo[] cpp_files = project.GetFiles("*.cpp");
    //var cpp_parsed_files = new Dictionary<FileInfo, Parser.File>();
    //foreach (FileInfo input_file in cpp_files) {
    //  Parser.File parser_file = Parser.ParseFile(input_file);
    //  cpp_parsed_files.Add(input_file, parser_file);
    //}

    //// Process the files in client projects.
    //foreach (DirectoryInfo client in clients) {
    //  FileInfo[] client_hpp_files = client.GetFiles("*.hpp");
    //  FileInfo[] client_cpp_files = client.GetFiles("*.cpp");
    //  FileInfo[] all_client_files =
    //      client_hpp_files.Union(client_cpp_files).ToArray();
    //  var all_parsed_client_files = new Dictionary<FileInfo, Parser.File>();
    //  foreach (FileInfo input_file in all_client_files) {
    //    all_parsed_client_files.Add(input_file, Parser.ParseFile(input_file));
    //  }
    //}
  }
}

}  // namespace renamespacer
}  // namespace principia
