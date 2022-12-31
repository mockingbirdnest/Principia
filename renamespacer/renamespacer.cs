using System;
using System.Collections.Generic;
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
        parent.children.Add(this);
      }
    }

    // Writes a single node.
    public abstract void WriteNode(string indent = "");

    public void WriteTree(string indent = "") {
      WriteNode(indent);
      foreach (Node child in children) {
        child.WriteTree(indent + "  ");
      }
    }

    public int line_number = 0;
    public Node parent = null;
    public List<Node> children = null;
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
    }

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "File " + file_info.FullName);
    }

    public FileInfo file_info = null;
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
      if (parent.GetType() == typeof(Namespace) &&
          ((Namespace)parent).is_internal) {
        is_internal = true;
      } else {
        is_internal = name.StartsWith("internal");
      }
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

    public override void WriteNode(string indent = "") {
      Console.WriteLine(indent + "UsingDeclaration " + declared_name +
                        (declared_in_namespace == null ?
                            "" : " from " + declared_in_namespace.name));
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

  public static File ParseFile(FileInfo input_file) {
    var file = new File(line_number: 0, input_file);
    Node current = file;

    StreamReader reader = input_file.OpenText();
    int line_number = 1;
    while (!reader.EndOfStream) {
      string line = reader.ReadLine();
      if (IsPrincipiaInclude(line) && !IsOwnHeaderInclude(line, input_file)) {
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
        if (current.GetType() != typeof(Namespace) ||
            ((Namespace)current).name != name) {
          Console.WriteLine("Bad closing namespace at line " + line_number +
                            " of " + input_file.Name);
        }
        current = current.parent;
      }
      ++line_number;
    }
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
}

class Renamespacer {
  // Usage:
  //   renamespacer --project:quantities --client:base --client:physics  --move
  // This will renamespace quantities and fix the references in the client
  // projects.  The files will be overwritten.
  static void Main(string[] args) {
    // Parse the arguments.
    DirectoryInfo project = null;
    var clients = new List<DirectoryInfo>();
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
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        }
      }
    }

    // Process the files in the project.
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] all_files = hpp_files.Union(cpp_files).ToArray();
    var all_parsed_files = new Dictionary<FileInfo, Parser.File>();
    foreach (FileInfo input_file in all_files) {
      Parser.File parser_file = Parser.ParseFile(input_file);
      all_parsed_files.Add(input_file, parser_file);
      Console.WriteLine(input_file.FullName);
      foreach (Parser.Node decl in
               Parser.CollectExportedDeclarations(parser_file)) {
        decl.WriteNode("  ");
      }
    }

    // Process the files in client projects.
    foreach (DirectoryInfo client in clients) {
      FileInfo[] client_hpp_files = client.GetFiles("*.hpp");
      FileInfo[] client_cpp_files = client.GetFiles("*.cpp");
      FileInfo[] all_client_files =
          client_hpp_files.Union(client_cpp_files).ToArray();
      var all_parsed_client_files = new Dictionary<FileInfo, Parser.File>();
      foreach (FileInfo input_file in all_client_files) {
        all_parsed_client_files.Add(input_file, Parser.ParseFile(input_file));
      }
    }
  }
}

}  // namespace renamespacer
}  // namespace principia
