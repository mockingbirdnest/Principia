from collections import defaultdict
import os;
import re;

# https://eel.is/c++draft/module.import#nt:module-import-declaration.
# We don’t use module partitions for now.
MODULE_IMPORT_DECLARATION = re.compile(
    r'''
      import \s+ (?P<module_name>[^\s;"<>]+) \s* ;
    | import \s* [<"](?P<header_name>[^"]*)[">] \s* ;
    ''',
    re.VERBOSE)

precompiled_units = []
compiled_units = []
linked_units = []


class TranslationUnit:
  def __init__(self, source_file: str):
    self.source_file = source_file.removeprefix("./")
    self.compiled = f"obj/{self.source_file}.o"
    self.binary = (f"bin/{self.source_file.removesuffix('.cpp')}"
                   if source_file.endswith("_test.cpp") or
                      source_file.rsplit("/", 1)[0] == "main.cpp"
                   else None)
    self.precompiled = None
    self.register()

  def register(self):
    if self.precompiled:
      precompiled_units.append(self)
    if self.compiled:
      compiled_units.append(self)
    if self.binary:
      linked_units.append(self)

  def __eq__(self, value):
    return isinstance(value, TranslationUnit) and self.source_file == value.source_file

  def __hash__(self):
    return self.source_file.__hash__()

  def __str__(self):
    return f"{self.__class__.__name__} {self.source_file}"


class ModuleInterfaceUnit(TranslationUnit):
  def __init__(self, module_name: str):
    self.module_name = module_name
    self.precompiled = f"pmu/{self.module_name}.pcm"
    self.binary = None
    if module_name == "std":
      # Handled in the manually-written Makefile.
      self.source_file = None
      self.compiled = None
    else:
      self.source_file = module_name.removeprefix("principia.").replace(
          ".", "/") + ".ixx"
      self.compiled = f"obj/{self.source_file}.o"
      self.register()

class HeaderUnit(TranslationUnit):
  def __init__(self, header_name: str):
    self.header_name = header_name
    self.source_file = header_name
    self.precompiled = f"phu/{self.source_file}.pcm"
    # Header units are precompiled, but not compiled.
    self.compiled = None
    self.binary = None
    self.register()

def find_all_cpp_files(path="."):
  result = []
  for name in os.listdir(path):
    p = path + "/" + name
    if os.path.isfile(p) and p.endswith(".cpp"):
      result.append(p)
    elif os.path.isdir(p) and p.count("/") < 2:
      result += find_all_cpp_files(p)
  return result

imports = defaultdict(dict)
transitive_dependencies = defaultdict(dict)

def add_transitive_imports(translation_unit: TranslationUnit):
  if translation_unit in imports:
    return
  with open(translation_unit.source_file, encoding="utf8") as f:
    for line in f.readlines():
      import_declaration = re.match(MODULE_IMPORT_DECLARATION, line)
      if import_declaration:
        if import_declaration.group("module_name"):
          imported_unit = ModuleInterfaceUnit(
              import_declaration.group("module_name"))
        else:
          imported_unit = HeaderUnit(import_declaration.group("header_name"))
        imports[translation_unit][imported_unit] = None
        transitive_dependencies[translation_unit][imported_unit] = None
        if imported_unit.source_file and os.path.exists(imported_unit.source_file):
          add_transitive_imports(imported_unit)
          transitive_dependencies[translation_unit].update(transitive_dependencies[imported_unit])

for path in find_all_cpp_files():
  add_transitive_imports(TranslationUnit(path))

for unit in sorted(precompiled_units, key=lambda u: u.source_file):
  if not os.path.exists(unit.source_file):
    # This must be a third-party dependency, let the makefile provide its path.
    continue
  rule = f"{unit.precompiled}: {unit.source_file}"
  for i in transitive_dependencies[unit]:
    rule += f" \\\n\t{i.precompiled}"
  print(rule)

for unit in sorted(compiled_units, key=lambda u: u.source_file):
  rule = f"{unit.compiled}: {unit.precompiled or unit.source_file}"
  for d in transitive_dependencies[unit]:
    if d.precompiled:
      rule += f" \\\n\t{d.precompiled}"
  print(rule)

for unit in sorted(linked_units, key=lambda u: u.source_file):
  if unit.binary:
    rule = f"{unit.binary}: {unit.compiled}"
    for d in transitive_dependencies[unit]:
      if d.compiled:
        rule += f" \\\n\t{d.compiled}"
    print(rule)
