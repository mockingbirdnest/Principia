import os
import re

NEITHER_WHITE_SPACE_NOR_SYNTAX = re.compile(
    r"[^\u0009-\u000D \u0085\u200E\u200F\u2028\u2029!-/\:-@\[-\^`\{-~]+"
)

ARITHMETIC_SYMBOLS = set((
    "Cube",
    "Derivative",
    "Difference",
    "Exponentiation",
    "Inverse",
    "Primitive",
    "Product",
    "Quotient",
    "Square",
    "Sum",
))

NAMED_QUANTITIES_SYMBOLS = set((
    "Acceleration",
    "Action",
    "AngularAcceleration",
    "AngularFrequency",
    "AngularMomentum",
    "Area",
    "ArealSpeed",
    "Capacitance",
    "CatalyticActivity",
    "Charge",
    "Concentration",
    "Conductance",
    "CubeRoot",
    "Degree2SphericalHarmonicCoefficient",
    "Degree3SphericalHarmonicCoefficient",
    "Density",
    "DynamicViscosity",
    "ElectricDisplacementField",
    "ElectricField",
    "Energy",
    "Entropy",
    "Force",
    "Frequency",
    "GravitationalParameter",
    "Illuminance",
    "Inductance",
    "Irradiance",
    "Jerk",
    "KinematicViscosity",
    "Luminance",
    "LuminousEfficacy",
    "LuminousEnergy",
    "LuminousExposure",
    "LuminousFlux",
    "MagneticField",
    "MagneticFlux",
    "MagneticFluxDensity",
    "MolarMass",
    "MolarVolume",
    "MomentOfInertia",
    "Momentum",
    "NthRoot",
    "Permeability",
    "Permittivity",
    "Power",
    "Pressure",
    "Radiance",
    "RadiantEnergy",
    "RadiantExposure",
    "RadiantFlux",
    "RadiantIntensity",
    "Resistance",
    "Snap",
    "SolidAngle",
    "SpecificAngularMomentum",
    "SpecificEnergy",
    "SpecificImpulse",
    "SpecificVolume",
    "SpectroscopicWavenumber",
    "Speed",
    "SquareRoot",
    "Stiffness",
    "Torque",
    "Variation",
    "Voltage",
    "Volume",
    "Wavenumber",
))

for bodies in False, True:
  for project in os.listdir():
    if not os.path.isdir(project):
      continue
    for file in os.listdir(project):
      file = f"{project}/{file}"
      if not file.endswith((".hpp", ".cpp")):
        continue
      if bodies != (file.endswith(".cpp") or
                    (file.endswith("_body.hpp") and
                     os.path.exists(file.removesuffix("_body.hpp")))):
        continue
      arithmetic_used_in_header = False
      named_quantities_used_in_header = False
      if bodies:
        header = file.removesuffix("_body.hpp").removesuffix(".cpp") + ".hpp"
        if os.path.exists(header):
          with open(header, encoding="utf8") as f:
            lines = f.readlines()
          for i, line in enumerate(lines):
            line = line.split("//", maxsplit=1)[0]
            identifiers : list[str] = NEITHER_WHITE_SPACE_NOR_SYNTAX.findall(line)
            named_quantities_used_in_header |= identifiers == [
              "using", "namespace", "principia", "quantities", "_arithmetic"]
            named_quantities_used_in_header |= identifiers == [
              "using", "namespace", "principia", "quantities", "_named_quantities"]
      if file in ("quantities/arithmetic.hpp",
                  "quantities/named_quantities.hpp"):
        continue
      with open(file, encoding="utf8") as f:
        lines = f.readlines()
      using_named_quantities_line = None
      arithmetic_symbols_used : set[str] = set()
      named_quantities_symbols_used : set[str] = set()
      for i, line in enumerate(lines):
        line = line.split("//", maxsplit=1)[0]
        identifiers : list[str] = NEITHER_WHITE_SPACE_NOR_SYNTAX.findall(line)
        arithmetic_symbols_used.update(
            identifier for identifier in identifiers
            if identifier in ARITHMETIC_SYMBOLS)
        named_quantities_symbols_used.update(
            identifier for identifier in identifiers
            if identifier in NAMED_QUANTITIES_SYMBOLS)
        if identifiers == [
          "using", "namespace", "principia", "quantities", "_named_quantities"]:
          using_named_quantities_line = i
      if ((arithmetic_symbols_used or
          named_quantities_symbols_used) and
          (not named_quantities_used_in_header and not arithmetic_used_in_header) and
          using_named_quantities_line is None):
        print("---", file, "uses", arithmetic_symbols_used or "",
              named_quantities_symbols_used or "",
              "but does not inherit nor have a using directive for "
              "_named_quantities. These may be homonyms in different "
              "namespaces.")
      if (not arithmetic_symbols_used and
          not named_quantities_symbols_used and
          using_named_quantities_line is not None):
        print("!!!", file,
              "has an unnecessary using directive for _named_quantities")
      if (not named_quantities_symbols_used and
          using_named_quantities_line is not None):
        print("    using directive for _named_quantities will be removed from",
              file)
        del lines[using_named_quantities_line]
      if (arithmetic_symbols_used and
          not arithmetic_used_in_header and
          using_named_quantities_line is not None):
        print("    using directive for _arithmetic       will be added to",
              file)
        lines.insert(using_named_quantities_line,
                     "using namespace principia::quantities::_arithmetic;\n")
      with open(file, mode="w", encoding="utf8") as f:
        f.writelines(lines)
