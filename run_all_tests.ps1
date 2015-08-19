$configurations = "Debug","Release"

ForEach ($configuration in $configurations) {
  Get-ChildItem ".\Google\glog\$configuration\*" -Include *test.exe | `
  Foreach-Object {
    & $_.FullName
  }

  Get-ChildItem ".\Google\googletest\msvc\gtest\$configuration\*" -Include *test.exe | `
  Foreach-Object {
    & $_.FullName
  }

  Get-ChildItem ".\Google\googlemock\msvc\2015\$configuration\*" -Include *test.exe | `
  Foreach-Object {
    & $_.FullName
  }

  cd ".\Google\protobuf"
  Get-ChildItem ".\vsprojects\$configuration\*" -Include *test.exe,tests.exe | `
  Foreach-Object {
    & $_.FullName
  }
  cd "..\.."

  Get-ChildItem ".\Google\benchmark\msvc\$configuration\*" -Include *test.exe | `
  Foreach-Object {
    & $_.FullName
  }
}