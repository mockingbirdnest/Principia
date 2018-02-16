$msbuild = join-path -path (Get-ItemProperty "HKLM:\software\Microsoft\MSBuild\ToolsVersions\14.0")."MSBuildToolsPath" -childpath "msbuild.exe"
$dependencies = @(".\Google\glog\google-glog.sln",
                  ".\Google\googletest\googletest\msvc\2015\gtest.sln",
                  ".\Google\googletest\googlemock\msvc\2015\gmock.sln",
                  ".\Google\protobuf\vsprojects\protobuf.sln",
                  ".\Google\benchmark\msvc\google-benchmark.sln")

function build_solutions($solutions) {
  foreach ($configuration in "Debug", "Release") {
    foreach ($platform in "x64") {
      foreach ($solution in $solutions) {
        &$msbuild /t:"Clean;Build" /m /property:VisualStudioVersion=14.0 /property:Configuration=$configuration /property:Platform=$platform $solution
        if (!$?) {
          exit 1
        }
      }
    }
  }
}

build_solutions($dependencies)
build_solutions(".\Principia\Principia.sln")
