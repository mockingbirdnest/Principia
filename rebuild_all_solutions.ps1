$msbuild = join-path -path (Get-ItemProperty "HKLM:\SOFTWARE\Wow6432Node\Microsoft\VisualStudio\SxS\VS7")."15.0" -childpath "MSBuild\15.0\Bin\msbuild.exe"
$dependencies = @(".\Google\glog\google-glog.sln",
                  ".\Google\googletest\googletest\msvc\2017\gtest.sln",
                  ".\Google\googletest\googlemock\msvc\2017\gmock.sln",
                  ".\Google\protobuf\vsprojects\protobuf.sln",
                  ".\Google\benchmark\msvc\google-benchmark.sln")

function build_solutions($solutions) {
  foreach ($configuration in "Debug", "Release") {
    foreach ($platform in "x64") {
      foreach ($solution in $solutions) {
        &$msbuild /t:"Build" /m /property:VisualStudioVersion=15.0 /property:Configuration=$configuration /property:Platform=$platform $solution
        if (!$?) {
          exit 1
        }
      }
    }
  }
}

build_solutions($dependencies)
build_solutions(".\Principia\Principia.sln")
