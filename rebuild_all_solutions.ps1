$ErrorActionPreference = "Stop"

$msbuild = &".\Principia\find_msbuild.ps1"

$dependencies = @(".\Google\glog\msvc\glog.sln",
                  ".\Google\googletest\googletest\msvc\gtest.sln",
                  ".\Google\googletest\googlemock\msvc\gmock.sln",
                  ".\Google\protobuf\vsprojects\protobuf.sln",
                  ".\Google\benchmark\msvc\google-benchmark.sln",
                  ".\Google\gipfeli\msvc\gipfeli.sln",
                  ".\Google\abseil-cpp\msvc\abseil-cpp.sln",
                  ".\Inria\core-math\msvc\core-math.sln",
                  ".\LLNL\zfp\msvc\zfp.sln")

foreach ($directory_and_repositories in @(
         @("Boost",  @("config", "multiprecision")),
         @("Google", @("glog", "googletest", "protobuf", "benchmark",
                       "gipfeli", "abseil-cpp", "chromium")),
         @("Inria",  @("core-math")),
         @("LLNL",   @("zfp")))) {
  $directory, $repositories = $directory_and_repositories
  New-Item -ItemType Directory -Force -Path $directory
  Push-Location -Path $directory
  foreach ($repository in $repositories) {
    if (!(Test-Path -Path $repository)) {
      git clone ("https://github.com/mockingbirdnest/" + $repository + ".git")
    }
    Push-Location $repository
    git fetch
    git checkout origin/HEAD
    if (!$?) {
      if ($args[0] -eq "--force") {
        git reset --hard origin/HEAD
        git clean -fdx
      } else {
        Pop-Location
        Pop-Location
        exit 1
      }
    }
    Pop-Location
  }
  Pop-Location
}

function build_solutions($solutions) {
  foreach ($configuration in "Debug", "Release") {
    foreach ($platform in "x64") {
      foreach ($solution in $solutions) {
        &$msbuild /t:"Clean;Build" /m /property:Configuration=$configuration /property:Platform=$platform $solution
        if (!$?) {
          exit 1
        }
      }
    }
  }
}

build_solutions($dependencies)
build_solutions(".\Principia\Principia.sln")
