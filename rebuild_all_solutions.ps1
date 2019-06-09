$ErrorActionPreference = "Stop"

$msbuild = &".\Principia\find_msbuild.ps1"

$dependencies = @(".\Google\glog\google-glog.sln",
                  ".\Google\googletest\googletest\msvc\gtest.sln",
                  ".\Google\googletest\googlemock\msvc\gmock.sln",
                  ".\Google\protobuf\vsprojects\protobuf.sln",
                  ".\Google\benchmark\msvc\google-benchmark.sln",
                  ".\Google\gipfeli\msvc\gipfeli.sln",
                  ".\Google\abseil-cpp\msvc\abseil-cpp.sln")

push-location -path "Google"

foreach ($repository in @("glog", "googletest", "protobuf", "benchmark",
                          "gipfeli", "abseil-cpp")) {
  if (!(test-path -path $repository)) {
    git clone ("https://github.com/mockingbirdnest/" + $repository + ".git")
  }
  push-location $repository
  git checkout master
  git pull
  if (!$?) {
    if ($args[0] -eq "--force") {
      git reset --hard origin/master
      git clean -fdx
    } else {
      pop-location
      pop-location
      exit 1
    }
  }
  pop-location
}
pop-location

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
