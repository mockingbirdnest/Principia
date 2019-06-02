$ErrorActionPreference = "Stop"

$msbuild = &".\find_msbuild.ps1"

&$msbuild            `
    /t:"Clean;Build" `
    /m               `
    .\debug_and_release_in_parallel.xml

if (!$?) {
  exit 1
}
