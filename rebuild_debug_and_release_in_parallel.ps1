$ErrorActionPreference = "Stop"

$msbuild = &".\find_msbuild.ps1"

&$msbuild /t:Restore Principia.sln
&$msbuild            `
    /t:"Clean;Build" `
    /m               `
    .\debug_and_release_in_parallel.xml

if (!$?) {
  exit 1
}
