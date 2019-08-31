$ErrorActionPreference = "Stop"

$msbuild = &".\find_msbuild.ps1"

&$msbuild                            `
    "/t:Clean;Build"                 `
    /m                               `
    /property:Configuration=Release  `
    /property:Platform=x64           `
    Principia.sln

if (!$?) {
  exit 1
}
