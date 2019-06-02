$version = "15.9.4"
$preview = "1"

$vswhere = "${Env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe";
$names = &$vswhere                     `
    -prerelease                           `
    -all                                  `
    -requires Microsoft.Component.MSBuild `
    -property installationName

$msbuildpaths = &$vswhere                 `
    -prerelease                           `
    -all                                  `
    -requires Microsoft.Component.MSBuild `
    -find MSBuild\**\Bin\MSBuild.exe

$i = 0;
foreach ($name in $names) {
  if ($name.startswith("VisualStudioPreview/$version-pre.$preview.")) {
    return $msbuildpaths[$i]
  }
  ++$i
}

write-error(
    "Could not find Visual Studio version $version preview $preview;" +
    " found the following versions:`n$([string]::join("`n", $names))")