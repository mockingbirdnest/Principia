$version = "16.2.2"
$preview = ""
if ($preview.length -gt 0) {
  $description = "version $version preview $preview"
  $path = "VisualStudioPreview/$version-pre.$preview."
} else {
  $description = "version $version"
  $path = "VisualStudio/$version+"
}

$vswhere = "${Env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$names = &$vswhere                        `
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
  if ($name.startswith("$path")) {
    return ($msbuildpaths | select-object -index $i)
  }
  ++$i
}

write-error(
    "Could not find Visual Studio $description;" +
    " found the following versions:`n$([string]::join("`n", $names))")
