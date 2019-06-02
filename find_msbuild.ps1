$version = "15.9.4"
$preview = "1"

$vswhere = "${Env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe";
$msbuild = &$vswhere                      `
    -prerelease                           `
    -version $version                     `
    -requires Microsoft.Component.MSBuild `
    -find MSBuild\**\Bin\MSBuild.exe

$name = &$vswhere                         `
    -prerelease                           `
    -version $version                     `
    -requires Microsoft.Component.MSBuild `
    -property installationName

if (!$msbuild -or
    !$name.startswith("VisualStudioPreview/$version-pre.$preview.")) {
  write-error (
      "could not find Visual Studio version $version preview $preview;" +
      " found the following versions:`n" +
      (&$vswhere -prerelease -all -requires Microsoft.Component.MSBuild `
                 -property installationName))
  exit 1
}

$msbuild
