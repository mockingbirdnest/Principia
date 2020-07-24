param([switch]$strict = $false)

$version = "16.6.2"
$preview = ""

if ($preview.length -gt 0) {
  $description = "version $version preview $preview"
  $path = "VisualStudioPreview/$version-pre.$preview.0+"
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

function version-tuple($name) {
  $tuple = [int[]]$name.split(@("/", "+"))[1].split([string[]]@(".", "-pre"), [StringSplitOptions]::none)
  if ($tuple.length -lt 5) {
    # Count non-previews as preview ∞.0.
    $tuple = $tuple + @((1.0 / 0.0), 0)
  }
  return $tuple
}

if ($strict) {
  write-error(
      "Could not find Visual Studio $description;" +
      " found the following versions:`n$([string]::join("`n", $names))")
} else {
  $earlier = $null
  $earlier_index = $null
  $later = $null
  $later_index = $null
  $i = 0
  foreach ($name in $names) {
    if ((version-tuple($name) -lt version-tuple($path)) -and
        (($earlier -eq $null) -or
         (version-tuple($name) -gt version-tuple($earlier)))) {
      $earlier = $name
      $earlier_index = $i
    }
    if ((version-tuple($name) -gt version-tuple($path)) -and
        (($later -eq $null) -or
         (version-tuple($name) -lt version-tuple($later)))) {
      $later = $name
      $later_index = $i
    }
    ++$i
  }
  if ($later -ne $null) {
    $best_match = $later
    $i = $later_index
  } elseif ($earlier -ne $null) {
    $best_match = $earlier
    $i = $earlier_index
  } else {
    write-error("Could not find Visual Studio")
    exit 1
  }
  write-warning(
      "Could not find Visual Studio $description;" +
      " falling back to $best_match from:`n$([string]::join("`n", $names))")
  return ($msbuildpaths | select-object -index $i)
}
