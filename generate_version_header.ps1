$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin"
$newversion = (git describe --tags --always --dirty --abbrev=40 --long)
$headerpath = (join-path $solutiondir "base/version.hpp")

$generateversionheader = {
  $text = [string]::format(
      "#pragma once`n"                                `
          + "`n"                                      `
          + "namespace principia {{`n"                `
          + "namespace base {{`n"                     `
          + "`n"                                      `
          + "char const kBuildDate[] = `"{0:O}`";`n"  `
          + "char const kVersion[] =`n"               `
          + "    `"{1}`";`n"                          `
          + "`n"                                      `
          + "}}  // namespace base`n"                 `
          + "}}  // namespace principia`n",           `
      (get-date).ToUniversalTime(),
      $newversion)
  [system.io.file]::writealltext(
      $headerpath,
      $text,
      [system.text.encoding]::utf8)
}

if (test-path -path $headerpath) {
  if ([system.io.file]::readalltext($headerpath) `
          -match '(?m)^\s+"([^"]+)";$.*') {
    $oldversion = $matches[1]
  }
  if ($oldversion.equals($newversion)) {
    echo "No change to git describe, leaving base/version.hpp untouched"
  } else {
    echo "Updating base/version.hpp, version is $newversion (was $oldversion)"
    &$generateversionheader
  }
} else {
  echo "Creating base/version.hpp, version is $newversion"
  &$generateversionheader
}
