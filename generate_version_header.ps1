$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin"
[IO.File]::ReadAllText((join-path $solutiondir "base/version.hpp"))  `
    -match '(?m)^\s+"([^"]+)";$.*'
$oldversion = $matches[1]
$newversion =  (git describe --tags --always --dirty --abbrev=40 --long)
if ($oldversion = $newversion) {
  echo "No change to git describe, leaving base/version.hpp unmodified"
} else {
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
      (join-path $solutiondir "base/version.hpp"),
      $text,
      [system.text.encoding]::utf8)
}
