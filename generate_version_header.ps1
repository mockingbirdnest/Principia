$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin"
$newversion = (git describe --tags --always --dirty --abbrev=40 --long)
$headerpath = (join-path $solutiondir "base/version.generated.h")

$generateversionheader = {
  $text = [string]::format(
      "`n"                                                                   +
      "#pragma once`n"                                                       +
      "`n"                                                                   +
      "namespace principia {{`n"                                             +
      "namespace base {{`n"                                                  +
      "`n"                                                                   +
      "char const BuildDate[] = `"{0:yyyy'-'MM'-'dd'T'HH':'mm':'ssK}`";`n"  +
      "char const Version[] =`n"                                             +
      "    u8`"{1}`";`n"                                                     +
      "`n"                                                                   +
      "}}  // namespace base`n"                                              +
      "}}  // namespace principia`n",
      (get-date).ToUniversalTime(),
      $newversion)
  [system.io.file]::writealltext(
      $headerpath,
      $text,
      [system.text.encoding]::utf8)
}

if (test-path -path $headerpath) {
  if ([system.io.file]::readalltext($headerpath) `
          -match '(?m)^\s+u8"([^"]+)";$.*') {
    $oldversion = $matches[1]
  }
  if ($oldversion.equals($newversion)) {
    echo "No change to git describe, leaving base/version.generated.h untouched"
  } else {
    echo ("Updating base/version.generated.h, " +
          "version is $newversion (was $oldversion)")
    &$generateversionheader
  }
} else {
  echo "Creating base/version.generated.h, version is $newversion"
  &$generateversionheader
}
