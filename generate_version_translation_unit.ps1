$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin;$env:localappdata\GitHub\Portab~1\mingw32\bin"
$newdate = [DateTime](git log -1 --format=%cd --date=iso-strict)
$newversion = (git describe --tags --always --dirty --abbrev=40 --long)
$headerpath = (join-path $solutiondir "base/version.generated.cc")

$versionheadertext = [string]::format(
    "`n"                                                                   +
    "#include `"base/version.hpp`"`n"                                      +
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
    $newdate.ToUniversalTime(),
    $newversion)

for(;;) {
  try {
    if ((test-path -path $headerpath) -and
        [system.io.file]::readalltext($headerpath).equals($versionheadertext)) {
      echo "No change to git describe, leaving base/version.generated.cc untouched"
      return
    }
    break
  } catch {
    start-sleep -m 10
  }
}

for(;;) {
  try {
    echo "Updating base/version.generated.cc, version is $newversion"
    [system.io.file]::writealltext(
          $headerpath,
          $versionheadertext,
          [system.text.encoding]::utf8)
    break
  } catch {
      start-sleep -m 10
  }
}
