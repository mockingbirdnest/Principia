$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin"

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
    (git describe --tags --always --dirty --abbrev=40 --long))
[system.io.file]::writealltext(
    (join-path $solutiondir "base/version.hpp"),
    $text,
    [system.text.encoding]::utf8)
