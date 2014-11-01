$solutiondir = resolve-path $args[0]
$encoding = [system.text.encoding]::utf8
$env:Path += ";$env:localappdata\GitHub\Portab~1\bin"

$text = [string]::format(
    "#pragma once`n"                       `
        + "`n"                             `
        + "namespace principia {{`n"       `
        + "namespace base {{`n"            `
        + "`n"                             `
        + "char const kVersion[] =`n"      `
        + "    `"{0}`";`n"                 `
        + "`n"                             `
        + "}}  // namespace base`n"        `
        + "}}  // namespace principia`n",  `
    (git describe --tags --always --dirty --abbrev=40 --long))
[system.io.file]::writealltext(
    (join-path $solutiondir "base/version.hpp"),
    $text,
    [system.text.encoding]::utf8)
