$solutiondir = resolve-path $args[0]
$env:Path += ";$env:localappdata\GitHub\Portab~1\bin"
[string]::format(
    "#pragma once`n"                                             `
        + "`n"                                                   `
        + "namespace principia {{`n"                             `
        + "namespace base {{`n"                                  `
        + "`n"                                                   `
        + "char const kVersion[] = `"{0}`";`n"                   `
        + "`n"                                                   `
        + "}}  // namespace base`n"                              `
        + "}}  // namespace principia`n",                        `
    (git describe --tags --always --dirty --abbrev=40 --long)) > `
        (join-path $solutiondir "base/version.hpp")
