cd $args[0]
$env:Path += ";$env:localappdata\GitHub\Portab~1\bin"
[string]::format("#pragma once`n#define PRINCIPIA_VERSION_ID ""{0}""`n",
                 (git describe --tags --always --dirty --abbrev=40 --long)) > "base/version.hpp"