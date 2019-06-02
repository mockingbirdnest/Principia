$ErrorActionPreference = "Stop"

$msbuild = &".\find_msbuild.ps1"

&$msbuild                                                                         `
    /t:"Clean;Build"                                                              `
    /m                                                                            `
    /property:Configuration=$configuration /property:Platform=$platform $solution `
    .\debug_and_release_in_parallel.xml
