$ErrorActionPreference = "Stop"

"Usage:"
"make_principia_release mathematician culture date primary_ksp_version"
"                       compatibility_ksp_versions"
""
"    ARGUMENT                         EXAMPLE VALUES"
"    mathematician                    陈景润"
"    culture                          zh-CN"
"    date                             2017111812"
"    primary_ksp_version              1.3.1"
"    compatibility_ksp_versions       @('1.2.1', '1.3.0')"


$mathematician = $args[0]
$culture = new-object CultureInfo $args[1]
$date = $args[2]
$primary_ksp_version = $args[3]
$compatibility_ksp_versions = $args[4]

if ($culture.lcid -eq 0x1000) {
  write-error ("Culture $culture has LCID 0x1000.")
}

if (!($date -match '^\d{10}$')) {
  write-error ("Argument $date should be the UT1 date of the new moon, " +
               "rounded to the nearest hour, formatted yyyyMMddHH.")
  exit 1
}

$msbuild = &".\find_msbuild.ps1"
$7zip = "${Env:ProgramW6432}\7-Zip\7z.exe"
if (!(test-path -path $7zip)) {
  write-error ("Could not find 7-Zip.")
}

$tag = "$date-$mathematician"

$remote = $(
    git remote -v |
    sls 'https://github.com/mockingbirdnest/Principia.git')[0].tostring().split()[0]

git fetch $remote

if ($tag -in $(git tag)) {
  write-error ("Tag $tag already exists.")
}

git checkout "$remote/master"
git tag $tag -m $mathematician

if (test-path .\Release) {
  rm .\Release -recurse -force
}

&$msbuild                           `
    /t:Clean                        `
    /property:Configuration=Release `
    /property:Platform=x64          `
    .\Principia.sln

&$msbuild                              `
    "/t:ksp_plugin;ksp_plugin_adapter" `
    /property:Configuration=Release    `
    /property:Platform=x64             `
    .\Principia.sln

foreach ($ksp_version in $compatibility_ksp_versions) {
  &$msbuild                                              `
      /t:ksp_plugin_adapter                              `
      "/property:Configuration=Release KSP $ksp_version" `
      /property:Platform=x64                             `
      .\Principia.sln
}

# Sanity check: the error message in the adapter should mention all supported
# versions.

if (!(sls ([regex]::escape(
               [text.encoding]::ascii.getstring(
                   [text.encoding]::unicode.getbytes(
                       $ksp_version)))) -encoding ASCII `
      ".\Release\GameData\Principia\ksp_plugin_adapter.dll")) {
  write-error ("Configuration Release does not target $primary_ksp_version.")
}

foreach ($ksp_version in $compatibility_ksp_versions) {
  if (!(sls ([regex]::escape(
                 [text.encoding]::ascii.getstring(
                     [text.encoding]::unicode.getbytes(
                         $ksp_version)))) -encoding ASCII `
        (".\Release\$ksp_version Compatibility\GameData\Principia\" +
             "ksp_plugin_adapter.dll"))) {
    write-error ("Configuration Release KSP $ksp_version does not target " +
                 "$ksp_version.")
  }
}

$lowercase_mathematician = $mathematician.tolower($culture)

&$7zip a "principia $lowercase_mathematician for $primary_ksp_version.zip" `
         ".\Release\GameData"

foreach ($ksp_version in $compatibility_ksp_versions) {
  cp "principia $lowercase_mathematician for $primary_ksp_version.zip" `
     "principia $lowercase_mathematician for $ksp_version.zip"
  &$7zip a "principia $lowercase_mathematician for $ksp_version.zip" `
           ".\Release\$ksp_version Compatibility\GameData"
}

if ($mathematician.contains("TEST")) {
  echo "Successfully built test release $tag.  Run"
  echo "  git tag --delete $tag"
  echo "to clean up."
} else {
  git push $remote --tags
}
