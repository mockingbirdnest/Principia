$OutputEncoding = [console]::InputEncoding = [console]::OutputEncoding =
    New-Object System.Text.UTF8Encoding
$solutiondir = resolve-path $args[0]
$assemblyinfopath = (join-path $solutiondir "ksp_plugin_adapter\Properties\AssemblyInfo.cs")

$datetime = (select-string -path $assemblyinfopath -pattern ".*\((20[0-9]{2}).*([0-9]{1,2}).*([0-9]{1,2}).*")
$datetime.Matches.Length

$assemblyinfotext = [string]::format(
    "using System.Reflection;`n" +
    "using System.Runtime.InteropServices;`n" +
    "[assembly: AssemblyTitle(`"ksp_plugin_adapter`")]`n" +
    "[assembly: AssemblyDescription(`"`")]`n" +
    "[assembly: AssemblyConfiguration(`"`")]`n" +
    "[assembly: AssemblyCompany(`"`")]`n" +
    "[assembly: AssemblyProduct(`"ksp_plugin_adapter`")]`n" +
    "[assembly: AssemblyCopyright(`"Copyright ©  2014`")]`n" +
    "[assembly: AssemblyTrademark("")]`n" +
    "[assembly: AssemblyCulture("")]`n" +
    "[assembly: ComVisible(false)]`n" +
    "[assembly: Guid(`"8e665b76-03f8-473f-89fc-b07d08df7fe6`")]`n" +
    "[assembly: AssemblyVersion(`"`")]`n" +
    "[assembly: AssemblyFileVersion(`"`")]`n",
        $newdate.ToUniversalTime(),
        $newversion)

for(;;) {
  try {
    if ((test-path -path $assemblyinfopath) -and
        [system.io.file]::readalltext($assemblyinfopath).equals($assemblyinfotext)) {
      echo "No change to plugin version, leaving ksp_plugin_adapter\Properties\AssemblyInfo.cs untouched"
      return
    }
    break
  } catch {
    start-sleep -m 10
  }
}

for(;;) {
  try {
    echo "Updating ksp_plugin_adapter\Properties\AssemblyInfo.cs, version is $datetime"
    [system.io.file]::writealltext(
          $assemblyinfopath,
          $assemblyinfotext,
          [system.text.encoding]::utf8)
    break
  } catch {
      start-sleep -m 10
  }
}
