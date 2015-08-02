#Instructions for building Principia on Windows

These instructions are for Visual Studio 2015, using the git
Powershell provided by [GitHub for Windows](https://windows.github.com/).

We assume a working installation of Kerbal Space Program version 1.0.4 is
found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

This project depends upon:
- the KSP assemblies `Assembly-CSharp.dll` and `Assembly-CSharp-firstpass.dll`,
  found in `<KSP directory>\KSP_Data\Managed`;
- the Unity assembly `UnityEngine.dll`, found in
  `<KSP directory>\KSP_Data\Managed`;
- our [fork](https://github.com/mockingbirdnest/benchmark) of the Google glog
  library;
- our [fork](https://github.com/mockingbirdnest/gmock) of the Google gmock
  library;
- our [fork](https://github.com/mockingbirdnest/benchmark) of the Google
  protobuf library;
- our [fork](https://github.com/mockingbirdnest/benchmark) of the Google
  benchmark library;
- parts of the Chromium codebase (for stack tracing support in glog on Windows),
  *modified according to the instructions below*.

##Installation steps.
###Dowloading Principia.
In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.
###KSP and Unity assemblies.
Copy these assemblies to the directory `<root>\KSP Assemblies`.
###Downloading the Google libraries.
In `<root>\Google`, run the following commands.
```powershell
git clone "https://github.com/mockingbirdnest/glog.git"
git clone "https://github.com/mockingbirdnest/gmock.git"
git clone "https://github.com/mockingbirdnest/protobuf.git"
git clone "https://github.com/mockingbirdnest/benchmark.git"
git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
$GitPromptSettings.RepositoriesInWhichToDisableFileStatus += join-path  (gi -path .).FullName chromium
cd chromium
git config core.sparsecheckout true
copy "..\..\Principia\documentation\setup files\chromium_sparse_checkout.txt" ".git/info/sparse-checkout"
git checkout
copy "..\..\Principia\documentation\setup files\chromium.patch"
git am "chromium.patch"
rm "chromium.patch"
cd ..
```

###Building the Google libraries.
In `<root>\Google`, run the following commands.  The replication of the protobuf
build commands is *not* a mistake.
```powershell
$msbuild = join-path -path (Get-ItemProperty "HKLM:\software\Microsoft\MSBuild\ToolsVersions\14.0")."MSBuildToolsPath" -childpath "msbuild.exe"
&$msbuild /t:Build /m /property:Configuration=Debug .\glog\google-glog.sln
&$msbuild /t:Build /m /property:Configuration=Release .\glog\google-glog.sln
&$msbuild /t:Build /m /property:Configuration=Debug .\gmock\msvc\2010\gmock.sln
&$msbuild /t:Build /m /property:Configuration=Release .\gmock\msvc\2010\gmock.sln
&$msbuild /t:Build /m /property:Configuration=Debug .\protobuf\vsprojects\protobuf.sln
&$msbuild /t:Build /m /property:Configuration=Debug .\protobuf\vsprojects\protobuf.sln
&$msbuild /t:Build /m /property:Configuration=Release .\protobuf\vsprojects\protobuf.sln
&$msbuild /t:Build /m /property:Configuration=Release .\protobuf\vsprojects\protobuf.sln
&$msbuild /t:Build /m /property:Configuration=Debug .\benchmark\msvc\google-benchmark.sln
&$msbuild /t:Build /m /property:Configuration=Release .\benchmark\msvc\google-benchmark.sln
```

###Building Principia.
In `<root>\Principia`, run the following commands.
```powershell
$msbuild = join-path -path (Get-ItemProperty "HKLM:\software\Microsoft\MSBuild\ToolsVersions\14.0")."MSBuildToolsPath" -childpath "msbuild.exe"
&$msbuild /t:Build /m /property:Configuration=Debug .\Principia.sln
&$msbuild /t:Build /m /property:Configuration=Release .\Principia.sln
```
