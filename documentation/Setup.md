# Instructions for building Principia on Windows

These instructions are for Visual Studio 2019 version 16.2.0 preview 1, using the git
Powershell provided by [GitHub for Windows](https://windows.github.com/).

We assume a working installation of Kerbal Space Program version 1.7.0 is
found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

This project depends upon:
- the KSP assembly `Assembly-CSharp.dll`, found in `<KSP directory>\KSP_Data\Managed`;
- the Unity assemblies `UnityEngine.dll` and `UnityEngine.UI.dll`, found in
  `<KSP directory>\KSP_Data\Managed`;
- our [fork](https://github.com/mockingbirdnest/glog) of the Google glog
  library;
- our [fork](https://github.com/mockingbirdnest/googletest) of the Google googletest
  and googlemock libraries;
- our [fork](https://github.com/mockingbirdnest/protobuf) of the Google
  protobuf library;
- our [fork](https://github.com/mockingbirdnest/benchmark) of the Google
  benchmark library;
- our [fork](https://github.com/mockingbirdnest/gipfeli) of the Google gipfeli library;
- our [fork](https://github.com/mockingbirdnest/abseil-cpp) of the Google Abseil C++ library;
- parts of the Chromium codebase (for stack tracing support in glog on Windows),
  *modified according to the instructions below*.

## Installation steps

### Dowloading Principia

In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.

### KSP and Unity assemblies

In order to build for KSP 1.3.1, copy the corresponding KSP 1.3.1 assemblies to `<root>\KSP Assemblies\1.3.1`.  For KSP 1.7.0 copy the corresponding KSP 1.7.0 assemblies to `<root>\KSP Assemblies\1.7.0`

### Downloading the Google libraries

In `<root>\Google`, run the following commands.
```powershell
git clone "https://github.com/mockingbirdnest/glog.git"
git clone "https://github.com/mockingbirdnest/googletest.git"
git clone "https://github.com/mockingbirdnest/protobuf.git"
git clone "https://github.com/mockingbirdnest/benchmark.git"
git clone "https://github.com/mockingbirdnest/gipfeli.git"
git clone "https://github.com/mockingbirdnest/abseil-cpp.git"
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
### Building

In `<root>`, run the following command.
```powershell
.\Principia\rebuild_all_solutions.ps1
```
