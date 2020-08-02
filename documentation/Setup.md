# Instructions for building Principia on Windows

Before starting, make sure the following are installed on your machine:
* Visual Studio 2019 version 16.4.2 with C# and C++ support;
* Windows SDK 10.0.18362;
* .NET Framework 4.7.2 SDK (and Targeting Pack).

The solution contains a C# project named `coverage_analyser` which requires
the Enterprise edition of Visual Studio to build.  It's only a development
tool for our own use, so if you use a different edition, just do 
`Project > Unload Project` on that project.  You will still be able to build 
and test the mod.

These instructions use the git Powershell provided by [GitHub for Windows](https://windows.github.com/).
We assume a working installation of Kerbal Space Program version 1.8.1 is found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

This project depends upon:
- the KSP assembly `Assembly-CSharp.dll`, found in `<KSP directory>\KSP_x64_Data\Managed`;
- the Unity assemblies `UnityEngine.CoreModule.dll`, `UnityEngine.dll`, `UnityEngine.ImageConversionModule.dll`, `UnityEngine.IMGUIMode.dll`, `UnityEngine.InputLegacyModule.dll`, `UnityEngine.PhysicsModule.dll`, `UnityEngine.TextRenderingModule.dll` and `UnityEngine.UI.dll`, found in
  `<KSP directory>\KSP_x64_Data\Managed`;
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
- our [fork]() of the LLNL zfp library;
- parts of the Chromium codebase (for stack tracing support in glog on Windows),
  *modified according to the instructions below*.

## Installation steps

### Dowloading Principia

In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.

### KSP and Unity assemblies

In order to build for KSP 1.8.1, copy the corresponding KSP 1.8.1 assemblies to `<root>\KSP Assemblies\1.8.1`

### Downloading the dependencies

In `<root>`, run the following commands.
```powershell
mkdir "Google"
mkdir "Third Party"
push-location -path "Google"
git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
$GitPromptSettings.RepositoriesInWhichToDisableFileStatus += join-path  (gi -path .).FullName chromium
push-location -path "chromium"
git config core.sparsecheckout true
copy "..\..\Principia\documentation\setup files\chromium_sparse_checkout.txt" ".git/info/sparse-checkout"
git checkout
copy "..\..\Principia\documentation\setup files\chromium.patch"
git am "chromium.patch"
rm "chromium.patch"
pop-location
pop-location
```
### Building

In `<root>`, run the following command:
```powershell
.\Principia\rebuild_all_solutions.ps1
```

# Instructions for building Principia on Linux and macOS

*Note that the released binaries for Linux and macOS are built using [Azure pipelines](https://dev.azure.com/mockingbirdnest/Principia/_build).  The instructions below are best effort.*

Before starting, make sure the following are installed on your machine:
* Build prerequisites: `build-essential`, `clang`, `libc++-dev`, `libc++abi-dev`, `monodevelop`, `subversion`, and `git`;
* Runtime dependencies: `libc++1`.

## Installing the dependencies

In `<root>/Principia`, run the following command:
```bash
./install_deps.sh
```
This will install and compile all the third-party components that Principia uses.  Don't proceed with the next step unless this step has completed without errors.

## Building Principia

In `<root>/Principia`, run the following command:
```bash
./principia_make.sh
```
If some of the unit tests fail, you may or may not be able to run the resulting version of Principia.
