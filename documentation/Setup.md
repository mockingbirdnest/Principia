# Instructions for building Principia on Windows

Before starting, make sure the following are installed on your machine:
* Visual Studio 2022 version 17.14.9 with C# and C++ support;
* Windows SDK 10.0.26100;
* .NET Framework 4.7.2 SDK (and Targeting Pack).

The solution contains a C# project named `coverage_analyser` which requires
the Enterprise edition of Visual Studio to build.  It's only a development
tool for our own use, so if you use a different edition, just do 
`Project > Unload Project` on that project.  You will still be able to build 
and test the mod.

These instructions use the git Powershell provided by [GitHub for Windows](https://windows.github.com/).
We assume a working installation of Kerbal Space Program version 1.12.5 is found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

This project depends upon:
- the KSP assemblies `Assembly-CSharp.dll` and `Assembly-CSharp-firstpass.dll`, found in `<KSP directory>\KSP_x64_Data\Managed`;
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
- our [fork](https://github.com/mockingbirdnest/zfp) of the LLNL zfp library;
- our [modified excerpts](https://github.com/mockingbirdnest/chromium) of the Chromium codebase (for stack tracing support in glog on Windows).

## Installation steps

### Dowloading Principia

In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.

### KSP and Unity assemblies

In order to build for KSP 1.12.5, copy the corresponding KSP 1.12.5 assemblies to `<root>\KSP Assemblies\1.12.5`

### Building

In `<root>`, run the following command:
```powershell
.\Principia\rebuild_all_solutions.ps1
```

*NOTE*: If the build step fails with error C1060, try changing the `/m` flag of the `msbuild` command in the above script to reduce the parallelism (for instance, `/m:1` or `/m:2`).  It may also be necessary to increase the size of the Windows paging file.  The script is optimized for machines with a sizeable amount of RAM.

# Instructions for building Principia on Linux and macOS

*Note that the released binaries for Linux and macOS are built using [Azure pipelines](https://dev.azure.com/mockingbirdnest/Principia/_build).  The instructions below are best effort.*

## Build Prerequisites
Start by cloning the repository; it will create a directory named `<root>/Principia`:
```bash
git clone https://github.com/mockingbirdnest/Principia.git
cd Principia
```

Before going further, make sure the following are installed on your machine:
### Linux
* Plugin build prerequisites: `unzip` `wget` `binutils` `make` `automake` `libtool` `curl` `cmake`;
* Adapter build prequisites: `msbuild`
* Clang version 20, which can be installed thus:
  ```bash
  wget https://apt.llvm.org/llvm.sh
  chmod +x llvm.sh
  sudo ./llvm.sh 20 all
  ```
### macOS
* Installation prerequisite: `brew`; to install:
  ```bash
  arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```
  The `arch` command ensures that this gets the right version when run on a Mac M1.
* Plugin build prerequisites: `cmake` `autoconf` `automake` `libtool` `python`;  to install:
  ```bash
  arch -x86_64 /usr/local/bin/brew install cmake autoconf automake libtool python
  export PATH="$(brew --prefix python3)/libexec/bin:$PATH"
  ```
* Clang version 20; to install: 
  ```bash
  arch -x86_64 /usr/local/bin/brew install llvm@20
  export PATH="/usr/local/opt/llvm@20/bin:$PATH"
  ```
* Adapter build prequisites: `msbuild`. It is included as part of [Mono](https://www.mono-project.com/download/stable/), and can alternatively be installed directly using a package manager.
You might be able to build it using [Visual Studio for Mac](https://visualstudio.microsoft.com/vs/mac/),
but this is untested.

Note that the resulting binary targets macOS 13 (Ventura); running with older versions of macOS may or may not work.
  
## Installing the dependencies

In `<root>/Principia`, run the following command:
```bash
./install_deps.sh
```
This will install and compile all the third-party components that Principia uses.  Don't proceed with the next step unless this step has completed without errors.

### KSP and Unity assemblies
In order to build the adapter, you will need references to the KSP and Unity assemblies.
Principia expects to find these in a directory adjacent to the principia directory.
For example, in order to build for KSP 1.12.5, the assembly directory should be linked to `../KSP Assemblies/1.12.5` (relative to the Principia directory).

A shell script is provided that will do this automatically on macOS:
```bash
./add_ksp_assemblies_macos.sh
```

## Building Principia

In `<root>/Principia`, run the following command:
```bash
make
```
See the Makefile for more options.
If some of the unit tests fail, you may or may not be able to run the resulting version of Principia.
