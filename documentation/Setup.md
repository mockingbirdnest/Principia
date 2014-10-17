#Instructions for building Principia

These instructions are for Visual Studio 2013.

We assume a working installation of Kerbal Space Program version 0.24.2 is
found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.

This project depends upon:
- The KSP assembly `Assembly-CSharp.dll`, found in 
`<KSP directory>\KSP_Data\Managed`;
- The Unity assembly `UnityEngine.dll`, found in
`<KSP directory>\KSP_Data\Managed`;
- The Google [glog 0.3.3](https://code.google.com/p/google-glog/downloads/list)
library, *modified according to the instructions below*;
- The Google [gmock/gtest 1.7.0](https://code.google.com/p/googlemock
/downloads/list) libraries, *modified according to the instructions below*.
- The Google benchmark library for Windows.

The following instructions should be followed before opening the repository, so
that all dependencies are found.
##KSP and Unity assemblies.
Those assemblies should be copied to the directory `<root>\KSP Assemblies`.
##Google assemblies.
1. Download [glog 0.3.3](https://code.google.com/p/google-glog/downloads/list),
  and unpack into `<root>\Google`.
  There should be a file at `<root>\Google\glog-0.3.3\README` if the unpacking
  was done correctly.
2. Download [gmock/gtest 1.7.0]
  (https://code.google.com/p/googlemock/downloads/list), and unpack into
  `<root>\Google`. There should be a file at `<root>\Google\gmock-1.7.0\README`
  if the unpacking was done correctly.
3. In `<root>\Google\glog-0.3.3`, run the following:
```bat
git init
copy "..\..\Principia\.gitattributes" ".gitattributes"
copy "..\..\Principia\.gitignore" ".gitignore"
git add ".gitattributes"
git add ".gitignore"
git commit -m "git files"
git add -A
git commit -m "add glog"
git am "..\..\Principia\documentation\setup files\glog.patch"
```
4. In `<root>\Google\gmock-1.7.0`, run the following:
```bat
git init
copy "..\..\Principia\.gitattributes" ".gitattributes"
copy "..\..\Principia\.gitignore" ".gitignore"
git add ".gitattributes"
git add ".gitignore"
git commit -m "git files"
git add -A
git commit -m "add gmock"
git am "..\..\Principia\documentation\setup files\gmock.patch"
```
4. Open `<root>\Google\glog-0.3.3\google-glog.sln` with Visual Studio 2013.
  Build for Debug and Release. Ignore any warnings. Close the solution.
5. Open `<root>\Google\gmock-1.7.0\msvc\2010\gmock.sln` with Visual Studio
  2013. Build for Debug and Release. Ignore any warnings. Close the solution.
6. In `<root>\Google`, run `git clone https://github.com/pleroy/benchmark.git`.
7. Open `<root>\Google\benchmark\msvc\google-benchmark.sln` with Visual Studio
  2013. Build for Debug and Release. Ignore any warnings. Close the solution.

You are now done with the setup of the dependencies.
Open `<root>\Principia\Principia.sln` and build.
