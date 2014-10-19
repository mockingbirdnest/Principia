#Instructions for building Principia

These instructions are for Visual Studio 2013, using the git
Powershell provided by git for Windows.

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
- Parts of the Chromium codebase (for stack tracing support in glog on Windows),
*modified according to the instructions below*;
- The Google [gmock/gtest 1.7.0](https://code.google.com/p/googlemock
/downloads/list) libraries, *modified according to the instructions below*;
- The Google benchmark library for Windows.

The following instructions should be followed before opening the repository, so
that all dependencies are found.
##KSP and Unity assemblies.
Those assemblies should be copied to the directory `<root>\KSP Assemblies`.
##Google projects.
0. In `<root>\Google`, run

  ```powershell
git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1
$GitPromptSettings.RepositoriesInWhichToDisableFileStatus += join-path  (gi -path .).FullName chromium
cd chromium
git config core.sparsecheckout true
copy "..\..\Principia\documentation\setup files\chromium_sparse_checkout" .git/info/sparse-checkout
git checkout master
copy "..\..\Principia\documentation\setup files\chromium.patch"
git am "chromium.patch"
rm "chromium.patch"
```
0. Download [glog 0.3.3](https://code.google.com/p/google-glog/downloads/list),
  and unpack into `<root>\Google`.
  There should be a file at `<root>\Google\glog-0.3.3\README` if the unpacking
  was done correctly.
0. Download [gmock/gtest 1.7.0]
  (https://code.google.com/p/googlemock/downloads/list), and unpack into
  `<root>\Google`. There should be a file at `<root>\Google\gmock-1.7.0\README`
  if the unpacking was done correctly.
0. In `<root>\Google\glog-0.3.3`, run the following:
  
  ```powershell
git init
copy "..\..\Principia\.gitattributes"
copy "..\..\Principia\.gitignore"
git add ".gitattributes"
git add ".gitignore"
git commit -m "git files"
git add -A
git commit -m "add glog"
copy "..\..\Principia\documentation\setup files\glog.patch"
git am "glog.patch"
rm "glog.patch"
  ```
0. In `<root>\Google\gmock-1.7.0`, run the following:
  
  ```powershell
git init
copy "..\..\Principia\.gitattributes"
copy "..\..\Principia\.gitignore"
git add ".gitattributes"
git add ".gitignore"
git commit -m "git files"
git add -A
git commit -m "add gmock"
copy "..\..\Principia\documentation\setup files\gmock.patch"
git am "gmock.patch"
rm "gmock.patch"
  ```
0. Open `<root>\Google\glog-0.3.3\google-glog.sln` with Visual Studio 2013.
  Build for Debug and Release. Ignore any warnings. Close the solution.
0. Open `<root>\Google\gmock-1.7.0\msvc\2010\gmock.sln` with Visual
  Studio 2013. Build for Debug and Release. Ignore any warnings. Close the
  solution.
0. In `<root>\Google`, run `git clone https://github.com/pleroy/benchmark.git`.
0. Open `<root>\Google\benchmark\msvc\google-benchmark.sln` with Visual
  Studio 2013. Build for Debug and Release. Ignore any warnings. Close the
  solution.

You are now done with the setup of the dependencies.
Open `<root>\Principia\Principia.sln` and build.
