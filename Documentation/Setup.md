#Instructions for building Principia

These instructions are for Visual Studio 2013.

We assume a working, up-to-date installation of Kerbal Space Program is found in
`<KSP directory>`.

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
- The Google [gmock 1.7.0 and gtest 1.7.0](https://code.google.com/p/googlemock
/downloads/list) libraries, *modified according to the instructions below*.

The following instructions should be followed before opening the repository, so
that all dependencies are found.
##KSP and Unity assemblies.
Those assemblies should be copied to the directory `<root>\KSP Assemblies`.
##Google assemblies.
1. Create a git repository (`git init`) in `<root>\Google`.
2. Copy the `.gitignore` and `.gitattributes` from `<root>\Principia` into 
`<root>\Google`, `git add -A` and commit.
3. Download [glog 0.3.3](https://code.google.com/p/google-glog/downloads/list),
and unpack into `<root>\Google`.
There should be a file at `<root>\Google\glog-0.3.3\README` if the unpacking was
done correctly.
4. Open `<root>\Google\glog-0.3.3\google-glog.sln` with Visual Studio 2013. You
will be prompted to upgrade the solution. Do so, save everything, close the
solution, then `git add -A` and commit.
5. Copy `<root>\Principia\Documentation\Setup Files\glog.patch` to
`<root>\Google\glog.patch`, then, in `<root>\Google`, run `git am glog.patch`.
6. Delete `<root>\Google\glog.patch`.
7. Open `<root>\Google\glog-0.3.3\google-glog.sln`, build for Debug and Release.
Close the solution.
8. Download [gmock 1.7.0 and gtest 1.7.0](https://code.google.com/p/googlemock/
downloads/list), and unpack into `<root>\Google`.
There should be a file at `<root>\Google\gmock-1.7.0\README` if the unpacking was
done correctly.
9. Open `<root>\Google\gmock-1.7.0\msvc\2010\gmock.sln` with Visual Studio 2013.
You will be prompted to upgrade, do so. Save everything, `git add -A`, commit.
10. Set the build configuration to Debug.
11. For every project, set the runtime library to Multi-threaded Debug DLL
(`/MDd`) as follows:
  1. Right-click on the project in the solution explorer, then click on
  `Properties`.
  2. In `Configuration Properties -> C/C++ -> Code Generation`, set
  `Runtime Library` to `Multi-threaded Debug DLL (/MDd)`.
12. Set the build configuration to Release.
13. For every project, set the runtime library to Multi-threaded DLL
(`/MD`), similarly to step 11.
14. In both build configurations, add glog to the include path of gmock_main by
prepending `..\..\..\glog-0.3.3\src\windows;` to `Additional Include Directories`
in `gmock_main Properties -> C/C++ -> General`.
15. Save everything, close the solution and commit.
16. Copy `<root>\Principia\Documentation\Setup Files\gmock.patch` to
`<root>\Google\gmock.patch`, then, in `<root>\Google`, run `git am gmock.patch`.
17. Delete `<root>\Google\gmock.patch`.
18. Open `<root>\Google\gmock-1.7.0\msvc\2010\gmock.sln`.
19. Build for Debug and Release.

You are now done with the setup of the dependencies.
Open `<root>\Principia\Principia.sln` and build.
