* Decide upon an organization name, typically the organisation name on GitHub, say Bar.
* Decide upon a component name, typically the repository name on GitHub, say foo.
* Create directory `Principia\Bar` (capitalized).
* In `Principia\Bar` run `git clone https://github.com/mockingbirdnest/foo.git`.
* In `Principia\Bar\foo` run `mkdir msvc`.
* Create solution `msvc\foo.sln`.
* Copy `common.props` into `msvc\foo.props`.
* `Add > New Project... > Empty Project C++`. name `foo`, created in `foo\msvc`.  This creates `foo\msvc\foo\...` which is good to avoid name clashes between the projects (existing components don't all have that structure).
* If there are benchmarks: `Add > New Project... > Empty Project C++`. name `benchmarks`, created in `foo\msvc`.
* If there are tests: `Add > New Project... > Empty Project C++`. name `tests`, created in `foo\msvc`.  (In complex cases there may be multiple test projects.)
* `Unload Project` for all the projects.  Edit them to read (don't touch the GUIDs):
```
<?xml version="1.0" encoding="utf-8"?>
<Project DefaultTargets="Build" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
  <PropertyGroup Label="Globals">
    <VCProjectVersion>17.0</VCProjectVersion>
    <ProjectGuid>{...}</ProjectGuid>
  </PropertyGroup>
  <Import Project="$(SolutionDir)foo.props" />
  <ItemGroup></ItemGroup>
</Project>
```
* Edit the `msvc\foo.props` file:
  * Replace `the-include-directory` with the name of the directory containing the headers to include.
  * Replace `the-benchmark-project` with the name of the benchmark project(s).  If there is none just use `xxx`.
  * Replace `the-test-project` with the name of the test project(s).  If there is none just use `xxx`.
  * If the benchmark project(s) don't use `Google\benchmark`, remove the imports under condition `PrincipiaBenchmarkProject` at the end of the file.
  * If the test project(s) don't use `Google\googletest`, remove the imports under condition `PrincipiaTestProject` at the end of the file.
* `Load All Projects`.
* Run `Add > Existing Item...` for all header/source files and the various projects.  If there are many, consider writing a script `Bar/foo/build_projects_helper.ps1`.  (See `abseil-cpp` for an example.)
* Make it compile.  This might require to ignore some warnings (change `foo.props`) or change the code.  The symbol `PRINCIPIA` is defined to condition the code changes if needed.
* Make it link.  This requires to add dependencies in the test and benchmark projects to link with the `foo` project's library.
* On the solution, `Save As Solution Filter`, name `msvc\foo.slnf`.
* On `Solution Items`, `Add > Existing Item...`, add the `msvc\foo.slnf` file.
* Edit `msvc\foo.slnf` to remove all the non-production projects (benchmarks, tests, etc.).
* In `Principia\Bar\foo` run `mkdir .github\workflows`.
* Create `build.yaml` in `.github\workflows`.  The simplest is to copy it from another project (e.g., `Google\re2`) and adapt it:
  * Change the `env` variables to have the correct name and to denote the right locations.
  * Change the artifact name to be `foo`.
* In the `build.yaml` file for Principia, add a step to download the `foo` artifact:
```
    - name: Download foo artifact
      uses: mockingbirdnest/actions/windows/download_artifact@main
      with:
        name: foo
        configuration: ${{ matrix.configuration }}
        directory: Bar
```