# Instructions for releasing Principia

All commands are to be run using the GitHub for Windows git Powershell.

We assume the remote corresponding to https://github.com/mockingbirdnest/Principia.git
is called `la-vache`.  `<Mathematician>` and `<mathematician>` stand for the version
names with the appropriate cases.
- [ ] Make sure that the fingerprints for detecting and fixing the stock system are up-to-date.
- [ ] If using profile-guided optimization (PGO), make sure that the `PlayerTest.Benchmark` test contains the latest journal file name and has been merged as the last pull request.
- [ ] Check that the release name, release date, and lunation number for the following update, used by the update reminder, are correctly set in the adapter.
- [ ] Make sure that the project can be built and tested under Linux.
- [ ] Create and push a tag for the current state of mockingbirdnest/master, where `"yyyyMMddHH"` is the date of the new moon in UT1:
```powershell
$mathematician = "<Mathematician>"
git fetch la-vache
git checkout la-vache/master
git tag ("yyyyMMddHH-" + $mathematician) -m $mathematician
git push la-vache --tags
```
- [ ] In Visual Studio, select configuration `Release`, platform `x64`.
- [ ] Remove all existing artifacts: `rm .\Release`.
- [ ] Run `Build > Rebuild Solution`.  This makes sure that the tools and generated code are up-to-date.  It also creates the `journal_test.exe` binary.
- [ ] If using PGO, do the following:
  - [ ] Select project `ksp_physics` and run `Build > Profile Guided Optimization > Instrument`.
  - [ ] Select project `journal` and run `Project > Set as StartUp Project`
  - [ ] With project `journal` selected, run `Project > Properties`.  For configuration `Release` select the `Debugging` page and set `Command Arguments` to `--gtest_filter=PlayerTest*Benchmarks --gtest_also_run_disabled_tests`.  Also set `Environment` to `PATH=$(VC_PGO_RunTime_Dir);$(LocalDebuggerEnvironment)`.
  - [ ] Run `Debug > Start Without Debugging`.
  - [ ] Answer `Do Not Continue With Build` to the dialog box that asks if you would want to rebuild `ksp_physics`.
  - [ ] Answer `Yes` to the dialog box that informs you that there were build errors.
  - [ ] The executable runs.  This takes a long time, 3-4 times longer than without instrumentation.  Have a cup of coffee.
  - [ ] Select project `ksp_physics` and run `Build > Profile Guided Optimization > Optimize`.  This creates the optimized `physics.dll`.
- [ ] Zip the contents of the `<root>\Principia\Release\GameData` folder into `principia <mathematician> for 1.x.y.zip`.
- [ ] If building releases for two versions of KSP, do the following:
  - [ ] In Visual Studio, select configuration `Release KSP 1.u.v`, platform `x64`.
  - [ ] Select project `ksp_plugin_adapter` and run `Build > Build ksp_plugin_adapter`
  - [ ] Copy `principia <mathematician> for 1.x.y.zip` into `principia <mathematician> for 1.u.v.zip`.
  - [ ] Add the contents of the `<root>\Principia\Release\1.u.v Compatibility\GameData` to `principia <mathematician> for 1.u.v.zip`, replacing files as necessary.
- [ ] Go to a Linux machine and checkout the newly created tag:
```shell
git fetch la-vache
git checkout $(git tag | tail -1)
```
- [ ] Build the Linux binary:
```shell
make test
make release
```
- [ ] Get the binary from the Linux machine and add the `Linux64` directory to the zip files (next to `x64`).
- [ ] Repeat the previous three steps on a Mac.
