# Instructions for releasing Principia

All commands are to be run using the GitHub for Windows git Powershell.

We assume the remote corresponding to https://github.com/mockingbirdnest/Principia.git
is called `la-vache`.  `<Mathematician>` and `<mathematician>` stand for the version
names with the appropriate cases.
The following will create and push a tag for the current state of mockingbirdnest/master.
```powershell
$mathematician = "<Mathematician>"
git fetch la-vache
git checkout la-vache/master
git tag ((get-date).ToString("yyyyMMddHH-") + $mathematician) -m $mathematician
git push la-vache --tags
```
After this is done, build the DLL using Profile Guided Optimization:

- [ ] Make sure that the fingerprints for detecting and fixing the stock system are up-to-date.
- [ ] Make sure that the `PlayerTest.Benchmark` test contains the latest journal file name and has been merged as the last pull request.
- [ ] Select configuration `Release`, platform `x64`.
- [ ] Remove all existing artifacts: `rm .\Release`.
- [ ] Run `Build > Rebuild Solution`.  This makes sure that the tools and generated code are up-to-date.  It also creates the `journal_test.exe` binary.
- [ ] Select project `ksp_plugin` and run `Build > Profile Guided Optimization > Instrument`.
- [ ] Select project `journal` and run `Project > Set as StartUp Project`
- [ ] With project `journal` selected, run `Project > Properties`.  For configuration `Release` select the `Debugging` page and set `Command Arguments` to `--gtest_filter=PlayerTest.Benchmarks --gtest_also_run_disabled_tests`.
- [ ] Run `Debug > Start Without Debugging`.
- [ ] Answer `Do Not Continue With Build` to the dialog box that asks if you would want to rebuild `ksp_plugin`.
- [ ] Answer `Yes` to the dialog box that informs you that there were build errors.
- [ ] The executable runs.  This takes a long time, 3-4 times longer than without instrumentation.  Have a cup of coffee.
- [ ] Select project `ksp_plugin` and run `Build > Profile Guided Optimization > Optimize`.  This creates the optimized `Principia.dll`.
- [ ] Zip the contents of the `<root>\Principia\Release\GameData` folder into `principia <mathematician>.zip`.
