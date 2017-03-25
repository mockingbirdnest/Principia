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

1. Select configuration `Release`, platform `x64`.
1. Remove all existing artifacts: `rm .\Release`.
1. Run `Build > Rebuild Solution`.  This makes sure that the tools and generated code are up-to-date.  It also creates the `journal_test.exe` binary.
1. Select project `ksp_plugin` and run `Build > Profile Guided Optimization > Instrument`.
1. Select project `journal` and run `Project > Set as StartUp Project`
1. With project `journal` selected, run `Project > Properties`.  For configuration `Release` select the `Debugging` page and set `Command Arguments` to `--gtest_filter=PlayerTest.Benchmarks`.
1. Run `Debug > Start Without Debugging`.
1. Answer `Do Not Continue With Build` to the dialog box that asks if you would want to rebuild `ksp_plugin`.
1. Answer `Yes` to the dialog box that informs you that there were build errors.
1. The executable runs.  This takes a long time, 3-4 times longer than without instrumentation.  Have a cup of coffee.
1. Select project `ksp_plugin` and run `Build > Profile Guided Optimization > Optimize`.  This creates the optimized `Principia.dll`.
1. Zip the contents of the `<root>\Principia\Release\GameData` folder into `principia <mathematician>.zip`.
