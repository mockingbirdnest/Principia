# Instructions for releasing Principia

All commands are to be run using the GitHub for Windows git Powershell.

We assume the remote corresponding to https://github.com/mockingbirdnest/Principia.git
is called `la-vache`.  `<Mathematician>` and `<mathematician>` stand for the version
names with the appropriate cases.
- [ ] Make sure that the fingerprints for detecting and fixing the stock system are up-to-date.
- [ ] Check that the release name, release date, and lunation number for the following update, used by the update reminder, are correctly set in the main window.
- [ ] Make sure that the project can be built and tested under Linux.
- [ ] Run the following command, where:
  - `yyyyMMddHH` is the date of the new moon in UT1,
  - `<language>` is the language code for the language in which `<Mathematician>` is written,
  - `1.x.y` is the primary KSP version (the latest version targeted by the `Release` configuration), and
  - `"1.u.v"…` is the (possibly empty) comma-separated list of supported KSP versions that require a separate build, where each version number is quoted.
```powershell
.\make_principia_release "<Mathematician>" "<language>" "yyyyMMddHH" "1.x.y" @("1.u.v"…)
```
  This creates `<root>\principia <mathematician> for 1.x.y.zip`, as well as `<root>\principia <mathematician> for 1.u.v.zip` as appropriate. 
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
