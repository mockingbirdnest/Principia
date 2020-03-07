# Instructions for releasing Principia

All commands are to be run using the GitHub for Windows git Powershell.

We assume the remote corresponding to https://github.com/mockingbirdnest/Principia.git
is called `la-vache`.  `<Mathematician>` and `<mathematician>` stand for the version
names with the appropriate cases.
- [ ] Make sure that the fingerprints for detecting and fixing the stock system are up-to-date.
- [ ] Check that the release name, release date, and lunation number for the following update, used by the update reminder, are correctly set in the main window.
- [ ] Make sure that the latest runs of the [Azure pipelines](https://dev.azure.com/mockingbirdnest/Principia/_build) for Ubuntu and MacOS are [green](https://www.youtube.com/watch?v=lFeLDc2CzOs&feature=youtu.be&t=61).
- [ ] Run the following command, where:
  - `yyyyMMddHH` is the date of the new moon in UT1,
  - `<language>` is the language code for the language in which `<Mathematician>` is written,
  - `1.x.y` is the primary KSP version (the latest version targeted by the `Release` configuration), and
  - `"1.u.v"…` is the (possibly empty) comma-separated list of supported KSP versions that require a separate build, where each version number is quoted.
```powershell
.\make_principia_release "<Mathematician>" "<language>" "yyyyMMddHH" "1.x.y" @("1.u.v"…)
```
  This creates `<root>\principia <mathematician> for 1.x.y.zip`, as well as `<root>\principia <mathematician> for 1.u.v.zip` as appropriate. 
- [ ] The Azure pipelines will start building the Ubuntu and MacOS releases as soon as the new tag is pushed.  This should take about 1 hour.
- [ ] Download the artifacts produced by each pipeline: go to the run for the tag, click on the published artifact, expand `Principia` and `Release`, then download `principia_Linux-yyyyMMddHH-Mathematician-0-g*.tar.gz` and `principia_Darwin-yyyyMMddHH-Mathematician-0-g*.tar.gz`, respectively.
- [ ] Open the `.tar.gz` files and move the `Linux64` and `MacOS64` directories to the `Principia` directory of the `<root>\principia <mathematician> for 1.x.y.zip` archive, and of the `<root>\principia <mathematician> for 1.u.v.zip` as appropriate.
