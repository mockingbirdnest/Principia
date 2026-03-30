# Instructions for catching up dependencies

- [ ] Create a tag for the current `master`:
```powershell
  git checkout master
  git tag master20210521
```
- [ ] Create a new branch for the catch-up:
```powershell
  git checkout master
  git checkout -b Ketchup
```
- [ ] Move `master` to the last commit at the time of the previous catch-up, i.e., the
      one on which we applied our changes; typically this will be `google/master`:
```powershell
  git checkout master
  git log --oneline
  git reset --hard <commit>
```
- [ ] Pull the changes made to `google/master`:
```powershell
  git pull google master
```
- [ ] Redo our changes on top of `google/master`:
```powershell
  git checkout Ketchup
  git rebase master
```
- [ ] At this point, check that the project compiles *all* the files and that the tests run.  More likely than not, some changes will be needed:
    * Update the `props` file for the project (and possibly dependents) to use the right language dialect.
    * Add new projects.  If there is a need to generate GUIDs, use https://guidgenerator.com/.
    * Remove unnecessary projects.
    * Update the contents of the projects.  Some solutions have a handy script `build_projects_helper.ps1` to help generate the XML for Visual Studio.  Make sure that the script gets updated to reflect the latest state of the files, including any handling of "exceptional" files (files that don't match a pattern and require special care).
- [ ] Update our `master`:
```powershell
  git checkout master
  git merge Ketchup
  git push -f
```
- [ ] Squash some of the commits on the `master` branch if needed.  Leave the `Ketchup` branch alone for future reference.

## Notes

In `absl` some of the tests related to the arithmetic of `Duration` fail in mysterious ways:

* `Time.FloorConversion`
* `Time.RoundtripConversion`
* `Time.ToChronoTime`
* `SleepFor.Bounded`
* `Duration.FactoryOverloads`
* `Duration.Range`
* `Duration.AbsoluteValue`
* `Duration.FormatDuration`
* `Duration.ParseDuration`

We should avoid `operator/` and `operator*` on `Duration`.

In `absl` the test `NotificationTest.SanityTest` is flaky (or slop-y).
