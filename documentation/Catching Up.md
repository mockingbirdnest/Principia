# Instructions for catching up dependencies

- [ ] Create a tag for the current master.
```powershell
  git checkout master
  git tag master20210521
```
- [ ] Create a new branch for the catch-up:
```powershell
  git checkout master
  git checkout -b Ketchup
```
- [ ] Move master to the last commit at the time of the previous catch-up (the
      one on which we applied our changes):
```powershell
  git checkout master
  git log --oneline
  git reset --hard <commit>
```
- [ ] Pull the changes made to google/master:
```powershell
  git pull google master
```
- [ ] Redo our changes on top of google/master:
```powershell
  git checkout Ketchup
  git rebase master
```
- [ ] Update our master
```powershell
  git checkout master
  git merge Ketchup
  git push -f
```
- [ ] Cleanup and squash some of the commits.
```powershell
  git branch -d Ketchup
```

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