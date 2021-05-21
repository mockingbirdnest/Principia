# Instructions for catching-up dependents

- [ ] Create a tag for the current master.
- [ ]  Create a new branch for the catch-up:
```powershell
  git checkout master
  git checkout -b Ketchup
```
- [ ] Move master to the last commit at the time of the previous catch-up (the
      one on which we applied our changes):
```powershell
  git checkout master
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
