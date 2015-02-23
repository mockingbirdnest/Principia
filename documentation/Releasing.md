#Instructions for releasing Principia

All commands are to be run using the GitHub for Windows git Powershell.

We assume the remote corresponding to https://github.com/mockingbirdnest/Principia.git
is called `la-vache`.  `<Mathematician>` and `<mathematician>` stand for the version
names with the appropriate cases.
The following will create and push a tag for the current state of mockingbirdnest/master.
```powershell
$mathematician = "<Mathematician>"
git fetch la-vache
git checkout la-vache/master
git tag ((get-date).ToString("yyyyMMddhh-") + $mathematician) -m $mathematician
git push la-vache --tags
```
After this is done, rebuild the solution for Release, and zip the contents of the
`<root>\Principia\Release\GameData` folder into `principia <mathematician>.zip`.
