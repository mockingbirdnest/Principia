# Support

This page explains where/how to get help with Principia. Please read through the
following guidelines.

## Discord

Questions regarding the usage of Principia are best directed to the [Discord
channel for
Principia](https://discord.com/channels/319857228905447436/480397772248580098)
(here is an [invite](https://discord.gg/Hz8sJx7K8e) if you don't have access already).
There are many knowledgeable users there who can help you with basic questions
and provide guidance on getting started with the mod.  The Discord channel is
*not* suitable for reporting bugs, however; a post like "when I do X it crashes"
will likely not been seen by the developers, and therefore the underlying bug
won’t be fixed.

## Forum

Historically, the [Principia thread](https://forum.kerbalspaceprogram.com/topic/162200-wip181-191-1101-1110%E2%80%932-1122%E2%80%935-principia%E2%80%94version-%E2%80%8E%E2%80%8E%D0%BA%D0%BE%D0%BB%D0%BC%D0%BE%D0%B3%D0%BE%D1%80%D0%BE%D0%B2-released-2024-10-02%E2%80%94n-body-and-extended-body-gravitation)
on the KSP forum was a good place to ask for help with the mod, and it contained
some interesting threads, although it has been less active in recent years.
Unfortunately, as of this writing (October 2024) the KSP forum seems to return
500 more often than not.  It is not recommended to post questions there, even
when it is up.  Use Discord instead.

## GitHub Issues

GitHub Issues should be used for reporting bugs (either crashes or
misbehaviours) and requesting new features or improvements.  GitHub is *not*
suitable for questions on how to use the mod ("how do I go to the Moon from Cape
Canaveral?", for instance).  Please use Discord for those.

### Crashes

There are two kinds of crashes that may be caused by Principia: crashes of 
Principia itself, or crashes of KSP or Unity.

When reporting a crash, please indicate whether the crash is reproducible.
There is no need to include reproduction steps.

#### Principia Crashes

These crashes may manifest themselves as a Windows dialog box (on Windows) or a 
SIGABRT (on Linux).  Or the game may just stop.

Go to `<KSP directory>\glog\Principia`, you should find a `FATAL` file with a
name that matches the time of the crash, for instance
`FATAL.20230920-223637.6656.log`. You should also find an `INFO` file created
shortly before that, for instance `INFO.20230920-223615.6656.log`.  Take these
two files (one is not enough), and
[attach](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/attaching-files)
them to the issue.  (Don't try to extract "interesting" snippets from these
files, we need them in their entirety.)

#### KSP or Unity Crashes

If you don’t have a `FATAL` file, you probably had a crash KSP or Unity that
Principia did not detect (Principia can still be the culprit).  Normally, you
should find a crash directory located in one the directories:
* Windows: either
  `%USERPROFILE%\AppData\LocalLow\Squad\Kerbal Space Program\Crashes` or
 `%USERPROFILE%\AppData\Local\Temp\Squad\Kerbal Space Program\Crashes`
 * macOS: possibly `~/Library/Logs/Unity`?
 * Linux: possibly `~/.config/unity3d/Squad/Kerbal Space Program`?

The crash directory has a name that starts with `Crash` and matches the time of
the crash, for instance `Crash_2024-08-22_143949780`.
[Attach](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/attaching-files)
an archive (`zip`, `7z`, `rar`, etc.) of the entire contents of this folder to
the issue (there should be at least 3 files, named `crash.dmp`,
`error.log`, and `Player.log`).  Also include the `<KSP directory>\KSP.log` and
the `INFO` file ([see above](#principia-crashes)) at the time of the crash.
(Don't try to extract "interesting" snippets from these files, we need them in
their entirety.)

### Unexpected Behaviour and Other Bugs

It is possible for the mod to not crash, but still behave in an unexpected
manner.  When that happens, please provide reproduction steps, possibly
including screenshots or movies showing the behaviour that you think is faulty.
Note that these issues may be caused by interactions between mods that are hard
for us to reproduce, so if you can give us a scenario that works in the stock
game that’s preferable.

### New Features and Improvements

If asking for a new feature, please explain the problem you are trying to solve,
not what you think the solution should be (that’s known as the 
[XY Problem](https://en.wikipedia.org/wiki/XY_problem)).  Most users are not in
a position to evaluate the complexity of implementing a particular solution, so
if you suggest a change that would be extremely onerous to implement, it’s
likely to be rejected, even if a cheap change would solve 80% of the problem.
