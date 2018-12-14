# GUI

The GUI is implemented the DialogGUI system that was introduced with KSP 1.1 as part of the Unity 5 upgrade.
This offers the following benefits:
* Styling, font size, font quality (using Text Mesh Pro rendering) are consistent with KSP at all UI scales.
* The GUI works well on high DPI screens that require 150-200% scale factors.
* The GUI is a gameobject, and thus not re-rendered unless neccesary, nor do you specifically have to decide when it is renderered.
* The GUI connects to the plugin in an event driven manner.

The alternative is doing full blown artwork design in Unity, but this was not done for several reasons:
* Creating a GUI via Unity is complex, requiring at least one additional C# assembly to be created.
* It requires special treatment after the fact to get the font rendering OK (i.e. use Text Mesh Pro), because compatible art assets are not available to community modders.
* Updating the GUI is not very appealing for people with an emphasis on programming, rather than artwork.

The design has one very explicit constraint, namely seperation of concerns.
Which means that the GUI is decoupled from the sum of data management, services and the principia plugin interface.
This is done for the sake of code readability.

Principia's GUI requirements are on-par with the most advanced DialogGUI examples out there, requiring tabs, scroll bars,
dynamicly growing and shrinking GUI's. The collective work of GUI design, as found on github, inspired the implementation.
The API documentation, like many parts of KSP API documentation, isn't that great, and thus could not be relied on as the 
sole source of information. In the rare case that a workaround is employed, there will be clear comments about this in the GUI code.
