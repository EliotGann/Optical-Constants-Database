# Optical-Constants-Database
Optical Constants database designed for Igor Pro

within "user procedures" (find your igor pro user file directory) create a shortcut of the "Optical Constants" directory, and make sure it is named exactly "Optical Constants" without shortcut or anything else in the shortcut name.  This will allow new .oc files which are synced to the github directory to automatically be available to your Igor Installation.

Usually this is a supporting library which is called from other programs as needed.  However, if you want to use Optical Constants independantly of other packages, it can be loaded by simply opening up your default procedure (ctrl+M) and adding the line
#include "Optical Constants"
to the top of the file, and compiling the macros.

A optical constants menu item will be added under Macros, where you can load and deposit optical constants.
