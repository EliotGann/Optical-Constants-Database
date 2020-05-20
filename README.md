# Optical-Constants-Database
Optical Constants database designed for Igor Pro


opticalconstants.ipf goes anywhere in user procedures (alias / shortcut from your repository location works fine)

within "user procedures" create a shortcut of the "Optical Constants" directory, and make sure it is names exactly "Optical Constants" without shortcut or anything else in the shortcut name.  This will allow new .oc files which are synced to the github directory will automatically be available to your Igor Installation.

Usually this is a supporting library which is called from other programs as needed.  However, if you want to use Optical Constants independantly of other packages, it can be loaded by simply opening up your default procedure (ctrl+M) and adding the line
#include "Optical Constants"
to the top of the file, and compiling the macros.

A optical constants menu item will be added under Macros, where you can load and deposit optical constants.
