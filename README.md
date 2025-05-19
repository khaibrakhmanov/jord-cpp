# jord-cpp
Numerical code to simulate the dynamics of a dust grain in an accretion disk

### Makefiles
Two Makefiles have been created for the program: one for Linux and one for Windows. Each Makefile supports building all four configurations: Release, Debug, Release-Steady, and Debug-Steady.

In addition to the listed targets, there is a standard Clean target for removing object and executable files.

All commands are run from the ```make_project``` directory:
```bash
# OS Linux:
make -f Makefile.linux release
make -f Makefile.linux debug
make -f Makefile.linux release_steady
make -f Makefile.linux debug_steady
make -f Makefile.linux clean
```
```bash
# OS Windows:
nmake /f Makefile.win release
nmake /f Makefile.linux debug
nmake /f Makefile.linux release_steady
nmake /f Makefile.linux debug_steady
nmake /f Makefile.linux clean
```