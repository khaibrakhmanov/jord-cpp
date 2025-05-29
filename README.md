# jord-cpp
Numerical code to simulate the dynamics of a dust grain in an accretion disk.

## Building and running the code

### Building the project using Makefiles
Two Makefiles are available for building the program: 
- Makefile.unix: building using make in UNIX-based shells (bash in Linux, powershell in Windows).
- Makefile.win: building using nmake Windows shell.

Each Makefile supports building four configurations (build targets). 
- `release`: main version for running simulations taking into account non-stationary coagulation. Optimised for speed. Produces executable *jord_release.exe* in _./bin/_ folder.
- `debug`: version for debugging the code. Full model taking into account non-stationary coagulation is used. Produces executable *jord_debug.exe* in _./bin/_ folder.
- `release_steady`: version for running the simulations within the approximation of steady drift. Optimised for speed. Produces executable *jord_release_steady.exe* in _./bin/_ folder.
- `debug_steady`: version for debugging the code. The approximation of steady drift is used. Produces executable *jord_debug_steady.exe* in _./bin/_ folder.

In addition, there is a standard `clean` target for removing object and executable files.

All commands are run from the ```make_project``` directory. Examples:
```bash
# UNIX-based shells (bash, powershell):
make -f Makefile.unix release
make -f Makefile.unix debug
make -f Makefile.unix release_steady
make -f Makefile.unix debug_steady
make -f Makefile.unix clean
```
```bash
# OS Windows:
nmake /f Makefile.win release
nmake /f Makefile.linux debug
nmake /f Makefile.linux release_steady
nmake /f Makefile.linux debug_steady
nmake /f Makefile.linux clean
```