.. highlight:: rst

Troubleshooting
===============

In recent versions of macOS, the SIP system may prevent python to access DYLD_LIBRARY_PATH environment
variable in sub-shells. This can be an issue in Q-Chem when compiled using mkl appearing the following
error message during execution ::

    dyld: Library not loaded: @rpath/libmkl_intel_lp64.dylib

The workaround without disabling the SIP system is to symbolic link the library *libmkl_intel_lp64.dylib* to ::

    /usr/local/lib

which is used by default to find the needed libraries.


