# LoAT-purrs

This is a slightly patched version of _The Parma University's Recurrence Relation Solver_ (PURRS), used in the complexity analysis tool LoAT.

The patches mostly involve bugfixes, adaptions to newer library versions, and conversion between PURRS and GiNaC expressions.

Information about PURRS: <http://www.cs.unipr.it/purrs>

## Build

Building LoAT-purrs is a bit painful. The following steps worked for me (on Debian Jessie):

* install the following packages from the Debian repositories:
  * `autoconf`
  * `automake`
  * `libtool`
  * `libginac-dev`
  * `libntl-dev`
  * `libreadline-dev`
* install `libgiac-dev`
  * available from Debian's unstable-repository or https://www-fourier.ujf-grenoble.fr/~parisse/install_en#packages
* comment the line `& !ctrl_c && !interrupted` in `/usr/include/giac/poly.h`, which unfortunately causes compilation errors
* compile and install LoAT-purrs by the following steps:
  * `autoreconf --install`
  * `autoconf` (I guess this step isn't necessary)
  * `automake`
  * `./configure`
  * `make`
  * `sudo checkinstall`
  * alternatively, `sudo make install` should work as well, but then you bypass your package manager...

## Problems

If you get compiler errors caused by the GiNaC headers (ginac.h), you might need to enable C++11 features:

```
./configure --with-cxxflags='-std=c++11'
```

This flag can also be used to specify other compiler flags, such as non standard include directories (via `-I/path/to/include`), since PURRS' Makefiles seem to ignore (or override) the usual environment variables (`CFLAGS`, `CXXFLAGS`).

## Demos

To reduce the number of dependencies, the `demos` target is currently disabled in the Makefile.
If you want to compile the demos (see the `demos` subdirectory), uncomment the corresponding `SUBDIRS` line in
the file `Makefile.am`.

