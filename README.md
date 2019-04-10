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
  * `libfltk1.3-dev`
* install `libgiac-dev`
  * available from Debian's unstable-repository or https://www-fourier.ujf-grenoble.fr/~parisse/install_en#packages
* compile and install LoAT-purrs by the following steps:
  * `autoreconf --install`
  * `automake`
  * `./configure`
  * `make`
  * `sudo checkinstall`
  * alternatively, `sudo make install` should work as well, but then you bypass your package manager...

## Known Issues

 * Adding _compiler flags_ (e.g. for non-standard include directories):

   PURRS' Makefiles seem to ignore (or override) the usual environment variables (`CFLAGS`, `CXXFLAGS`).
   Instead, one can specify flags to `configure` with the `--with-cxxflags` option (see below for examples).

 * Depending on how `libntl` was configured, it might need `pthread`. So if you see _linker errors_ similar to these:

   ```
   libntl.so: undefined reference to `pthread_key_create'
   libntl.so: undefined reference to `pthread_setspecific'
   ```

   You might have to add the `-pthread` flag to the compiler:

   ```
   ./configure --with-cxxflags='-pthread'
   ```

 * If you get _compiler errors caused by the GiNaC headers_ (ginac.h), you might need to enable C++11 features:

   ```
   ./configure --with-cxxflags='-std=c++11'
   ```

## Demos

To reduce the number of dependencies, the `demos` target is currently disabled in the Makefile.
If you want to compile the demos (see the `demos` subdirectory), uncomment the corresponding `SUBDIRS` line in
the file `Makefile.am`.

