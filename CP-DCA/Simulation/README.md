# CP-DCA simulation against SE-SPECK and IF-SPECK

This is a script for CP-DCA simulation against SE-SPECK and IF-SPECK.

# Experiment Environment
Windows 10

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

cmake version 3.21.1

# Build

```
$ mkdir build
$ cd build
$ cmake -DBUILD_TARGET=win32 -G "MinGW Makefiles" ..
$ mingw32-make
```

## Run

```
$ .\CPDCA.exe
```