# CP-ADCA against SE-SPECK32-K64

This is a script for CP-ADCA against a self-equivalence white-box SPECK implementation. The SE-SPECK is constructed with block size 32 and key size 64. The encodings are affine self-equivalences.

# Experiment Environment
Windows 10

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

# Build

```
$ gcc default_white_box_speck_32_K64.c
```

## Run

```
$ .\a.exe
```