# NHWAVE_SWAN
This is a coupled wave-current model. Coupling non hydrostatic circulation model NHWAVE with the wave model SWAN.

## Changes to the build environment

On 2019-01-03, the `Makefile` for this project was restructured.  Compiler-specific items were shifted to external files that can be found in the `Makefile.d` directory:

```
$ ls -1 Makefile.d
gfortran.inc
ifort.inc
```

The `Makefile` itself includes a specific compiler-specific file after setting initial values for some of the variables in the build environment.  The person building the program need only alter which compiler-specific file is included to adapt to different compilers (or systems).

Care has also been taken to augment (not overwrite) standard build variables like `CPPFLAGS` and `LDFLAGS`.  This allows environment configuration software (like modules or VALET) to setup those variables and influence the build process as expected.  Even lacking that, the user is free to set `CPPFLAGS` and `LDFLAGS` by hand to exercise the same level of control over the build.

### Less-verbose build output

The build has been made quieter by default, not displaying the commands being issued by make:

```
$ make
[CPP ] sw_swmod1.fsw => sw_swmod1.f
[FC  ] sw_swmod1.f
[CPP ] sw_swmod2.fsw => sw_swmod2.f
[FC  ] sw_swmod2.f
[CPP ] sw_m_fileio.fsw90 => sw_m_fileio.f90
[FC  ] sw_m_fileio.f90
  :
```

To revert to displaying the commands being issued, set the `VERBOSE` variable when issuing the `make` command:

```
$ make VERBOSE=1
mpif90 nspcg/nspcg.f -c -o nspcg.o -O3 -mcmodel=medium -mtune=generic  -I/opt/shared/openmpi/1.6.5-intel2013-sp1/include -I/lustre/work/kirby/sw/hypre/2.9.0b-intel/include
   :
```

### Debug builds

The previous structure of the `Makefile` required coordination in the setting of the `OPT` and `DEBFLGS` variables.  Enabling debugging required the `OPT` flags to be empty and the `DEBFLGS` set, and vice versa for optimized builds.  The updated structure uses the `DEBUG_ENABLE` variable set at the top of the `Makefile` to set `OPT` to optimization versus debugging flags, and only `OPT` is passed to the compiler.

If `DEBUG_ENABLE` is unset or an empty string, an optimized build is performed; if `DEBUG_ENABLE` is set to any non-empty value, a debugging build is performed.  A debug build can also be requested from the command line itself:

```
$ make VERBOSE=1 DEBUG_ENABLE=1
cpp -I/opt/shared/openmpi/1.6.5-intel2013-sp1/include -I/lustre/work/kirby/sw/hypre/2.9.0b-intel/include -P -C -traditional  -DDOUBLE_PRECISION -DPARALLEL -DVFMC -DINTEL sw_swmod1.fsw > sw_swmod1.f
mpif90  sw_swmod1.f -c -check all -g -O0  -I/opt/shared/openmpi/1.6.5-intel2013-sp1/include -I/lustre/work/kirby/sw/hypre/2.9.0b-intel/include
cpp -I/opt/shared/openmpi/1.6.5-intel2013-sp1/include -I/lustre/work/kirby/sw/hypre/2.9.0b-intel/include -P -C -traditional  -DDOUBLE_PRECISION -DPARALLEL -DVFMC -DINTEL sw_swmod2.fsw > sw_swmod2.f
mpif90  sw_swmod2.f -c -check all -g -O0  -I/opt/shared/openmpi/1.6.5-intel2013-sp1/include -I/lustre/work/kirby/sw/hypre/2.9.0b-intel/include
  :
```
