##
## Compiler:          Intel Fortran
##
## This file contains variable definitions that are specific to the Intel
## compiler suite (or an MPI built atop it).
##

# Enable INTEL conditional sections in the source code:
CODE_PATH_FLAGS    += -DINTEL

# Use the system's standard cpp to pre-process the Fortran source:
CPP                = cpp

# Fortan compiler:
ifeq (,$(findstring -DPARALLEL,$(CODE_PATH_FLAGS)))
	FC         = ifort
else
        FC         = mpif90
endif

# Additional flags for F90 code:
F90FLAGS           =

# Additional flags for non-F90 code:
F77FLAGS           =

# Any additional files to remove on clean/clobber?  The .i and .i90
# are interprocedural optimization helpers the compiler generates.
CLEAN_ADDITIONAL   = *.i *.i90

ifeq ($(DEBUG_ENABLE),)
	#
	# Optimization flags tend to be specific to the compiler and the machine
	# on which you're building the code.  We first set the ones that we want
	# on any machine, then augment with machine-specific bits:
	#
	#   -O3	                   perform level-3 code optimizations; compile time
	#                          will be longer and restructuring of the code
	#                          will likely depart from as-written form, so
	#                          the executable will be harder to debug
	#   -mcmodel=medium        make no assumptions about the runtime data
	#                          size (allocatables, etc.) of the program; the
	#                          executable itself is restricted to < 2 GB and
	#                          a more efficient addressing mode is used inside
	#                          the program logic (= marginally faster)
	#
	OPT            = -O3 -mcmodel=medium
	#
	# Mills
	#
	#   -mtune=generic         no processor-specific optimizations
	#
	OPT            += -mtune=generic
	# Farber
	#
	#   -mtune=corei7-avx      enable AVX/SSE SIMD optimizations
	#
	#OPT            += -mtune=corei7-avx
else
	#
	# For debuggining:
	#
	#   -check all             add code to check various conditions at runtime
	#                          (like array bounds, variable initialization)
	#   -g                     generate full debug profile in object files
	#   -O0                    perform no code optimizations
	#
	OPT            = -check all -g -O0
endif
