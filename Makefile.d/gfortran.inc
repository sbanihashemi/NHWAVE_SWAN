##
## Compiler:          GNU Fortran
##
## This file contains variable definitions that are specific to the GNU
## gfortran compiler (or an MPI built atop it).
##

# Use the system's standard cpp to pre-process the Fortran source:
CPP                = cpp

# Fortan compiler:
ifeq (,$(findstring -DPARALLEL,$(CODE_PATH_FLAGS)))
	FC         = gfortran
else
        FC         = mpif90
endif

# Additional flags for F90 code:
F90FLAGS           = -ffree-form -ffree-line-length-none

# Additional flags for non-F90 code:
F77FLAGS           =

# Any additional files to remove on clean/clobber?
CLEAN_ADDITIONAL   =

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
	#
	OPT            = -O3
else
	#
	# For debuggining:
	#
	#   -fbounds-check         add array bounds checking to the runtime
	#   -fcheck-array-temporaries
	#                          add code to warn when temporary scratch arrays
	#                          must be created to pass a scalar as an array
	#   -g                     generate full debug profile in object files
	#   -O0                    perform no code optimizations
	#
	OPT            = -fbounds-check -fcheck-array-temporaries -g -O0
endif
