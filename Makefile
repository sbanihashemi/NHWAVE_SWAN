#-----------BEGIN MAKEFILE-------------------------------------------------

            SHELL           = /bin/sh
            # Default flags to the pre-processor:
            DEF_FLAGS       = -P -C -traditional
            # What should we name the target executable produced?
            EXEC            = nhwave_swan
            # Where should we copy the finished executable on a
            # "make install" by the builder?
            INSTALL_DIR     = ../gfortran/bin
            # Set DEBUG_ENABLE to any non-empty value to produce a
            # non-optimized, more easily debugged executable:
            DEBUG_ENABLE    =

#==========================================================================

#
# Code path flags -- these affect the features that will be conditionally
# enabled inside the source when it is pre-processed into Fortran source
# files.  We start with an empty value, then augment with features:
CODE_PATH_FLAGS    =
#
# Uncommenting lines below enables that feature:
#
CODE_PATH_FLAGS    += -DDOUBLE_PRECISION
CODE_PATH_FLAGS    += -DPARALLEL
#CODE_PATH_FLAGS    += -DLANDSLIDE
#CODE_PATH_FLAGS    += -DSALINITY
#CODE_PATH_FLAGS    += -DTEMPERATURE
#CODE_PATH_FLAGS    += -DBUBBLE
#CODE_PATH_FLAGS    += -DSEDIMENT
#CODE_PATH_FLAGS    += -DVEGETATION
#CODE_PATH_FLAGS    += -DSTATIONARY
#CODE_PATH_FLAGS    += -DCOUPLING
#CODE_PATH_FLAGS    += -DBALANCE2D
CODE_PATH_FLAGS    += -DVFMC

ifeq ($(VERBOSE),)
	VERBOSE_PFX=@
else
	VERBOSE_PFX=
endif

#
# Include compiler-specific stuff:
#
include Makefile.d/gfortran.inc

#==========================================================================

#
# FFLAGS is a standard build environment variable that is often set by
# the user or automated environment initialization systems (modules, VALET)
# so we should augment it, not overwrite its existing value:
#
FFLAGS             += $(OPT)

#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
CPPARGS            = $(CPPFLAGS) $(DEF_FLAGS) $(CODE_PATH_FLAGS)

#--------------------------------------------------------------------------
#  Libraries
#--------------------------------------------------------------------------

#
# LDFLAGS is a standard build environment variable that is often set by
# the user or automated environment initialization systems (modules, VALET)
# and our linking rule(s) use that, so we need only mention library names
# here:
#
LIBS               = -lHYPRE

#
# In theory, our environment initialization system (modules, VALET) setup
# CPPFLAGS for us with the necessary header search paths:
#
INCS               = $(CPPFLAGS)

#--------------------------------------------------------------------------
#  NHWAVE-SWAN Source Code.
#--------------------------------------------------------------------------

MODS               = sw_swmod1.F sw_swmod2.F sw_m_fileio.F sw_m_constants.F \
                     sw_serv_xnl4v5.F sw_mod_xnl4v5.F  sw_ocpcre.F sw_ocpids.F \
                     sw_ocpmix.F nw_mod_param.F nw_mod_global.F mod_pass.F \
                     sw_common.F nw_mod_util.F

MAIN               = nhwave.F coupler.F main_pass.F sw_init.F sw_swanmain.F \
                     sw_cycle.F sw_swancom1.F sw_swancom2.F sw_swancom3.F \
                     sw_swancom4.F sw_swancom5.F sw_swanout1.F sw_swanout2.F \
                     sw_swanparll.F sw_swanpre1.F sw_swanpre2.F sw_swanser.F \
                     nw_nhwavemain.F

SRCS               = $(MODS) $(MAIN)

OBJS               = $(SRCS:.F=.o) nspcg.o

#--------------------------------------------------------------------------
#  Linking Directives (also the default rule)
#--------------------------------------------------------------------------

$(EXEC): $(OBJS)
ifneq ($(VERBOSE_PFX),)
	@echo "[LD  ] $(EXEC)"
endif
	$(VERBOSE_PFX) $(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

#--------------------------------------------------------------------------
#  Install the target
#--------------------------------------------------------------------------

install: $(EXEC) $(INSTALL_DIR)
ifneq ($(VERBOSE_PFX),)
	@echo "[CP  ] $(EXEC) $(INSTALL_DIR)"
endif
	$(VERBOSE_PFX) cp $(EXEC) $(INSTALL_DIR)

$(INSTALL_DIR):
	$(VERBOSE_PFX) mkdir -p $(INSTALL_DIR)

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
	$(VERBOSE_PFX) $(RM) $(EXEC) *.o *.mod *.f90 *.f $(CLEAN_ADDITIONAL)

clobber: clean
	$(VERBOSE_PFX) $(RM) Errfile* PRINT* run_* log.txt swaninit nh_sw_dk.qs.* output_nw/* output_sw/*

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .fsw .fsw90 .fnw

%.o: %.fsw
ifneq ($(VERBOSE_PFX),)
	@echo "[CPP ] $*.fsw => $*.f"
endif
	$(VERBOSE_PFX) $(CPP) $(CPPARGS) $*.fsw > $*.f
ifneq ($(VERBOSE_PFX),)
	@echo "[FC  ] $*.f"
endif
	$(VERBOSE_PFX) $(FC)  $*.f -c $(FFLAGS) $(F77FLAGS) $(INCS)

%.o: %.fsw90
ifneq ($(VERBOSE_PFX),)
	@echo "[CPP ] $*.fsw90 => $*.f90"
endif
	$(VERBOSE_PFX) $(CPP) $(CPPARGS) $*.fsw90 > $*.f90
ifneq ($(VERBOSE_PFX),)
	@echo "[FC  ] $*.f90"
endif
	$(VERBOSE_PFX) $(FC)  $*.f90 -c $(FFLAGS) $(F90FLAGS) $(INCS)

%.o: %.fnw
ifneq ($(VERBOSE_PFX),)
	@echo "[CPP ] $*.fnw => $*.f90"
endif
	$(VERBOSE_PFX) $(CPP) $(CPPARGS) $*.fnw > $*.f90
ifneq ($(VERBOSE_PFX),)
	@echo "[FC  ] $*.f90"
endif
	$(VERBOSE_PFX) $(FC)  $*.f90 -c $(FFLAGS) $(F90FLAGS) $(INCS)

#--------------------------------------------------------------------------
#  Stuff specific to individual source files
#--------------------------------------------------------------------------

nspcg.o: nspcg/nspcg.f
ifneq ($(VERBOSE_PFX),)
	@echo "[FC  ] nspcg/nspcg.f"
endif
	$(VERBOSE_PFX) $(FC) $< -c -o $@ $(FFLAGS) $(F77FLAGS) $(INCS)
