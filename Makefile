#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            DEF_FLAGS     = -P -C -traditional 
            EXEC          = nhwave_swan
#==========================================================================
#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------

            FLAG_1 = -DDOUBLE_PRECISION 
            FLAG_2 = -DPARALLEL
#            FLAG_3 = -DLANDSLIDE
#           FLAG_4 = -DSALINITY
#           FLAG_5 = -DTEMPERATURE
#             FLAG_6 = -DBUBBLE
#             FLAG_7 = -DSEDIMENT
#             FLAG_8 = -DVEGETATION
             FLAG_9 = -DINTEL
#             FLAG_10 = -DSTATIONARY
#             FLAG_11 = -DCOUPLING
#             FLAG_12 = -DBALANCE2D
             FLAG_13 = -DVFDK
#              
#--------------------------------------------------------------------------
#  mpi defs 
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp 
         CPPFLAGS = $(DEF_FLAGS)
#         FC       = ifort
         FC       = mpif90
         DEBFLGS  = 
         OPT      = #-g
#         OPT      = -O3  -qstrict -qtune=auto -qcache=auto -qalign=4k -w -qfixed
         CLIB     = 
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib
#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
           CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) \
			$(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6) \
			$(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10)  \
			$(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
			$(FLAG_15) $(FLAG_16) $(FLAG_17) $(FLAG_18) \
			$(FLAG_19) $(FLAG_20) $(FLAG_21) $(FLAG_22) \
			$(FLAG_23) $(FLAG_24)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------

            LIBS  = -L/home/1370/hypre-2.8.0b/src/hypre/lib -lHYPRE
            INCS  = -L/home/1370/hypre-2.8.0b/src/hypre/include/

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .fsw .fsw90 .fnw

.fsw.o:
	$(CPP) $(CPPARGS) $*.fsw > $*.f
	$(FC)  $*.f -c $(FFLAGS) $(INCS)

.fsw90.o:
	$(CPP) $(CPPARGS) $*.fsw90 > $*.f90
	$(FC)  $*.f90 -c $(FFLAGS) $(INCS)
.fnw.o:
	$(CPP) $(CPPARGS) $*.fnw > $*.f90
	$(FC)  $*.f90 -c $(FFLAGS) $(INCS)
#	\rm $*.f90
#--------------------------------------------------------------------------
#  NHWAVE-SWAN Source Code.
#--------------------------------------------------------------------------

MODS  = sw_swmod1.F sw_swmod2.F sw_m_fileio.F sw_m_constants.F sw_serv_xnl4v5.F sw_mod_xnl4v5.F  sw_ocpcre.F sw_ocpids.F sw_ocpmix.F nw_mod_param.F nw_mod_global.F mod_pass.F sw_common.F nw_mod_util.F\

MAIN  = nhwave.F coupler.F main_pass.F sw_init.F sw_swanmain.F sw_cycle.F sw_swancom1.F sw_swancom2.F sw_swancom3.F sw_swancom4.F sw_swancom5.F sw_swanout1.F sw_swanout2.F sw_swanparll.F sw_swanpre1.F sw_swanpre2.F sw_swanser.F nw_nhwavemain.F\

SRCS = $(MODS)  $(MAIN)

OBJS = $(SRCS:.F=.o) nspcg.o

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90 *.f

clobber:	clean
		/bin/rm -f  nhwave_swan Errfile* PRINT* run_* log.txt swaninit nh_sw_dk.qs.* output_nw/* output_sw/*






