##############################################################################
# MAKEFILE: compile SPH programs.                                            #
# Copyright (c) 2009, Joshua E. Barnes, Honolulu, Hawaii                     #
##############################################################################

BINFILES = sphcode_e sphcode_u sphcode_ui sphcode_ux

INCFILES = sphcode.h sphdefs.h kdtree.h smooth.h

SRCFILES = sphcode.c sphio.c kdtree.c smooth.c sphload.c sphgrav.c

# Optimization flags; use second version for SunOS.
OPT = -O3
# OPT = -xO3

# Uncomment next line for correct compilation under Linux.
FSF = -ffloat-store

########################################################################
# sphcode_e: entropy is thermodynamic variable.

sphcode_e: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=e \
		FLAGS="-DENTROPY -DGRAVITY"

sphcode_eq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=eq \
		FLAGS="-DENTROPY -DGRAVITY -DQUICKSCAN"

sphcode_ex: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ex \
		FLAGS="-DENTROPY -DEXTGRAV"

########################################################################
# sphcode_ea: entropy is thermodynamic parameter.

sphcode_ea: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ea \
		FLAGS="-DENTROPY -DADIABATIC -DGRAVITY"

sphcode_eaq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=eaq \
		FLAGS="-DENTROPY -DADIABATIC -DGRAVITY -DQUICKSCAN"

sphcode_eax: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=eax \
		FLAGS="-DENTROPY -DADIABATIC -DEXTGRAV"

########################################################################
# sphcode_u: energy is thermodynamic variable.

sphcode_u: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=u \
		FLAGS="-DGRAVITY"

sphcode_uq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uq \
		FLAGS="-DGRAVITY -DQUICKSCAN"

sphcode_ux: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ux \
		FLAGS="-DEXTGRAV"

########################################################################
# sphcode_ui: energy is thermodynamic parameter.

sphcode_ui: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ui \
		FLAGS="-DISOTHERMAL -DGRAVITY"

sphcode_uiq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uiq \
		FLAGS="-DISOTHERMAL -DGRAVITY -DQUICKSCAN"

sphcode_uix: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uix \
		FLAGS="-DISOTHERMAL -DEXTGRAV"

########################################################################
# sphcode_ur: energy is thermodynamic variable, cooling included.

sphcode_ur: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ur \
		FLAGS="-DRADIATING -DGRAVITY"

sphcode_urq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=urq \
		FLAGS="-DRADIATING -DGRAVITY -DQUICKSCAN"

sphcode_urx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=urx \
		FLAGS="-DRADIATING -DEXTGRAV"

########################################################################
# sphcode_uo: energy is thermodynamic variable, cooling from optically
# thin regions included.

sphcode_uo: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uo \
		FLAGS="-DRADIATING -DOPAQUE -DGRAVITY"

sphcode_uoq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uoq \
		FLAGS="-DRADIATING -DOPAQUE -DGRAVITY -DQUICKSCAN"

sphcode_uox: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uox \
		FLAGS="-DRADIATING -DOPAQUE -DEXTGRAV"

########################################################################
# sphcode_ud: energy is thermodynamic variable, cooling from optically
# thin regions and diffusion in optically thick regions included.

sphcode_ud: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ud \
		FLAGS="-DRADIATING -DDIFFUSING -DGRAVITY"

sphcode_udq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=udq \
		FLAGS="-DRADIATING -DDIFFUSING -DGRAVITY -DQUICKSCAN"

sphcode_udx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udx \
		FLAGS="-DRADIATING -DDIFFUSING -DEXTGRAV"

sphcode_udz: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udz \
		FLAGS="-DRADIATING -DDIFFUSING -DNOACCEL"

########################################################################
# sphcode_uc: energy is thermodynamic variable, conduction included.

sphcode_uc: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uc \
		FLAGS="-DCONDUCTING -DGRAVITY"

sphcode_ucq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ucq \
		FLAGS="-DCONDUCTING -DGRAVITY -DQUICKSCAN"

sphcode_ucx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ucx \
		FLAGS="-DCONDUCTING -DEXTGRAV"

########################################################################
# sphcode_us: energy is thermodynamic variable, starformation included.

sphcode_us: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=us \
		FLAGS="-DSTARFORM -DGRAVITY"

sphcode_usq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=usq \
		FLAGS="-DSTARFORM -DGRAVITY -DQUICKSCAN"

sphcode_usx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=usx \
		FLAGS="-DSTARFORM -DEXTGRAV"

########################################################################
# sphcode_uis: energy is thermodynamic parameter, starformation included.

sphcode_uis: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uis \
		FLAGS="-DISOTHERMAL -DSTARFORM -DGRAVITY"

sphcode_uisq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uisq \
		FLAGS="-DISOTHERMAL -DSTARFORM -DGRAVITY -DQUICKSCAN"

sphcode_uiSq: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uiSq \
		FLAGS="-DISOTHERMAL -DSTARFORM -DGRAVITY -DQUICKSCAN \
		       -DCOMPVISC"

sphcode_uisx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uisx \
		FLAGS="-DISOTHERMAL -DSTARFORM -DEXTGRAV"

########################################################################
# sphcode_uir: energy is thermodynamic parameter, recycling included.

sphcode_uir: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uir \
		FLAGS="-DISOTHERMAL -DSTARDEATH -DGRAVITY"

sphcode_uisr: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uisr \
		FLAGS="-DISOTHERMAL -DSTARFORM -DSTARDEATH -DGRAVITY"

########################################################################
# Build sphcode with self-gravity.

code_con_grav: sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	       sphload_$V.o sphgrav_$V.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V \
	   sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	   sphload_$V.o sphgrav_$V.o -lNBody -lClib -lm

########################################################################
# Build sphcode with external gravity.

code_sin_grav: sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V \
	   sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	   -lNBody -lGSP -lClib -lm

########################################################################
# Build individual object files.

sphcode_$V.o: sphcode.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphcode_$V.o -c sphcode.c

sphio_$V.o: sphio.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphio_$V.o -c sphio.c

kdtree_$V.o: kdtree.c kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o kdtree_$V.o -c kdtree.c

smooth_$V.o: smooth.c smooth.h kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FSF) $(FLAGS) \
	    -o smooth_$V.o -c smooth.c

sphload_$V.o: sphload.c sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphload_$V.o -c sphload.c

sphgrav_$V.o: sphgrav.c sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphgrav_$V.o -c sphgrav.c

########################################################################
