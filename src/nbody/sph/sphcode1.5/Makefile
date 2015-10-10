##############################################################################
# MAKEFILE: compile SPH programs.                                            #
# Copyright (c) 2012, Joshua E. Barnes, Honolulu, Hawaii                     #
##############################################################################

BINFILES = sphcode_e sphcode_ea sphcode_u sphcode_ui sphcode_ux\
	   sphcode_ur sphcode_uo sphcode_ud sphcode_uc sphcode_us\
	   sphcode_uis sphcode_uir sphdensity

INCFILES = sphcode.h sphdefs.h kdtree.h smooth.h

SRCFILES = sphcode.c sphio.c kdtree.c smooth.c sphload.c sphgrav.c

# Optimization flags; use second version for SunOS.
OPT = -O3
# OPT = -xO3

# Uncomment next line for correct compilation under Linux.
FSF = -ffloat-store

make_all: $(BINFILES)

########################################################################
# sphcode_e: entropy is thermodynamic variable.

sphcode_e: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=e \
		FLAGS="-DENTROPY -DGRAVITY"

sphcode_ex: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ex \
		FLAGS="-DENTROPY -DEXTGRAV"

########################################################################
# sphcode_ea: entropy is thermodynamic parameter.

sphcode_ea: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ea \
		FLAGS="-DENTROPY -DADIABATIC -DGRAVITY"

sphcode_eax: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=eax \
		FLAGS="-DENTROPY -DADIABATIC -DEXTGRAV"

########################################################################
# sphcode_u: energy is thermodynamic variable.

sphcode_u: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=u \
		FLAGS="-DGRAVITY -DCOMPVISC"

sphcode_ux: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ux \
		FLAGS="-DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_ui: energy is thermodynamic parameter.

sphcode_ui: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ui \
		FLAGS="-DISOTHERMAL -DGRAVITY -DCOMPVISC"

sphcode_uix: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uix \
		FLAGS="-DISOTHERMAL -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_ur: energy is thermodynamic variable, cooling included.

sphcode_ur: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ur \
		FLAGS="-DRADIATING -DGRAVITY -DCOMPVISC"

sphcode_urx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=urx \
		FLAGS="-DRADIATING -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_uo: energy is thermodynamic variable, cooling from optically
# thin regions included.

sphcode_uo: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uo \
		FLAGS="-DRADIATING -DOPAQUE -DGRAVITY -DCOMPVISC"

sphcode_uox: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uox \
		FLAGS="-DRADIATING -DOPAQUE -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_ud: energy is thermodynamic variable, cooling from optically
# thin regions and diffusion in optically thick regions included.

sphcode_ud: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ud \
		FLAGS="-DRADIATING -DDIFFUSING -DGRAVITY -DCOMPVISC"

sphcode_udx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udx \
		FLAGS="-DRADIATING -DDIFFUSING -DEXTGRAV -DCOMPVISC"

sphcode_udz: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udz \
		FLAGS="-DRADIATING -DDIFFUSING -DNOACCEL"

########################################################################
# sphcode_uc: energy is thermodynamic variable, conduction included.

sphcode_uc: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uc \
		FLAGS="-DCONDUCTING -DGRAVITY -DCOMPVISC"

sphcode_ucx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ucx \
		FLAGS="-DCONDUCTING -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_us: energy is thermodynamic variable, starformation included.

sphcode_us: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=us \
		FLAGS="-DSTARFORM -DGRAVITY -DCOMPVISC"

sphcode_usx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=usx \
		FLAGS="-DSTARFORM -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_uis: energy is thermodynamic parameter, starformation included.

sphcode_uis: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uis \
		FLAGS="-DISOTHERMAL -DSTARFORM -DGRAVITY -DCOMPVISC"

sphcode_uisx: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uisx \
		FLAGS="-DISOTHERMAL -DSTARFORM -DEXTGRAV -DCOMPVISC"

########################################################################
# sphcode_uir: energy is thermodynamic parameter, recycling included.

sphcode_uir: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uir \
		FLAGS="-DISOTHERMAL -DSTARFORM -DMASSLOSS -DGRAVITY -DCOMPVISC"

########################################################################
# Build sphcode with self-gravity.

code_con_grav: sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	       sphload_$V.o sphgrav_$V.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V \
	   sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	   sphload_$V.o sphgrav_$V.o -lNBody -lClib -lgsl -lgslcblas -lm

########################################################################
# Build sphcode with external gravity.

code_sin_grav: sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V \
	   sphcode_$V.o sphio_$V.o kdtree_$V.o smooth_$V.o \
	   -lNBody -lGSP -lClib -lgsl -lgslcblas -lm

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

sphdensity: sphdensity.o kdtree_ss.o smooth_ss.o
	$(ZCC) $(ZLDFLAGS) -o sphdensity \
	   sphdensity.o kdtree_ss.o smooth_ss.o \
	   -lNBody -lGSP -lClib -lgsl -lgslcblas -lm

sphdensity.o: sphdensity.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphdensity.o -c sphdensity.c

kdtree_ss.o: kdtree.c kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o kdtree_ss.o -c kdtree.c

smooth_ss.o: smooth.c smooth.h kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FSF) $(FLAGS) \
	    -o smooth_ss.o -c smooth.c

# Generate source code listing.

zeno_sphcode.ps: kdtree.c kdtree.h smooth.c smooth.h sphcode.c sphcode.h \
	         sphdefs.h sphdensity.c sphgrav.c sphio.c  sphload.c
	enscript -M Letterdj -o zeno_sphcode.ps \
	         sphcode.c sphcode.h sphdefs.h kdtree.c kdtree.h \
	         smooth.c smooth.h sphgrav.c sphio.c sphload.c sphdensity.c
