# ~/nbody/graphics/Makefile: mantain N-body graphics utilities.

BINFILES = snapppm snapsmooth snapview viewfile

# Location of X11 libraries.

X11LIB = /opt/X11/lib
# X11LIB = /usr/X11/lib

# Targets used to build Zeno utilities.

snapppm: snapppm.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapppm snapppm.c \
	    -lNBody -lClib -lm

snapsmooth: snapsmooth.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsmooth snapsmooth.c \
	    -lNBody -lClib -lm

snapview: snapview.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -L$(X11LIB) \
	    -o snapview snapview.c -lNBody -lClib \
	    -lglut -lGLU -lGL -lXmu -lX11 -lm

viewfile: viewfile.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o viewfile viewfile.c \
	    -lClib

# Obsolete programs.

snapinterpolate: snapinterpolate.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapinterpolate snapinterpolate.c \
	    -lNBody -lClib -lm
