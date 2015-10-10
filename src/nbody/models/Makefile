# ~/nbody/models/Makefile: build programs for model construction.

BINFILES = kepsnap plumsnap polysnap snaptestdisk

# Targets used to build Zeno utilities.

kepsnap: kepsnap.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o kepsnap kepsnap.c \
	      -lNBody -lClib -lm

plumsnap: plumsnap.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o plumsnap plumsnap.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

polysnap: polysnap.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o polysnap polysnap.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

snaptestdisk: snaptestdisk.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaptestdisk snaptestdisk.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

# Under development

extforce: extforce.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED -o extforce extforce.c \
	      -lNBody -lClib -lm

exptest: exptest.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED -o exptest exptest.c \
	      -lNBody -lClib -lm

# Generate listing of source code.

zeno_models.ps: kepsnap.c plumsnap.c polysnap.c snaptestdisk.c
	enscript -M Letterdj -o zeno_models.ps \
	          kepsnap.c plumsnap.c polysnap.c snaptestdisk.c
