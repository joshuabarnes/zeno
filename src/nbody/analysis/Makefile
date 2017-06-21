##  ~/nbody/analysis/Makefile: mantain N-body analysis utilities.
##  Copyright (c) 2017, Joshua E. Barnes, Honolulu, Hawaii.
##  _____________________________________________________________

BINFILES = snapavg$X snaphist$X snapobject$X snaprholist$X \
	   snapshape$X snapstat$X

##  Compilation options are selected by changing OPTIONS (and X):
##	make binfiles OPTIONS="-DDOUBLEPREC" X=_dp
##  compiles BINFILES in double precision, with "_dp" suffixes.
##  _____________________________________________________________

OPTIONS =

binfiles:
	$(MAKE) $(BINFILES)

# Targets used to build Zeno utilities.

snapavg$X: snapavg.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapavg$X \
	   snapavg.c -lNBody -lClib -lm

snaphist$X: snaphist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaphist$X \
	   snaphist.c -lNBody -lClib -lm

snapobject$X: snapobject.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapobject$X \
	   snapobject.c -lNBody -lClib -lm

snaprholist$X: snaprholist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaprholist$X \
	   snaprholist.c -lNBody -lClib -lm

snapshape$X: snapshape.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapshape$X \
	   snapshape.c -lNBody -lClib -lgsl -lgslcblas -lm

snapstat$X: snapstat.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapstat$X \
	   snapstat.c -lNBody -lClib -lm

# Development code for dynamic-loading.

snapstat_new: snapstat_new.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapstat_new \
	   snapstat_new.c -lNBody -lClib -lm -ldl

sm_x.so: sm_x.c
	$(ZCC) $(ZCCFLAGS) -shared $(OPTIONS) -o sm_x.so -fPIC sm_x.c

# Generate listing of source code.

zeno_analysis.ps: snapavg.c snaphist.c snapobject.c snaprholist.c \
	          snapshape.c snapstat.c snapstat_new.c
	enscript -M Letterdj -o zeno_analysis.ps \
	          snapavg.c snaphist.c snapobject.c snaprholist.c \
	          snapshape.c snapstat.c snapstat_new.c
