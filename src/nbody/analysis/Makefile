# ~/nbody/analysis/Makefile: mantain N-body analysis utilities.

BINFILES = snapavg snaphist snapobject snaprholist snapshape snapstat

# Targets used to build Zeno utilities.

snapavg: snapavg.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapavg snapavg.c \
	   -lNBody -lClib -lm

snaphist: snaphist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaphist snaphist.c \
	   -lNBody -lClib -lm

snapobject: snapobject.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapobject snapobject.c \
	   -lNBody -lClib -lm

snaprholist: snaprholist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaprholist snaprholist.c \
	   -lNBody -lClib -lm

snapshape: snapshape.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapshape snapshape.c \
	   -lNBody -lClib -lgsl -lgslcblas -lm

snapstat: snapstat.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstat snapstat.c \
	   -lNBody -lClib -lm

# Development code for dynamic-loading.

snapstat_new: snapstat_new.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstat_new snapstat_new.c \
	   -lNBody -lClib -lm -ldl

sm_x.so: sm_x.c
	$(ZCC) $(ZCCFLAGS) -shared -o sm_x.so -fPIC sm_x.c

# Generate listing of source code.

zeno_analysis.ps: snapavg.c snaphist.c snapobject.c snaprholist.c \
	          snapshape.c snapstat.c snapstat_new.c
	enscript -M Letterdj -o zeno_analysis.ps \
	          snapavg.c snaphist.c snapobject.c snaprholist.c \
	          snapshape.c snapstat.c snapstat_new.c
