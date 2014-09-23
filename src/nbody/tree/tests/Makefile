# ~/zeno/src/nbody/tree/tests/Makefile: treecode tests, focusing on
# errors related to softening (March 2012).

EPS016 = 0.0625
EPS064 = 0.015625
EPS256 = 0.00390625
EPSZZZ = 0.0

########################################################################

# cusp10.dat: Mass points in cusp with slope=-1

cusp10.dat:
	makesnap out=- shape=cusp nbody=262144 seed=456 | \
	  snapset - cusp10.dat x="x/rsqrt(r)" y="y/rsqrt(r)" z="z/rsqrt(r)"

# cusp20.dat: Mass points in cusp with slope=-2

cusp20.dat:
	makesnap out=cusp20.dat shape=cusp nbody=262144 seed=123

# test10.dat: Test points in cusp with slope=-1

test10.dat:
	makesnap out=- shape=cusp mass=0 nbody=4096 seed=456 | \
	  snapset - - x="x/rsqrt(r)" y="y/rsqrt(r)" z="z/rsqrt(r)" | \
	    snapsort - test10.dat value=r

# test20.dat: Test points in cusp with slope=-2

test20.dat:
	makesnap out=- shape=cusp mass=0 nbody=4096 seed=456 | \
	  snapsort - test20.dat value=r

# testLx.dat: Test points uniformly distributed in log x along x axis

testLx.dat:
	makesnap out=- mass=0 radius=0 speed=0 nbody=25 | \
	  snapset in=- out=testLx.dat \
	          x="rdex(-3.0 * ((real) i) / (n - 1))"

# testLr.dat: Test points uniformly distributed in log r

testLr.dat:
	makesnap out=- shape=shell mass=0 radius=1 speed=0 nbody=4096 | \
	  snapset in=- out=testLr.dat \
	          x="x * rdex(-3.0 * ((real) i) / (n - 1))" \
	          y="y * rdex(-3.0 * ((real) i) / (n - 1))" \
	          z="z * rdex(-3.0 * ((real) i) / (n - 1))"

########################################################################
# Examine gravity evaluation on a cell-by-cell basis.

test_gravLx20_$(EPS)$X.txt: testLx.dat cusp20.dat
	cat testLx.dat cusp20.dat | \
	  snapcons - - nbody=262169 | \
	    treecode_tf - tstop=0 nstatic=-262144 log=/dev/null options=$(OPT) \
	                theta=0.8 usequad=true eps=$(EPS$(EPS)) > test_gravLx20_$(EPS)$X.txt

extract_each_test: test_gravLx20_$(EPS)$X.txt
	gawk '{ print $$1; }' test_gravLx20_$(EPS)$X.txt | uniq | cat -n | \
	  gawk \
	    '{ printf("make extract_test N=%02d R=%s EPS=$(EPS) X=$X\n",$$1,$$2); }' | sh

extract_test: test_gravLx20_$(EPS)$X.txt
	gawk '{ if ($$1 == $R) print $0; }' \
	     test_gravLx20_$(EPS)$X.txt > test_gravLx20_$(EPS)$X_$N.txt

########################################################################

grav2020_ensemble: cusp20.dat test20.dat
	$(MAKE) grav2020_ds EPS=$(EPS)
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8te OPT=theta-eff
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6te OPT=theta-eff
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4te OPT=theta-eff
	$(MAKE) grav2020_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2te OPT=theta-eff

grav2020_ds: cusp20.dat test20.dat
	rm -f grav2020_$(EPS)_ds.dat
	cat test20.dat cusp20.dat | \
	  snapcons - - nbody=260k | \
	    snapforce in=- out=grav2020_$(EPS)_ds.dat nforce=4096 eps=$(EPS$(EPS))

grav2020_tc: cusp20.dat test20.dat
	rm -f grav2020_$(EPS)_$(ACC).dat
	cat test20.dat cusp20.dat | \
	  snapcons - - nbody=260k | \
	    treecode - trace=- tstop=0 log=/dev/null nstatic=-262144 options=$(OPT) \
	             usequad=$(QUAD) theta=$(THETA) eps=$(EPS$(EPS)) | \
	      snapsift - grav2020_$(EPS)_$(ACC).dat sieve="i<4096"

grav10Lr_ensemble: cusp10.dat testLr.dat
	$(MAKE) grav10Lr_ds EPS=$(EPS)
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8te OPT=theta-eff
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6te OPT=theta-eff
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4te OPT=theta-eff
	$(MAKE) grav10Lr_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2te OPT=theta-eff

grav10Lr_ds: cusp10.dat testLr.dat
	rm -f grav10Lr_$(EPS)_ds.dat
	cat testLr.dat cusp10.dat | \
	  snapcons - - nbody=260k | \
	    snapforce in=- out=grav10Lr_$(EPS)_ds.dat nforce=4096 eps=$(EPS$(EPS))

grav10Lr_tc: cusp10.dat testLr.dat
	rm -f grav10Lr_$(EPS)_$(ACC).dat
	cat testLr.dat cusp10.dat | \
	  snapcons - - nbody=260k | \
	    treecode - trace=- tstop=0 log=/dev/null nstatic=-262144 options=$(OPT) \
	             usequad=$(QUAD) theta=$(THETA) eps=$(EPS$(EPS)) | \
	      snapsift - grav10Lr_$(EPS)_$(ACC).dat sieve="i<4096"

grav20Lr_ensemble: cusp20.dat testLr.dat
	$(MAKE) grav20Lr_ds EPS=$(EPS)
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.8 ACC=q8te OPT=theta-eff
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.6 ACC=q6te OPT=theta-eff
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.4 ACC=q4te OPT=theta-eff
	$(MAKE) grav20Lr_tc EPS=$(EPS) QUAD=t THETA=0.2 ACC=q2te OPT=theta-eff

grav20Lr_ds: cusp20.dat testLr.dat
	rm -f grav20Lr_$(EPS)_ds.dat
	cat testLr.dat cusp20.dat | \
	  snapcons - - nbody=260k | \
	    snapforce in=- out=grav20Lr_$(EPS)_ds.dat nforce=4096 eps=$(EPS$(EPS))

grav20Lr_tc: cusp20.dat testLr.dat
	rm -f grav20Lr_$(EPS)_$(ACC).dat
	cat testLr.dat cusp20.dat | \
	  snapcons - - nbody=260k | \
	    treecode - trace=- tstop=0 log=/dev/null nstatic=-262144 options=$(OPT) \
	             usequad=$(QUAD) theta=$(THETA) eps=$(EPS$(EPS)) | \
	      snapsift - grav20Lr_$(EPS)_$(ACC).dat sieve="i<4096"

########################################################################

CON = 2020

compare_accel:
	(cat grav$(CON)_$(EPS)_ds.dat; \
	 snapset grav$(CON)_$(EPS)_$(ACC).dat - passall=f \
	         require=Acceleration produce=AuxVec \
	         auxvx=ax auxvy=ay auxvz=az) | \
	  snapcollect - - produce=Position,Acceleration,AuxVec | \
	    snapset - - produce=Aux \
                    aux=r \
	            auxvx="dotvp(acc,pos)/r" \
	            auxvy="distv(auxv,acc) / absv(acc)" \
	            auxvz="(dotvp(auxv,pos)/r - dotvp(acc,pos)/r) / \
	                    distv(auxv,acc)" | \
	      snaplist - Aux,AuxVec auxhead="r,a_r,da/a,cos(xi)" \
	        > cmpacc$(CON)_$(EPS)_$(ACC).txt

avg_accel_error: cmpacc$(CON)_$(EPS)_$(ACC).txt
	gawk '{ \
	  if (NR == 1) \
	    printf("#%9s  %10s\n", "ravg", "erravg"); \
	  if (NR > 1) { \
	    ravg += log($$1); eavg += log($$3); n++; \
	    if (n == 16) { \
	      printf("%10.4e  %10.4e\n", exp(ravg / 16), exp(eavg / 16)); \
	      ravg = eavg = n = 0; \
	    } \
	  } \
	}' cmpacc$(CON)_$(EPS)_$(ACC).txt > cmpacc$(CON)_$(EPS)_$(ACC)_avg.txt
