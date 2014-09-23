# ~/zeno/src/nbody/sph/tests/Makefile: test SPH codes.

############################################################################################
# Composite hot gas + dark matter model.

i001.dat:
	halogsp out=mod01d.dat m_a=0.32116400 a=1.0 b=10.0 \
	        taper=sw rrange=1/2048:32
	halogsp out=mod01g.dat m_a=0.10705466 a=1.0 b=10.0 \
	        taper=sw rrange=1/2048:32
	gspadd in=mod01g.dat,mod01d.dat out=mod01.dat
	gspsphere gsp=mod01g.dat out=- grav=mod01.dat \
	          nbody=16384 zerocm=false | \
	  snapset in=- out=i001g.dat produce=BodyType type=0x60
	gsprealize gsp=mod01d.dat out=- grav=mod01.dat nbody=49152 | \
	  snapset in=- out=i001d.dat \
	          produce=Uinternal,EntropyFunc,BodyType type=0x40
	cat i001g.dat i001d.dat | \
	  snapcons in=- out=i001.dat nbody=65536 \
	           produce=Position,Velocity,Mass,Uinternal,EntropyFunc,BodyType

############################################################################################
# Runs 100-108: run test cases with default parameters (which may not
# work right in all cases...).

run100:
	sphcode_e out=r100.dat save=r100_state%1d.dat log=r100.log \
	          outputs=EntropyFunc,Position,Velocity,Density,UdotInternal

run101:
	sphcode_ea out=r101.dat save=r101_state%1d.dat log=r101.log \
	           outputs=Position,Velocity,Density,UdotInternal

run102:
	sphcode_u out=r102.dat save=r102_state%1d.dat log=r102.log \
	          outputs=Uinternal,Position,Velocity,Density,UdotInternal

run103:
	sphcode_ui out=r103.dat save=r103_state%1d.dat log=r103.log \
	           outputs=Position,Velocity,Density,UdotInternal

run104:
	sphcode_uc out=r104.dat save=r104_state%1d.dat log=r104.log \
	           outputs=Uinternal,Position,Velocity,Density,UdotInternal

run105:
	sphcode_ur out=r105.dat save=r105_state%1d.dat log=r105.log \
	           outputs=Uinternal,Position,Velocity,Density,UdotInternal

run106:
	sphcode_uo out=r106.dat save=r106_state%1d.dat log=r106.log \
	           outputs=Uinternal,Position,Velocity,Density,UdotInternal

run107:
	sphcode_us out=r107.dat save=r107_state%1d.dat log=r107.log \
	           outputs=Uinternal,Position,Velocity,Density,UdotInternal

############################################################################################
# Runs 110-118: run test cases with adjusted parameters, saving images at each timestep.

NTEST = 16384

run110:
	mkdir -p r110$X
	rm -f r110$X/*
	sphcode_e out=r110$X/r110.dat save=r110$X/r110_state%1d.dat log=r110$X/r110.log \
	          outputs=EntropyFunc,Position,Velocity,Density,UdotInternal \
	          testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="0.2 * udot / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(entf)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r110$X/r110_%03x%s.pgm value=rgb threedim=false

run111:
	mkdir -p r111$X
	rm -f r111$X/*
	sphcode_ea out=r111$X/r111.dat save=r111$X/r111_state%1d.dat log=r111$X/r111.log \
	           outputs=Position,Velocity,Density,UdotInternal \
	           testuint=0.2 testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="0.2 * udot / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(entf)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r111$X/r111_%03x%s.pgm value=rgb threedim=false

run112:
	mkdir -p r112$X
	rm -f r112$X/*
	sphcode_u out=r112$X/r112.dat save=r112$X/r112_state%1d.dat log=r112$X/r112.log \
	          outputs=Uinternal,Position,Velocity,Density,UdotInternal,UdotViscosity \
	          testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="0.2 * udotvis / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r112$X/r112_%03x%s.pgm value=rgb threedim=false

run113:
	mkdir -p r113$X
	rm -f r113$X/*
	sphcode_ui out=r113$X/r113.dat save=r113$X/r113_state%1d.dat log=r113$X/r113.log \
	           outputs=Position,Velocity,Density,UdotInternal,UdotViscosity \
	           testuint=0.4 testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/8" \
	          auxvx="MAX(0, 0.005 * (udot + 20)) / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r113$X/r113_%03x%s.pgm value=rgb threedim=false

run114:
	mkdir -p r114$X
	rm -f r114$X/*
	sphcode_uc out=r114$X/r114.dat save=r114$X/r114_state%1d.dat log=r114$X/r114.log \
	           outputs=Uinternal,Position,Velocity,Density,UdotInternal,UdotViscosity \
	           conduct=0.4 testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="0.2 * udotvis / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r114$X/r114_%03x%s.pgm value=rgb threedim=false

run115:
	mkdir -p r115$X
	rm -f r115$X/*
	sphcode_ur out=r115$X/r115.dat save=r115$X/r115_state%1d.dat log=r115$X/r115.log \
	           outputs="Uinternal,Position,Velocity,Density,UdotInternal, \
	                    UdotViscosity,UdotRadiation" \
	           uradpk=10.0 lambpk=0.1 testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="-0.2 * udotrad / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r115$X/r115_%03x%s.pgm value=rgb threedim=false

run116:
	mkdir -p r116$X
	rm -f r116$X/*
	sphcode_uo out=r116$X/r116.dat save=r116$X/r116_state%1d.dat log=r116$X/r116.log \
	           outputs="Uinternal,Position,Velocity,Density,UdotInternal, \
	                    UdotViscosity,UdotRadiation,OpticalDepth" \
	           uradpk=10.0 lambpk=0.1 opacity=0.25 testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" tau=0 \
	          auxvx="-1.0 * udotrad / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r116$X/r116_%03x%s.pgm value=rgb threedim=false

run117:
	mkdir -p r117$X
	rm -f r117$X/*
	sphcode_ud out=r117$X/r117.dat save=r117$X/r117_state%1d.dat log=r117$X/r117.log \
	           outputs="Uinternal,Position,Velocity,Density,UdotInternal, \
	                    UdotViscosity,UdotRadiation,OpticalDepth" \
	           uradpk=10.0 lambpk=0.1 sigmastar=1.0 opacity=0.25 \
	           testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" tau=0 \
	          auxvx="-1.0 * udotrad / n" \
	          auxvy="rpow(MAX(0, 0.75 + rlog10(uint)/3), 4.0) / n" \
	          auxvz="0.25 / n" | \
	    snapsmooth - pgmout=r117$X/r117_%03x%s.pgm value=rgb threedim=false

run118:
	mkdir -p r118$X
	rm -f r118$X/*
	sphcode_us out=r118$X/r118.dat save=r118$X/r118_state%1d.dat log=r118$X/r118.log \
	           outputs="BodyType,Uinternal,Position,Velocity,Density, \
	                    UdotInternal,UdotViscosity" \
	           starprob=0.25 starlog=r118$X/starlog118.txt \
	           testbody=$(NTEST) dtout=1/32 trace=- | \
	  snapset - - produce=SmoothLength,AuxVec smooth=0.01 \
	          x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" \
	          auxvx="type==0x60 ? 0.2 * udotvis / n : 1.0/n" \
	          auxvy="type==0x60 ? 0.0 : 1.0/n" \
	          auxvz="type==0x60 ? 0.25 / n : 0.0" | \
	    snapsmooth - pgmout=r118$X/r118_%03x%s.pgm value=rgb threedim=false

##############################################################################################

list_diagnostics_11X:
	$(MAKE) -f Makefile ID=110 list_energy list_tstep list_rad list_rho \
	                          list_udot list_entf
	$(MAKE) -f Makefile ID=111 list_energy list_tstep list_rad list_rho \
	                          list_udot
	$(MAKE) -f Makefile ID=112 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint
	$(MAKE) -f Makefile ID=113 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis
	$(MAKE) -f Makefile ID=114 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint
	$(MAKE) -f Makefile ID=115 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint list_udotrad
	$(MAKE) -f Makefile ID=116 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint list_udotrad list_tau
	$(MAKE) -f Makefile ID=117 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint list_udotrad list_tau
	$(MAKE) -f Makefile ID=118 list_energy list_tstep list_rad list_rho \
	                          list_udot list_udotvis list_uint

list_energy: r$(ID)/r$(ID).log
	gawk \
	  '{ \
	     if (out == 1) { \
	       print $$0; \
	       out = 0; \
	     }  \
	     if ($$1 == "Time") { \
	       out = 1; \
	       if (head == 0) { \
	         printf("#%s\n", substr($$0, 2)); \
	         head = 1; \
	       } \
	     } \
	   }' r$(ID)/r$(ID).log > r$(ID)/r$(ID)_energy.txt

list_tstep: r$(ID)/r$(ID).log
	gawk \
	  '{ \
	     if ($$1 == "time:") \
	       tnow = $$2; \
	     if (out == 1) { \
	       printf("%9.6f %9.4f %9.4f %9.4f %9.2f " \
	              "%5.1f %5.1f %5.1f %5.1f %5.1f\n", \
	              tnow, $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9); \
	       out = 0; \
	     }  \
	     if ($$1 == "log" && $$2 == "hmin") { \
	       out = 1; \
	       if (head == 0) { \
	         printf("#%8s %s\n", "Time", substr($$0, 1, 69)); \
	         head = 1; \
	       } \
	     } \
	   }' r$(ID)/r$(ID).log > r$(ID)/r$(ID)_tstep.txt

list_rad: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=r options=time,OCT > r$(ID)/r$(ID)_rad.txt

list_rho: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=rho options=time,OCT > r$(ID)/r$(ID)_rho.txt

list_entf: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=entf options=time,OCT > r$(ID)/r$(ID)_entf.txt

list_uint: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=uint options=time,OCT > r$(ID)/r$(ID)_uint.txt

list_udot: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=udot options=time,OCT > r$(ID)/r$(ID)_udot.txt

list_udotvis: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=udotvis options=time,OCT > r$(ID)/r$(ID)_udotvis.txt

list_udotrad: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=udotrad options=time,OCT > r$(ID)/r$(ID)_udotrad.txt

list_tau: r$(ID)/r$(ID).dat
	snapstat r$(ID)/r$(ID).dat value=tau options=time,OCT > r$(ID)/r$(ID)_tau.txt

##############################################################################################

CRF_LAST = 0x200

combine_rgb_frames:
	gawk \
	  'BEGIN { \
	     for (i = 0; i <= $(CRF_LAST); i++) \
	       printf("rgb3toppm r$(ID)/r$(ID)_%03xr.pgm r$(ID)/r$(ID)_%03xg.pgm " \
	              "r$(ID)/r$(ID)_%03xb.pgm > r$(ID)/r$(ID)_%03x.ppm\n",i,i,i,i); \
	   }' | sh	     

render_frames: r$(ID)/r$(ID).dat
	snapset r$(ID)/r$(ID).dat - produce=SmoothLength,Aux smooth=0.01 \
	        aux="1.0/n" x="0.625 * (rlog10(r) + 1.2)" y="0.875 * vr/4" | \
	  snapsmooth - pgmout=r$(ID)/r$(ID)_%03x.pgm value=bright threedim=false
