OBJ = zeroout.o fstop.o newint.o fit.o interp.o trltog.o spin_trafo.o fitpot.o rwave.o \
      initwave.o dfv.o dfv_new.o bsstep.o mmid.o rzextr.o zsphbesjh.o zsphbesj.o cinterp.o cgaunt_c.o \
      gaunt.o gauleg.o plglmax.o ifacts_c.o clock_time.o cmtruni.o \
      constraint.o v_plus_minus.o congauss_c.o u_sigma_u.o fnpi.o wrtmtx.o mbeqa.o \
      cnewint.o getclm.o ylag.o ricbes.o newder.o \
      bulirsch_stoer.o mod_midpoint.o \
      calculateGauntCoeficients.o \
      readLastLine.o #rtc.o

all: libMisc.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

#rtc.o: rtc.c
#	$(CC) -c rtc.c

libMisc.a : $(OBJ)
	$(AR) -r libMisc.a $(OBJ) 
	cp libMisc.a $(TOP_DIR)/lib
