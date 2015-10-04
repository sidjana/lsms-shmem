OBJ = calculateTauMatrix.o makegij_c.o setgij.o block_inverse.o zblock_lu.o wasinv.o zmar1.o wasinv_p.o \
      zuqmx.o zutfx.o zucpx.o zaxpby.o zrandn.o tau_inv_postproc.o trgtol.o green_function.o gf_local.o \
      int_zz_zj.o mdosms_c.o mgreen_c.o write_kkrmat.o

all: libMultipleScattering.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.F
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libMultipleScattering.a : $(OBJ)
	$(AR) -r libMultipleScattering.a $(OBJ) 
	cp libMultipleScattering.a $(TOP_DIR)/lib
