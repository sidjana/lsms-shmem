OBJ = SingleSiteScattering.o readSingleAtomData_hdf5.o readSingleAtomData_bigcell.o \
      F_readSingleAtomData_bigcell.o single_scatterer_nonrel.o semrel_c.o scalar_m.o \
      writeSingleAtomData_hdf5.o writeSingleAtomData_bigcell.o \
      F_writeSingleAtomData_bigcell.o 

all: libSingleSite.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.f90
	$(F77) -c -o $@ $<

libSingleSite.a : $(OBJ)
	$(AR) -r libSingleSite.a $(OBJ) 
	cp libSingleSite.a $(TOP_DIR)/lib
