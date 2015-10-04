OBJ = distributeAtoms.o LSMSCommunication.o 

all: libCommunication.a

clean:
	rm -f *.a *.o

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<  #../../instr.o ../../rtc.o  /home/siddhart/opt/x86_adapt/build/libx86_adapt.so


libCommunication.a : $(OBJ)
	$(AR) -r libCommunication.a $(OBJ)  #../../instr.o ../../rtc.o   /home/siddhart/opt/x86_adapt/build/libx86_adapt.so
	cp libCommunication.a $(TOP_DIR)/lib
