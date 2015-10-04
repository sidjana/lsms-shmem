// lsms class to encapsulate a version of LSMS_1.9 for use in gWL etc.

#ifndef LSMS_CLASS_H
#define LSMS_CLASS_H

//#include <mpi.h>

#include <vector>

#include "Real.hpp"
#include "../Potential/PotentialShifter.hpp"
#include "mixing.hpp"


extern "C"{
unsigned long long get_rtc();
}

class LSMS {
public:
 // LSMS(MPI_Comm comm, const char * i_lsms, const char * out_prefix);
  LSMS(SHMEM_activeset comm_shmem, const char * i_lsms, const char * out_prefix);
  ~LSMS();

  int version(){return LSMS_version;}
  int numSpins(){return lsms.num_atoms;}
  void setEvec(std::vector<std::vector<Real> > &);
  void setEvec(Real *);
  void setEvecAndSpinPotentialShift(Real *);
  void getEvec(std::vector<std::vector<Real> > &);
  void getEvec(Real *);
  void getMag(std::vector<std::vector<Real> > &);
  void getMag(Real *);
  Real oneStepEnergy();
  Real oneStepEnergy(Real *eb);
  Real oneStepEnergy(std::vector<std::vector<Real> > &ev){setEvec(ev); return oneStepEnergy();}
  Real multiStepEnergy();
  Real scfEnergy(Real *eb);
  Real scfEnergy() {Real eb; return scfEnergy(&eb);}
  void setEfTol(Real e) {efTol=e;}
  Real getEf(void) {return lsms.chempot;}
  void writePot(char *name);

  Real energyDifference;
private:
  char prefix[256];
  int LSMS_version;
  static int max_num_local;
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  static LocalTypeInfo local;
  MixingParameters mix;

  PotentialShifter potentialShifter;

  Mixing *mixing;

  Real efTol;
  Real energyTolerance;
};


#endif

