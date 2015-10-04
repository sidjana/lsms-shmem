#ifndef LSMS_SYSTEM_PARAM_H
#define LSMS_SYSTEM_PARAM_H
#include <stdio.h>
#include <string.h>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "SingleSite/AtomData.hpp"

#include "Misc/Indices.hpp"

class LSMSGlobals {
public:
  void setIstop(const char *c){strncpy(istop,c,32); for(int i=strlen(c); i<32; i++) istop[i]=' ';}
  bool checkIstop(const char *c){return (strncmp(istop,c,32)==0);}
  int iprpts,ipcore;
  int iprint;
  int print_node,default_iprint;
  char istop[32];

// for GPU only
  int GPUThreads;
};

class EnergyContourParameters {
public:
  int grid, npts;
  Real ebot,etop,eitop,eibot;
// Grouping of energies for single site solver
  int maxGroupSize;
  int groupSize() {int nume=npts; if(grid==0) nume=1; else if(grid==2) nume++; return std::min(nume,maxGroupSize);}
};

class LSMSSystemParameters {
public:
  char systemid[80];
  char title[80];
  char potential_file_in[128];
  char potential_file_out[128];
  int pot_in_type,pot_out_type;
  int mixing; // combines LSMS_1's mix_quant & mix_algor : -1 don't mix. mix_quant=mixing%4; mix_algor=mixing>>2;
              // mix_quant  0: charge, 1: potential
              // mix_algor  0: simple (linear) mixing; 1: broyden
              // charge, simple: 0; broyden: 4
              // potential, simple: 1; boyden: 5
  Real alphaDV; // mixing parameter for density or potential
  int num_atoms;
  int nspin;
  int nrel_rel;
  int nrelc,nrelv;
  int n_spin_cant;
  int n_spin_pola;
  int nscf;
  int writeSteps;
  int mtasa;
  int vSpinShiftFlag; // if !=0 : shift the spin up and down potentials according to atom.vSpinShift
                      // this is used in WL-LSMS with moment magnitude fluctuations
  int fixRMT; // n_fix_mt from LSMS_1:
              //   0 -> set rmt to calculated inscribed sphere in atomic volume
              //   1 -> set it to rmt read in from potential file
  Real clight;
  int maxlmax;
  LSMSGlobals global;
  AngularMomentumIndices angularMomentumIndices;
  EnergyContourParameters energyContour;

// no. of Gaussian points for volume integration
  int ngaussr,ngaussq;

// Properties of the whole system:
  Real chempot;                // Chemical potential
  Real zvaltss;                // Total valence charge
  Real volumeTotal;            // Total cell volume
  Real volumeNorm;             // Volume renormalization factor
  Real volumeInterstitial;     // Total interstitial volume
  Real u0;                     // Contribution of the Muffin-tin zero potential to the Coulomb energy
  Real totalEnergy;            // Total energy
  //Real pressure;               // Pressure

};

extern const char *potentialTypeName[];

class AtomType {
public:
  AtomType() : pot_in_idx(-1), store_id(-1) {}
  char name[4];
  int lmax,Z,Zc,Zs,Zv;
  int first_instance, number_of_instances;
  Real rsteps[4];
  Real rLIZ, rad;
  int node,local_id;
  int store_id;   // position in tmatStore
  int pot_in_idx;
};

class CrystalParameters {
public:
  int maxlmax;
  CrystalParameters() : bravais(3,3) {}
  void resize(size_t n) {type.resize(n); position.resize(3,n); evecs.resize(3,n);}
  void resizeTypes(size_t n) {types.resize(n);}
  Matrix<Real> bravais;
  Real omega; // bravais lattice volume
  int num_atoms,num_types;
  Matrix<Real> position,evecs;
  std::vector<int> type;
  std::vector<AtomType> types;
};



class LocalTypeInfo {
public:
  //LocalTypeInfo() : num_local(0), atom(0) {}
  //~LocalTypeInfo() { //if(num_local) delete[] atom;}
  // void setNumLocal(int n) {num_local=n; atom = new AtomData[n]; global_id.resize(n);}
  void setNumLocal(int n)
  {
    num_local=n; atom.resize(n); global_id.resize(n); n_per_type.resize(n);
    for(int i=0; i<num_local; i++) atom[i].reset();
  }
  void setGlobalId(int rank,CrystalParameters &crystal)
  {
    for(int i=0; i<crystal.num_types; i++)
    {
      if(rank==crystal.types[i].node)
      {
        global_id[crystal.types[i].local_id]=i;
        n_per_type[crystal.types[i].local_id]=crystal.types[i].number_of_instances;
      }
    }
  }
  void setMaxPts(int n) {for(int i=0; i<num_local; i++) atom[i].resizePotential(n);}
  void setMaxCore(int n) {for(int i=0; i<num_local; i++) atom[i].resizeCore(n);}
  int maxNrmat(void) {int v=0; for(int i=0; i<num_local; i++) if(atom[i].nrmat>v) v=atom[i].nrmat; return v;}
  int num_local;
  std::vector<int> global_id;
  std::vector<AtomData> atom;
  std::vector<int> n_per_type;

  int lDimTmatStore,blkSizeTmatStore;
  Matrix<Complex> tmatStore;
  std::vector<int> tmatStoreGlobalIdx;

  Real qrms[2];
};

void printLSMSGlobals(FILE *f,LSMSSystemParameters &lsms);
void printLSMSSystemParameters(FILE *f,LSMSSystemParameters &lsms);
void printCrystalParameters(FILE *f, CrystalParameters &crystal);
void printLocalTypeInfo(FILE *f, LocalTypeInfo &local);
void printLIZInfo(FILE * f, AtomData &atom);

#endif
