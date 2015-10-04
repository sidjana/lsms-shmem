// lsms class to encapsulate a version of LSMS_1.9 for use in gWL etc.

//#include <mpi.h>
#include <iostream>
#include <vector>

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Madelung/Madelung.hpp"
#include "VORPOL/VORPOL.hpp"
#include "Potential/calculateChargesPotential.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "calculateDensities.hpp"
#include "calculateEvec.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"

#include "lsmsClass.hpp"

/* redefining static member variable */ 
int LSMS::max_num_local;

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
#endif

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

#ifdef BUILDKKRMATRIX_GPU
void *allocateDStore(void);
void freeDStore(void *);
void *allocateDConst(void);
void freeDConst(void *);

std::vector<void *> deviceConstants;
// std::vector<void *> deviceStorage;
void * deviceStorage;
#endif

void initLSMSLuaInterface(lua_State *L);
int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix);
void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local);
void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local,
                 SphericalHarmonicsCoeficients &shc);
void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);


// oneStepEnergy calculates the frozen potential energy without converging the Fermi energy
Real LSMS::oneStepEnergy(Real *eb)
{
  static Real eband;

  calculateCoreStates(comm,lsms,local);
  energyContourIntegration(comm,lsms,local);
  calculateChemPot(comm,lsms,local,eband);
  *eb=eband;
  return eband;
}

Real LSMS::oneStepEnergy()
{
  Real eband;
  return oneStepEnergy(&eband);
}

Real LSMS::multiStepEnergy()
{
  static Real eband,ef;
  int iterationCount=1;

  ef=lsms.chempot;

  if(lsms.vSpinShiftFlag!=0)
  {
    potentialShifter.applyShifts(local);
  }

  calculateCoreStates(comm,lsms,local);

  energyContourIntegration(comm,lsms,local);
  calculateChemPot(comm,lsms,local,eband);
  while(std::abs(ef-lsms.chempot)>efTol && iterationCount<lsms.nscf)
  {
    ef=lsms.chempot;
    iterationCount++;
    energyContourIntegration(comm,lsms,local);
    calculateChemPot(comm,lsms,local,eband);
  }

  // calculate the Zeeman contribution from the spin shift and adjust the band energy acordingly
  if(1==0)
  if(lsms.vSpinShiftFlag!=0)
  {
    calculateEvec(lsms,local);
    calculateAllLocalChargeDensities(lsms,local);
    calculateLocalCharges(lsms,local);
    static Real eZeeman=0.0;
    for(int i=0; i<local.num_local; i++)
    {
      Real s=1.0;
      printf(" local.atom[%d].qvalws = %lf\n",i,local.atom[i].qvalws);
      printf(" local.atom[%d].mvalws = %lf\n",i,local.atom[i].mvalws);
      printf(" local.atom[%d].xvalws = %lf  %lf\n",i,local.atom[i].xvalwsNew[0],local.atom[i].xvalwsNew[1]);
      printf(" local.atom[%d].evec = %lf  %lf  %lf\n",i,local.atom[i].evec[0],local.atom[i].evec[1],local.atom[i].evec[2]);
      printf(" local.atom[%d].evecNew = %lf  %lf  %lf\n",i,local.atom[i].evecNew[0],local.atom[i].evecNew[1],local.atom[i].evecNew[2]);
      printf(" local.atom[%d].afm = %d\n", i,local.atom[i].afm);
/*
      if(local.atom[i].evec[0]*local.atom[i].evecNew[0]+
         local.atom[i].evec[1]*local.atom[i].evecNew[1]+
         local.atom[i].evec[2]*local.atom[i].evecNew[2] < 0.0) s=-1.0;
      eZeeman+=-s*2.0*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                   +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
*/
      eZeeman+=-s*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                   +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
    }
    globalSum_double(comm,eZeeman);
    // if(lsms.global.iprint>=0)
    {

      printf("LSMS::multiStepEnergy(): eZeeman=%lf\n",eZeeman);
    }
    eband-=eZeeman;
  }

  if(lsms.global.iprint>=0)
  {
    if(iterationCount<lsms.nscf) printf("LSMS::multiStepEnergy() converged in %d steps.\n",iterationCount);
    else printf("LSMS::multiStepEnergy() did not converge in %d steps.\n",lsms.nscf);
  }
  return eband;
}

Real LSMS::scfEnergy(Real *eb)
{
  FILE *kFile=NULL;
  Real oldTotalEnergy=0.0;
  
  static Real eband;



  int iterationCount=0;
  if(lsms.global.iprint >= 0) printf("Total number of iterations:%d\n",lsms.nscf);

  // double timeScfLoop=get_rtc();
  // double timeCalcChemPot = 0.0;
  if(lsms.vSpinShiftFlag!=0)
  {
    potentialShifter.applyShifts(local);
  }

  calculateCoreStates(comm,lsms,local);

  // FILE *kFile=fopen("k.out","w");

  for(iterationCount=0; iterationCount<lsms.nscf; iterationCount++)
  {
    //if(lsms.global.iprint>=0)
    if(comm.comm.rank==0)
      printf("Iteration %d:\n",iterationCount);

    energyContourIntegration(comm,lsms,local);
    // double dTimeCCP = get_rtc();
    // if(!lsms.global.checkIstop("buildKKRMatrix"))
    calculateChemPot(comm,lsms,local,eband);
    // dTimeCCP=get_rtc() - dTimeCCP;
    // timeCalcChemPot += dTimeCCP;
    calculateEvec(lsms,local);
    calculateAllLocalChargeDensities(lsms,local);
    calculateChargesPotential(comm,lsms,local,crystal,0);
    calculateTotalEnergy(comm,lsms,local,crystal);

    mixing -> updateChargeDensity(lsms,local.atom);

    // If charge is mixed, recalculate the potential  (need a flag for this from input)
    calculateChargesPotential(comm,lsms,local,crystal,1);
    mixing -> updatePotential(lsms,local.atom);

    energyDifference=oldTotalEnergy-lsms.totalEnergy;
    if(comm.comm.rank==0)
    {
      printf("Band Energy = %lf Ry %10s", eband, "");
      printf("Fermi Energy = %lf Ry\n", lsms.chempot);
      printf("Total Energy = %lf Ry\n", lsms.totalEnergy);
      printf("Energy Change = %lg Ry\n", energyDifference);
    }

    if(kFile!=NULL)
      fprintf(kFile,"%3d %20.12lf %12.6lf\n",iterationCount,lsms.totalEnergy,lsms.chempot);

    // calculate core states for new potential if we are performing scf calculations
    calculateCoreStates(comm,lsms,local);

// check for convergence
    if(std::abs(energyDifference)<energyTolerance)
      break;

    oldTotalEnergy=lsms.totalEnergy;
  }

  //calculateChemPot(comm,lsms,local,*eband);
  *eb=eband;

  return lsms.totalEnergy;
}

// still need LSMS::scfEnergy

void LSMS::writePot(char * name)
{
  printf("*******    lsms::writePot not implemented yet! *********\n");
}

//LSMS::LSMS(MPI_Comm _comm, const char * i_lsms, const char * out_prefix)
LSMS::LSMS(SHMEM_activeset _comm, const char * i_lsms, const char * out_prefix)
{
  char h[]="NO OUTPUT";

  if(out_prefix==NULL) strncpy(prefix,h,255);
  else strncpy(prefix,out_prefix,255);

/*
  if(rank==0)
    {
      std::cout<<"initializing lsms::lsms with i_lsms="<<i_lsms<<" and out_prefix="
	       <<prefix<<std::endl;
    }
*/

  efTol=1.0e-6;
  energyTolerance=1.0e-8;

  lua_State *L=lua_open();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  //TODO: check whether _comm will always be MPI_COMM_WORLD or 
  //not
  initializeCommunication(comm, _comm);

  lsms.global.iprpts=1051;
  lsms.global.ipcore=15;
  lsms.global.setIstop("main");
  lsms.global.iprint=-1;
  lsms.global.default_iprint=-1;
  lsms.global.print_node=0;
  lsms.global.GPUThreads=1;
#ifdef _OPENMP
  lsms.global.GPUThreads=std::min(16,omp_get_max_threads());
#else
  lsms.global.GPUThreads=1;
#endif
  lsms.ngaussr=10;
  lsms.ngaussq=40;
  // if(comm.rank==comm.size-1) lsms.global.iprint=0;

  lsms.vSpinShiftFlag=0;

  if(comm.comm.rank==0)
  {
    if(strncmp(prefix,"0_",2)==0)
    {
      printf("LSMS_3: Program started\n");
      //printf("Using %d MPI processes\n",comm.size);
      printf("Using %d OpenSHMEM PEs\n",comm.comm.size);
      printf("Using %d OpenMP threads\n",omp_get_max_threads());
      acceleratorPrint();
      printf("Reading input file '%s'\n",i_lsms);
    }

    if(luaL_loadfile(L, i_lsms) || lua_pcall(L,0,0,0))
    {
      fprintf(stderr,"!! Cannot run input file '%s'!!\n",i_lsms);
      exit(1);
    }

    if(readInput(L,lsms,crystal,mix))
    {
      fprintf(stderr,"!! Something wrong in input file!!\n");
      exit(1);
    }
  }

  communicateParameters(comm,lsms,crystal,mix);
  if(comm.comm.rank!=lsms.global.print_node) lsms.global.iprint=lsms.global.default_iprint;

  local.setNumLocal(distributeTypes(crystal, comm));
  max_num_local=local.num_local;
  globalMax_int(comm,max_num_local);
  local.setGlobalId(comm.comm.rank,crystal);

  lsms.angularMomentumIndices.init(2*crystal.maxlmax);
  sphericalHarmonicsCoeficients.init(2*crystal.maxlmax);

  gauntCoeficients.init(lsms,lsms.angularMomentumIndices,sphericalHarmonicsCoeficients);
  iFactors.init(lsms,crystal.maxlmax);

  buildLIZandCommLists(comm, lsms, crystal, local);

// initialize the potential accelerators (GPU)
// we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
// which is only available after building the LIZ

  acceleratorInitialize(lsms.n_spin_cant*local.maxNrmat(),lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.resize(local.num_local);
  for(int i=0; i<local.num_local; i++) deviceConstants[i]=allocateDConst();
  // deviceStorage.resize(1);
  deviceStorage=allocateDStore();
#endif

  for(int i=0; i<local.num_local; i++)
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());

// set maximal number of radial grid points and core states if reading from bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if(lsms.global.iprint>=0)
  {
    printLSMSSystemParameters(stdout,lsms);
    printCrystalParameters(stdout,crystal);
  }
  if(lsms.global.iprint>=1)
  {
    fprintf(stdout,"LIZ for atom 0 on this node\n");
    printLIZInfo(stdout,local.atom[0]);
    printCommunicationInfo(stdout, comm);
  }
 
 
  loadPotentials(comm,lsms,crystal,local);

  setupVorpol(lsms,crystal,local,sphericalHarmonicsCoeficients);

  calculateVolumes(comm,lsms,crystal,local);

// need to calculate madelung matrices
  calculateMadelungMatrices(lsms,crystal,local);

  if(lsms.global.iprint>=0)
  {
    printLocalTypeInfo(stdout,local);
  }

// initialize Mixing
  setupMixing(mix, mixing);

  potentialShifter.resize(local.num_local);
  potentialShifter.resetPotentials(local);

  LSMS_version=3000;
//  if(comm.rank==0)
//    std::cout<<prefix<<" lsms::lsms - LSMS version is "<<LSMS_version<<std::endl; 

  lua_close(L);
}

LSMS::~LSMS()
{
  local.tmatStore.unpinMemory();
#ifdef BUILDKKRMATRIX_GPU
  for(int i=0; i<local.num_local; i++) freeDConst(deviceConstants[i]);
  freeDStore(deviceStorage);
#endif
  acceleratorFinalize();
//  finalizeCommunication();
}

void LSMS::setEvec(std::vector<std::vector<Real> > &ev)
{
  //MPI_Status status;
  //distribute the evecs to the respective processes
  if(comm.comm.rank==0)
  {
    //std::vector<MPI_Request> request(crystal.num_types);
    int n_req=0;
    for(int p=0; p<crystal.num_types; p++)
    {
      int l=crystal.types[p].local_id;
      int n=crystal.types[p].node;
      if(n==0)
      {
        local.atom[l].evec[0]=(ev[p])[0]; local.atom[l].evec[1]=(ev[p])[1]; local.atom[l].evec[2]=(ev[p])[2];
      } else {
        //MPI_Isend(&(ev[p])[0],3,MPI_DOUBLE,n,l,comm.comm,&request[n_req++]);
        shmem_double_put(&local.atom[p].evec[0], &(ev[p])[0], 3, comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*n);
      }
    }
    //for(int i=0; i<n_req; i++)
      //MPI_Wait(&request[i],&status);
  } else {
    /*
    std::vector<MPI_Request> request(local.num_local);
    for(int p=0; p<local.num_local; p++)
    {
      MPI_Irecv(&local.atom[p].evec[0],3,MPI_DOUBLE,0,p,comm.comm,&request[p]);
    }
    for(int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i],&status);
    */
  }
  shmem_barrier_all();
}

void LSMS::setEvec(Real *ev)
{
  //MPI_Status status;
  //distribute the evecs to the respective processes
  if(comm.comm.rank==0)
  {
    //std::vector<MPI_Request> request(crystal.num_types);
    int n_req=0;
    for(int p=0; p<crystal.num_types; p++)
    {
      int l=crystal.types[p].local_id;
      int n=crystal.types[p].node;
      if(n==0)
      {
        local.atom[l].evec[0]=ev[3*p+0]; local.atom[l].evec[1]=ev[3*p+1]; local.atom[l].evec[2]=ev[3*p+2];
      } else {
        //MPI_Isend(&ev[3*p],3,MPI_DOUBLE,n,l,comm.comm,&request[n_req++]);
        shmem_double_put(&local.atom[p].evec[0], &ev[3*p], 3, comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*n);
      }
    }
    //for(int i=0; i<n_req; i++)
      //MPI_Wait(&request[i],&status);
  } else {
    /*
    std::vector<MPI_Request> request(local.num_local);
    for(int p=0; p<local.num_local; p++)
    {
      MPI_Irecv(&local.atom[p].evec[0],3,MPI_DOUBLE,0,p,comm.comm,&request[p]);
    }
    for(int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i],&status);
      */
  }
  shmem_barrier_all();

}

void LSMS::setEvecAndSpinPotentialShift(Real *ev)
{
  //MPI_Status status;
  lsms.vSpinShiftFlag=1;
  //distribute the evecs to the respective processes
  if(comm.comm.rank==0)
  {
    //std::vector<MPI_Request> request(crystal.num_types);
    int n_req=0;
    for(int p=0; p<crystal.num_types; p++)
    {
      int l=crystal.types[p].local_id;
      int n=crystal.types[p].node;
      if(n==0)
      {
        local.atom[l].evec[0]=ev[3*p+0]; local.atom[l].evec[1]=ev[3*p+1]; local.atom[l].evec[2]=ev[3*p+2];
      } else {
        //MPI_Isend(&ev[3*p],3,MPI_DOUBLE,n,l,comm.comm,&request[n_req++]);
        shmem_double_put(&local.atom[p].evec[0], &ev[3*p], 3, comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*n);
      }
    }
    //for(int i=0; i<n_req; i++)
      //MPI_Wait(&request[i],&status);
  } else {
    /*
    std::vector<MPI_Request> request(local.num_local);
    for(int p=0; p<local.num_local; p++)
    {
      MPI_Irecv(&local.atom[p].evec[0],3,MPI_DOUBLE,0,p,comm.comm,&request[p]);
    }
    for(int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i],&status);
      */
  }
  shmem_barrier_all();

  for(int p=0; p<local.num_local; p++)
  {
    Real m=std::sqrt(local.atom[p].evec[0]*local.atom[p].evec[0]+
		  local.atom[p].evec[1]*local.atom[p].evec[1]+
		  local.atom[p].evec[2]*local.atom[p].evec[2]);
    local.atom[p].vSpinShift=m-1.0;
    local.atom[p].evec[0]=local.atom[p].evec[0]/m;
    local.atom[p].evec[1]=local.atom[p].evec[1]/m;
    local.atom[p].evec[2]=local.atom[p].evec[2]/m;
  }
}


void LSMS::getEvec(std::vector<std::vector<Real> > &ev)
{
  Array3d<Real> r_buf(4,max_num_local,comm.comm.size);
  //Matrix<Real> s_buf(4,max_num_local);
  Real* s_buf = (Real*)shmalloc(4*max_num_local*sizeof(Real));
  int row_size=max_num_local;

  for(int i=0; i<max_num_local; i++)
  {
    //s_buf(0,i)=-1.0;
    *(s_buf+(i))=-1.0;
  }

  for(int i=0; i<local.num_local; i++)
  {
    //s_buf(0,i)=Real(local.global_id[i]);
    *(s_buf+(i))=Real(local.global_id[i]);
    //s_buf(1,i)=local.atom[i].evec[0]; s_buf(2,i)=local.atom[i].evec[1]; s_buf(3,i)=local.atom[i].evec[2];
    *(s_buf+(1*row_size+i))=local.atom[i].evec[0]; *(s_buf+(2*row_size+i))=local.atom[i].evec[1]; *(s_buf+(3*row_size+i))=local.atom[i].evec[2];
  }

  //MPI_Gather(&s_buf(0,0),4*max_num_local,MPI_DOUBLE,&r_buf(0,0,0),4*max_num_local,MPI_DOUBLE,0,comm.comm);
  //Naive Gather implementation
  if(comm.comm.rank=0)
    for(int i=0;i<comm.comm.size;i++)
       shmem_double_get(&r_buf(0,0,0),s_buf,4*max_num_local,comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*i);

  shmem_barrier_all();

  if(comm.comm.rank==0)
  {
    for(int p=0; p<comm.comm.size; p++)
    {
      for(int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
        {
          int j=r_buf(0,i,p);
          (ev[j])[0]=r_buf(1,i,p); (ev[j])[1]=r_buf(2,i,p); (ev[j])[2]=r_buf(3,i,p);
        }
    }
  }
}

void LSMS::getEvec(Real *ev)
{
  Array3d<Real> r_buf(4,max_num_local,comm.comm.size);
  //Matrix<Real> s_buf(4,max_num_local);
  Real* s_buf = (Real*)shmalloc(4*max_num_local*sizeof(Real));
  int row_size = max_num_local;

  for(int i=0; i<max_num_local; i++)
  {  //s_buf(0,i)=-1.0;
     *(s_buf+(i))=-1.0;
  }

  for(int i=0; i<local.num_local; i++)
  {
    //s_buf(0,i)=Real(local.global_id[i]);
    *(s_buf+(i))=Real(local.global_id[i]);
    //s_buf(1,i)=local.atom[i].evec[0]; s_buf(2,i)=local.atom[i].evec[1]; s_buf(3,i)=local.atom[i].evec[2];
    *(s_buf+(1*row_size+i))=local.atom[i].evec[0]; *(s_buf+(2*row_size+i))=local.atom[i].evec[1]; *(s_buf+(3*row_size+i))=local.atom[i].evec[2];
  }

  //MPI_Gather(&s_buf(0,0),4*max_num_local,MPI_DOUBLE,&r_buf(0,0,0),4*max_num_local,MPI_DOUBLE,0,comm.comm);


  //Naive Gather implementation
  if(comm.comm.rank=0)
     for(int i=0;i<comm.comm.size;i++)
         shmem_double_get(&r_buf(0,0,0),s_buf,4*max_num_local,comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*i);
  
  shmem_barrier_all();
  

  if(comm.comm.rank==0)
  {
    for(int p=0; p<comm.comm.size; p++)
    {
      for(int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
        {
          int j=r_buf(0,i,p);
          ev[3*j+0]=r_buf(1,i,p); ev[3*j+1]=r_buf(2,i,p); ev[3*j+2]=r_buf(3,i,p);
        }
    }
  }

}

void LSMS::getMag(std::vector<std::vector<Real> > &ev)
{
  Array3d<Real> r_buf(4,max_num_local,comm.comm.size);
  //Matrix<Real> s_buf(4,max_num_local);
  double* s_buf = (double*) shmalloc(4*max_num_local*sizeof(Real));
  int row_size=max_num_local;
  Real mag;

  for(int i=0; i<max_num_local; i++)
  {
    //s_buf(0,i)=-1.0;
    *(s_buf+i)=-1.0;

  }

  for(int i=0; i<local.num_local; i++)
  {
    //s_buf(0,i)=Real(local.global_id[i]);
    *(s_buf+(i))=Real(local.global_id[i]);
    //s_buf(0,i)=Real(local.global_id[i]);
    *(s_buf+(i))=Real(local.global_id[i]);

/* The following comment has been left untouched.
    s_buf(1,i)=local.atom[i].dosckint[1] + local.atom[i].evec[0] * local.atom[i].mcpsc_mt;
    s_buf(2,i)=local.atom[i].dosckint[2] + local.atom[i].evec[1] * local.atom[i].mcpsc_mt;
    s_buf(3,i)=local.atom[i].dosckint[3] + local.atom[i].evec[2] * local.atom[i].mcpsc_mt;
*/
    mag=local.atom[i].mtotws;
    //s_buf(1,i)=local.atom[i].evec[0]*mag;
    *(s_buf+(1*row_size+i))=local.atom[i].evec[0]*mag;
    //s_buf(2,i)=local.atom[i].evec[1]*mag;
    *(s_buf+(2*row_size+i))=local.atom[i].evec[1]*mag;
    //s_buf(3,i)=local.atom[i].evec[2]*mag;
    *(s_buf+(3*row_size+i))=local.atom[i].evec[2]*mag;

  }

  //MPI_Gather(&s_buf(0,0),4*max_num_local,MPI_DOUBLE,&r_buf(0,0,0),4*max_num_local,MPI_DOUBLE,0,comm.comm);
  //Naive Gather implementation
  if(comm.comm.rank=0)
    for(int i=0;i<comm.comm.size;i++)
       shmem_double_get(&r_buf(0,0,0),s_buf,4*max_num_local,comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*i);
  
  shmem_barrier_all();
  

  if(comm.comm.rank==0)
  {
    for(int p=0; p<comm.comm.size; p++)
    {
      for(int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
        {
          int j=r_buf(0,i,p);
          (ev[j])[0]=r_buf(1,i,p); (ev[j])[1]=r_buf(2,i,p); (ev[j])[2]=r_buf(3,i,p);
        }
    }
  }
}

void LSMS::getMag(Real *ev)
{
  // Array3d<Real> r_buf(4,max_num_local,comm.size);
  Real r_buf[4*max_num_local*comm.comm.size];
  // Matrix<Real> s_buf(4,max_num_local);

  Real mag;
  //Real s_buf[4*max_num_local];
  Real* s_buf = (Real*) shmalloc(4*max_num_local*sizeof(Real));

  for(int i=0; i<max_num_local; i++)
    *(s_buf+(4*i))=-1.0;
    // s_buf(0,i)=-1.0;

  for(int i=0; i<local.num_local; i++)
  {
    *(s_buf+(4*i))=Real(local.global_id[i]);
/*
    s_buf[1+4*i]=local.atom[i].dosckint[1] + local.atom[i].evec[0] * local.atom[i].mcpsc_mt;
    s_buf[2+4*i]=local.atom[i].dosckint[2] + local.atom[i].evec[1] * local.atom[i].mcpsc_mt;
    s_buf[3+4*i]=local.atom[i].dosckint[3] + local.atom[i].evec[2] * local.atom[i].mcpsc_mt;
*/
    mag=local.atom[i].mtotws;
    *(s_buf+(1+4*i))=local.atom[i].evec[0]*mag;
    *(s_buf+(2+4*i))=local.atom[i].evec[1]*mag;
    *(s_buf+(3+4*i))=local.atom[i].evec[2]*mag;

  }

  //MPI_Gather(s_buf,4*max_num_local,MPI_DOUBLE,r_buf,4*max_num_local,MPI_DOUBLE,0,comm.comm);
  //Naive Gather implementation
  if(comm.comm.rank=0)
    for(int i=0;i<comm.comm.size;i++)
       shmem_double_get(r_buf,s_buf,4*max_num_local,comm.comm.start_pe+pow(2,comm.comm.logPE_stride)*i);
  
  shmem_barrier_all();
  

  // MPI_Gather(&s_buf(0,0),4*max_num_local,MPI_DOUBLE,&r_buf(0,0,0),100,MPI_DOUBLE,0,comm.comm);

  if(comm.comm.rank==0)
  {
    for(int p=0; p<comm.comm.size; p++)
    {
      // for(int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
      for(int i=0; i<max_num_local && r_buf[4*(i+max_num_local*p)]>=0; i++)
        {
          // int j=r_buf(0,i,p);
          int j=r_buf[4*(i+max_num_local*p)];
          // ev[3*j+0]=r_buf(1,i,p); ev[3*j+1]=r_buf(2,i,p); ev[3*j+2]=r_buf(3,i,p);
          ev[3*j+0]=r_buf[1+4*(i+max_num_local*p)]; ev[3*j+1]=r_buf[2+4*(i+max_num_local*p)];
          ev[3*j+2]=r_buf[3+4*(i+max_num_local*p)];
        }
    }
  }
}

