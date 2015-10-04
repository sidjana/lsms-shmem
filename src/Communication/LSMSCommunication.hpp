#ifndef LSMSCOMMUNICATION_H
#define LSMSCOMMUNICATION_H

#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include <vector>

#include "Main/SystemParameters.hpp"
#include "Main/mixing.hpp"
#include "SingleSite/AtomData.hpp"
#include <shmem.h>



typedef double Real;

// 'pSync' variables with 'long' data-type are numbered
// All other variables for 'pSync's and 'pWrk's end with // _<data-type_abbrv>

extern long pSync1[_SHMEM_BCAST_SYNC_SIZE];
extern long pSync2[_SHMEM_BCAST_SYNC_SIZE];

extern long long pWrk_ll[_SHMEM_REDUCE_SYNC_SIZE];
extern int pWrk_i[_SHMEM_REDUCE_SYNC_SIZE];
extern double pWrk_d[_SHMEM_REDUCE_SYNC_SIZE];


typedef double Real;

class TmatCommType {
public:
  int remoteNode;
  int numTmats;
  std::vector<int> tmatStoreIdx;
  std::vector<int> globalIdx;
 // std::vector<MPI_Request> communicationRequest;
};

struct SHMEM_activeset 
{
  int rank;
  int size;
  int start_pe; 
  int logPE_stride;
};

class LSMSCommunication {
public:
  //MPI_Comm comm;
  SHMEM_activeset comm;
  int numTmatTo, numTmatFrom;
  std::vector<TmatCommType> tmatTo, tmatFrom;
};

void initializeCommunication(LSMSCommunication &comm);
//void initializeCommunication(LSMSCommunication &comm, MPI_Comm mpiCommunicator);
void initializeCommunication(LSMSCommunication &comm, SHMEM_activeset comm_shmem);
void finalizeCommunication(void);
void exitLSMS(LSMSCommunication &comm, int errorCode);


void allocate_symm_buffers();
void deallocate_symm_buffers();


void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                           CrystalParameters &crystal, MixingParameters &mix);

void communicateSingleAtomData(LSMSCommunication &comm, int from, int to,
                               int &local_id, AtomData &atom, int tag=0);

void expectTmatCommunication(LSMSCommunication &comm, LocalTypeInfo &local);
void sendTmats(LSMSCommunication &comm, LocalTypeInfo &local);
void finalizeTmatCommunication(LSMSCommunication &comm);

void printCommunicationInfo(FILE *f, LSMSCommunication &comm);


/*inline*/ void globalMax_int(LSMSCommunication &comm,int &a);
/*inline*/ void globalMax_double(LSMSCommunication &comm,double &a);
/*inline*/ void globalSum_int(LSMSCommunication &comm,int &a);
/*inline*/ void globalSum_double(LSMSCommunication &comm,double &a);
/*inline*/ void globalSum_double(LSMSCommunication &comm,double &a, int n);
/*inline*/ void globalSum_real(LSMSCommunication &comm,double *a, int n);

// 
// /*inline*/
// void globalMax_int(LSMSCommunication &comm,int &a)
// {
//   static int r_i;  
//   r_i=a;
//   shmem_int_max_to_all(&(a), &r_i, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_i, pSync1);
// }
// 
// /*inline*/
// void globalMax_double(LSMSCommunication &comm,double &a)
// {
//   static double r_d;  
//   r_d=a;
//   shmem_double_max_to_all(&(a), &r_d, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync1);
// }
// 
// /*inline*/
// void globalSum_int(LSMSCommunication &comm,int &a)
// {
//   static int r_i;  
//   r_i=a;
//   shmem_int_sum_to_all(&a, &r_i, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_i, pSync2);
// }
// 
// /*inline*/
// void globalSum_double(LSMSCommunication &comm,double &a)
// {
//   static double r_d;  
//   r_d=a;
//   shmem_double_sum_to_all(&a, &r_d, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
// }
// 
// /*inline*/
// void globalSum_double(LSMSCommunication &comm,double &a, int n)
// {
//   static double r_d;  
//   r_d=a;
//   shmem_double_sum_to_all(&a, &r_d, n,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
// }
// 
// /*inline*/
// void globalSum_real(LSMSCommunication &comm,double *a, int n)
// {
//   double* r_d = (double*)shmalloc(n*sizeof(double));  
//   memcpy(r_d,a,n*sizeof(double));
//   shmem_double_sum_to_all(a, r_d, n,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
// }
// 
// 

#endif
