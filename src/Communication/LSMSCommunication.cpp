//#include <mpi.h>
#include "LSMSCommunication.hpp" // this includes <shmem.h>

extern "C" {
#include "../../instr.h"
}

int const MAXPTS=3051;
int const MAXCORE=30;


int char_size=sizeof(char);
int int_size=sizeof(int);
int double_size=sizeof(double);
int* sync_send_flag;
int* sync_recv_flag;
char *p2p_buf;

enum shmemType { INT, DOUBLE };


 long pSync1[_SHMEM_BCAST_SYNC_SIZE];
 long pSync2[_SHMEM_BCAST_SYNC_SIZE];

 long long pWrk_ll[_SHMEM_REDUCE_SYNC_SIZE];
 int pWrk_i[_SHMEM_REDUCE_SYNC_SIZE];
 double pWrk_d[_SHMEM_REDUCE_SYNC_SIZE];



void initializeCommunication(LSMSCommunication &comm)
{
  //MPI_Init(NULL,NULL);
  //comm.comm=MPI_COMM_WORLD;
  //MPI_Comm_rank(comm.comm, &comm.rank);
  //MPI_Comm_size(comm.comm, &comm.size);
  
  int i;
  start_pes(0);
  allocate_symm_buffers();
  comm.comm.rank = _my_pe();
  comm.comm.size = _num_pes();
  comm.comm.start_pe = 0;
  comm.comm.logPE_stride = 0;
  sync_send_flag=(int*)shmalloc(comm.comm.size*sizeof(int));
  sync_recv_flag=(int*)shmalloc(comm.comm.size*sizeof(int));
  memset(sync_send_flag,0,comm.comm.size*sizeof(int));
  memset(sync_recv_flag,0,comm.comm.size*sizeof(int));


  for (i=0;i<comm.comm.size;i++)
  {
    sync_send_flag[i]=0;
    sync_recv_flag[i]=0;
  }

  shmem_barrier_all();

  for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i += 1)
  {
        pSync1[i] = _SHMEM_SYNC_VALUE;
        pSync2[i] = _SHMEM_SYNC_VALUE;
  }

}

void initializeCommunication(LSMSCommunication &comm, SHMEM_activeset comm_shmem)
{
  //comm.comm=mpiCommunicator;
  //MPI_Comm_rank(comm.comm, &comm.rank);
  //MPI_Comm_size(comm.comm, &comm.size);
  comm.comm.rank = _my_pe();
  comm.comm.size = _num_pes();
  comm.comm.start_pe = comm_shmem.start_pe;
  comm.comm.logPE_stride = comm_shmem.logPE_stride;
}

void finalizeCommunication(void)
{
  deallocate_symm_buffers();
  //MPI_Finalize();
}

void exitLSMS(LSMSCommunication &comm, int errorCode)
{
  //MPI_Abort(comm.comm, errorCode);
  exit(1);
}

void allocate_symm_buffers()
{
  int s=sizeof(AtomData)+sizeof(Real)*(2*3*MAXPTS+2*MAXCORE)+sizeof(int)*3*2*MAXCORE+sizeof(int)*2;
  p2p_buf = (char*)shmalloc(s);
}

void deallocate_symm_buffers()
{
  shfree(p2p_buf);
}

void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms, 
                           CrystalParameters &crystal, MixingParameters &mix)
{
  int const s=sizeof(LSMSSystemParameters)+9*sizeof(Real)+sizeof(int)+10
    +sizeof(MixingParameters)+5*sizeof(int);
  int rem=0,ele=0;
  int tot_bufsize=s;
  rem=s%32;
  ele=s/32;
  if  (rem!=0)
  {
    tot_bufsize=s-rem+32;
    ele++;
  }
  // TODO fine-tune this size
  tot_bufsize=65536;
  char* buf=(char*)shmalloc(tot_bufsize);
  int pos=0;
  int sec_id;

  if(comm.comm.rank==0)
  {
    
    //MPI_Pack(lsms.systemid,80,MPI_CHAR,buf,s,&pos,comm.comm);
    //MPI_Pack(lsms.title,80,MPI_CHAR,buf,s,&pos,comm.comm);
    //MPI_Pack(lsms.potential_file_in,128,MPI_CHAR,buf,s,&pos,comm.comm);
    //MPI_Pack(lsms.potential_file_out,128,MPI_CHAR,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.pot_in_type,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.pot_out_type,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.num_atoms,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.nspin,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.nrel_rel,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.nrelc,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.nrelv,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.n_spin_cant,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.n_spin_pola,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.mtasa,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.fixRMT,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.nscf,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.writeSteps,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.clight,1,MPI_DOUBLE,buf,s,&pos,comm.comm);

    //MPI_Pack(&lsms.energyContour.grid,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.npts,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.ebot,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.etop,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.eibot,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.eitop,1,MPI_DOUBLE,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.energyContour.maxGroupSize,1,MPI_INT,buf,s,&pos,comm.comm);

    //MPI_Pack(&lsms.mixing,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.alphaDV,1,MPI_DOUBLE,buf,s,&pos,comm.comm);

    //MPI_Pack(&lsms.global.iprint,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.global.print_node,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.global.default_iprint,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.global.istop,32,MPI_CHAR,buf,s,&pos,comm.comm);
    //MPI_Pack(&lsms.global.GPUThreads,32,MPI_INT,buf,s,&pos,comm.comm);

    //MPI_Pack(&crystal.num_types,1,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&crystal.bravais(0,0),9,MPI_DOUBLE,buf,s,&pos,comm.comm);

    //************  MemCpying  ***************
    memcpy(&buf[pos],&lsms.systemid,80*char_size); pos = pos+80*char_size;
    memcpy(&buf[pos],&lsms.title,80*char_size); pos = pos+80*char_size;
    memcpy(&buf[pos],&lsms.potential_file_in,128*char_size); pos = pos+128*char_size;
    memcpy(&buf[pos],&lsms.potential_file_out,128*char_size); pos = pos+128*char_size;
    memcpy(&buf[pos],&lsms.pot_in_type,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.pot_out_type ,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.num_atoms,int_size); pos = pos+int_size;

    memcpy(&buf[pos],&lsms.nspin,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.nrel_rel,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.nrelc,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.nrelv,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.n_spin_cant,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.n_spin_pola,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.mtasa,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.fixRMT,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.nscf,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.writeSteps,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.clight,double_size); pos = pos+double_size;

    memcpy(&buf[pos],&lsms.energyContour.grid,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.energyContour.npts,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.energyContour.ebot,double_size); pos = pos+double_size;
    memcpy(&buf[pos],&lsms.energyContour.etop,double_size); pos = pos+double_size;
    memcpy(&buf[pos],&lsms.energyContour.eibot,double_size); pos = pos+double_size;
    memcpy(&buf[pos],&lsms.energyContour.eitop,double_size); pos = pos+double_size;
    memcpy(&buf[pos],&lsms.energyContour.maxGroupSize,int_size); pos = pos+int_size;

    memcpy(&buf[pos],&lsms.mixing,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.alphaDV,double_size); pos = pos+double_size;

    memcpy(&buf[pos],&lsms.global.iprint,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.global.print_node,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.global.default_iprint,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&lsms.global.istop,32*char_size); pos = pos+32*char_size;
    memcpy(&buf[pos],&lsms.global.GPUThreads,32*int_size); pos = pos+32*int_size;

    memcpy(&buf[pos],&crystal.num_types,int_size); pos = pos+int_size;
    memcpy(&buf[pos],&crystal.bravais(0,0),9*double_size); pos = pos+9*double_size;


// MixingParameters
    // MPI_CXX_BOOL is not always available
    // MPI_Pack(&mix.quantity[0],mix.numQuantities,MPI_CXX_BOOL,buf,s,&pos,comm.comm);
    // copy to temporary int array and send this
    int tmpQuantity[mix.numQuantities];
    for(int i=0; i<mix.numQuantities; i++)
      if(mix.quantity[i])
        tmpQuantity[i] = 1;
      else
        tmpQuantity[i] = 0; 
    //MPI_Pack(&tmpQuantity[0],mix.numQuantities,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&mix.algorithm[0],mix.numQuantities,MPI_INT,buf,s,&pos,comm.comm);
    //MPI_Pack(&mix.mixingParameter[0],mix.numQuantities,MPI_DOUBLE,buf,s,&pos,comm.comm);
    memcpy(&buf[pos],&tmpQuantity[0],mix.numQuantities*int_size); pos = pos+mix.numQuantities*int_size;
    memcpy(&buf[pos],&mix.algorithm[0],mix.numQuantities*int_size); pos = pos+mix.numQuantities*int_size;
    memcpy(&buf[pos],&mix.mixingParameter[0],mix.numQuantities*double_size); pos = pos+mix.numQuantities*double_size;

  }
  //MPI_Bcast(buf,s,MPI_PACKED,0,comm.comm);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  shmem_broadcast32(&buf[0], &buf[0], tot_bufsize, 0, 0, 0, comm.comm.size,pSync1);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  if(comm.comm.rank!=0)
  {
    int pos=0;
    //MPI_Unpack(buf,s,&pos,lsms.systemid,80,MPI_CHAR,comm.comm);
    //MPI_Unpack(buf,s,&pos,lsms.title,80,MPI_CHAR,comm.comm);
    //MPI_Unpack(buf,s,&pos,lsms.potential_file_in,128,MPI_CHAR,comm.comm);
    //MPI_Unpack(buf,s,&pos,lsms.potential_file_out,128,MPI_CHAR,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.pot_in_type,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.pot_out_type,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.num_atoms,1,MPI_INT,comm.comm);
    memcpy(&lsms.systemid,&buf[pos],80*char_size); pos = pos+80*char_size;
    memcpy(&lsms.title,&buf[pos],80*char_size); pos = pos+80*char_size;
    memcpy(&lsms.potential_file_in,&buf[pos],128*char_size); pos = pos+128*char_size;
    memcpy(&lsms.potential_file_out,&buf[pos],128*char_size); pos = pos+128*char_size;
    memcpy(&lsms.pot_in_type,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.pot_out_type,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.num_atoms,&buf[pos],int_size); pos = pos+int_size;
    crystal.num_atoms=lsms.num_atoms;
    //MPI_Unpack(buf,s,&pos,&lsms.nspin,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.nrel_rel,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.nrelc,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.nrelv,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.n_spin_cant,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.n_spin_pola,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.mtasa,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.fixRMT,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.nscf,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.writeSteps,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.clight,1,MPI_DOUBLE,comm.comm);

    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.grid,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.npts,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.ebot,1,MPI_DOUBLE,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.etop,1,MPI_DOUBLE,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.eibot,1,MPI_DOUBLE,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.eitop,1,MPI_DOUBLE,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.energyContour.maxGroupSize,1,MPI_INT,comm.comm);

    //MPI_Unpack(buf,s,&pos,&lsms.mixing,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.alphaDV,1,MPI_DOUBLE,comm.comm);

    //MPI_Unpack(buf,s,&pos,&lsms.global.iprint,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.global.print_node,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.global.default_iprint,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.global.istop,32,MPI_CHAR,comm.comm);
    //MPI_Unpack(buf,s,&pos,&lsms.global.GPUThreads,32,MPI_INT,comm.comm);

    //MPI_Unpack(buf,s,&pos,&crystal.num_types,1,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&crystal.bravais(0,0),9,MPI_DOUBLE,comm.comm);

    memcpy(&lsms.nspin,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.nrel_rel,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.nrelc,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.nrelv,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.n_spin_cant,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.n_spin_pola,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.mtasa,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.fixRMT,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.nscf,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.writeSteps,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.clight,&buf[pos],double_size); pos = pos+double_size;

    memcpy(&lsms.energyContour.grid,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.energyContour.npts,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.energyContour.ebot,&buf[pos],double_size); pos = pos+double_size;
    memcpy(&lsms.energyContour.etop,&buf[pos],double_size); pos = pos+double_size;
    memcpy(&lsms.energyContour.eibot,&buf[pos],double_size); pos = pos+double_size;
    memcpy(&lsms.energyContour.eitop,&buf[pos],double_size); pos = pos+double_size;
    memcpy(&lsms.energyContour.maxGroupSize,&buf[pos],int_size); pos = pos+int_size;

    memcpy(&lsms.mixing,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.alphaDV,&buf[pos],double_size); pos = pos+double_size;

    memcpy(&lsms.global.iprint,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.global.print_node,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.global.default_iprint,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&lsms.global.istop,&buf[pos],32*char_size); pos = pos+32*char_size;
    memcpy(&lsms.global.GPUThreads,&buf[pos],32*int_size); pos = pos+32*int_size;

    memcpy(&crystal.num_types,&buf[pos],int_size); pos = pos+int_size;
    memcpy(&crystal.bravais(0,0),&buf[pos],9*double_size); pos = pos+9*double_size;

    crystal.resize(crystal.num_atoms);
    crystal.resizeTypes(crystal.num_types);


// MixingParameters
    // MPI_CXX_BOOL is not always available
    // MPI_Unpack(buf,s,&pos,&mix.quantity[0],mix.numQuantities,MPI_CXX_BOOL,comm.comm);
    // recieve temporary int array and copy
    int tmpQuantity[mix.numQuantities];
    //MPI_Unpack(buf,s,&pos,&tmpQuantity[0],mix.numQuantities,MPI_INT,comm.comm);
    memcpy(&tmpQuantity[0],&buf[pos],mix.numQuantities*int_size); pos = pos+mix.numQuantities*int_size;

    for(int i=0; i<mix.numQuantities; i++)
      if(tmpQuantity[i]==1)
        mix.quantity[i] = true;
      else
        mix.quantity[i] = false; 
    //MPI_Unpack(buf,s,&pos,&mix.algorithm[0],mix.numQuantities,MPI_INT,comm.comm);
    //MPI_Unpack(buf,s,&pos,&mix.mixingParameter[0],mix.numQuantities,MPI_DOUBLE,comm.comm);
    memcpy(&mix.algorithm[0],&buf[pos],mix.numQuantities*int_size); pos = pos+mix.numQuantities*int_size;
    memcpy(&mix.mixingParameter[0],&buf[pos],mix.numQuantities*double_size); pos = pos+mix.numQuantities*double_size;
  }

 for(int i=0; i<mix.numQuantities; i++)
      printf("mix.quantity[%d]=%d\n", i,mix.quantity[i]);

  // Allocate buffer for transmitting Crystal params
  int buff_size;

  if((crystal.num_types*sizeof(AtomType)) > (3*crystal.num_atoms*double_size))
     buff_size = crystal.num_types*sizeof(AtomType);
  else 
     buff_size = 3*crystal.num_atoms*double_size;  
 
  shfree(buf);
  // TODO finetune buff-size
  buff_size=1048576; //sizeof(LSMSSystemParameters)+9*sizeof(Real);
  rem=buff_size%64;
  ele=buff_size/64;
  if(rem != 0)
  {
     buff_size=buff_size-rem+64;
     ele++;
  }

  double* temp_buff=(double*) shmalloc(buff_size);
  int*    temp_intbuff=(int*) shmalloc(buff_size);

  //MPI_Bcast(&crystal.position(0,0),3*crystal.num_atoms,MPI_DOUBLE,0,comm.comm);
//TODO check if a barrier is neededa after broadcast ... data not updated otherwise
  if(comm.comm.rank == 0)
      memcpy(temp_buff,&crystal.position(0,0),3*crystal.num_atoms*double_size);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  shmem_broadcast64(temp_buff, temp_buff,3*crystal.num_atoms, 0, 0, 0, comm.comm.size,pSync1);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  if(comm.comm.rank != 0)
      memcpy(&crystal.position(0,0),temp_buff,3*crystal.num_atoms*double_size);

  //MPI_Bcast(&crystal.evecs(0,0),3*crystal.num_atoms,MPI_DOUBLE,0,comm.comm);
  if(comm.comm.rank == 0){
      memcpy(temp_buff,&crystal.evecs(0,0),3*crystal.num_atoms*double_size);
}
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  shmem_broadcast64(temp_buff, temp_buff, 3*crystal.num_atoms, 0, 0, 0, comm.comm.size,pSync1);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  if(comm.comm.rank != 0){
      memcpy(&crystal.evecs(0,0),temp_buff,3*crystal.num_atoms*double_size);
}

  //MPI_Bcast(&crystal.type[0],crystal.num_atoms,MPI_INT,0,comm.comm);
  if(comm.comm.rank == 0){
      memcpy(temp_intbuff,&crystal.type[0],crystal.num_atoms*int_size);
  }
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  shmem_broadcast32(temp_intbuff, temp_intbuff, crystal.num_atoms, 0, 0, 0, comm.comm.size,pSync1);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  if(comm.comm.rank != 0){
      memcpy(&crystal.type[0],temp_intbuff,crystal.num_atoms*int_size);
  }

// This is dangerous and assumes homogeneous nodes:
  //MPI_Bcast(&crystal.types[0],crystal.num_types*sizeof(AtomType),MPI_BYTE,0,comm.comm);
  if(comm.comm.rank == 0)
      memcpy(temp_buff,&crystal.types[0],crystal.num_types*sizeof(AtomType));
  // having to use the smallest possible broadcast:"32"-type
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  shmem_broadcast32(temp_buff,temp_buff,crystal.num_types*sizeof(AtomType)/4,0,0,0,comm.comm.size,pSync1);
  shmem_barrier(0, 0, comm.comm.size,pSync2);
  if(comm.comm.rank != 0)
      memcpy(&crystal.types[0],temp_buff,crystal.num_types*sizeof(AtomType));

  shmem_barrier(0, 0, comm.comm.size,pSync1);
  shfree(temp_buff);
  shfree(temp_intbuff);

// get maximum lmax
  crystal.maxlmax=0;
  for(int i=0; i<crystal.num_types; i++)
    if(crystal.types[i].lmax>crystal.maxlmax) crystal.maxlmax=crystal.types[i].lmax; 
  lsms.maxlmax=crystal.maxlmax;
}

void communicateSingleAtomData(LSMSCommunication &comm, int from, int to, int &local_id, AtomData &atom, int tag)
{
  //The buffers used in this func are pre-allocated within initializeCommunication() of size 's' below 
  //int s=sizeof(AtomData)+sizeof(Real)*(2*3*MAXPTS+2*MAXCORE)+sizeof(int)*3*2*MAXCORE+sizeof(int);
  // 304 bytes transferred in each of the ITER_MAX iterations
  const int maxPts=MAXPTS;
  const int maxCore=MAXCORE;
  int t,i;
  static int count=0;
  const int ITER_MAX=1;
  int sec_id;
  char instr_name[100];
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);

  sprintf(instr_name,"SA_%s_%d_%d",hostname,from,to);
  if(comm.comm.rank != 0) 
    INITFREQ("2901000");
  if(comm.comm.rank >= PE_ON_FPGA)
    sec_id = START(instr_name,ITER_MAX); // RAPL instr
  

  if(comm.comm.rank==from)
  {

   for (i=0;i<ITER_MAX;i++){
    int pos=0;

    memcpy(&p2p_buf[pos],&local_id,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.jmt,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.jws,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.xstart,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],&atom.rmt,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],atom.header,80*char_size); pos+=80*char_size;
    memcpy(&p2p_buf[pos],&atom.alat,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],&atom.efermi,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],&atom.vdif,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],&atom.ztotss,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],&atom.zcorss,double_size); pos+=double_size;
    memcpy(&p2p_buf[pos],atom.evec,3*double_size); pos+=3*double_size;
    memcpy(&p2p_buf[pos],&atom.nspin,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.numc,int_size); pos+=int_size;

    t=atom.vr.n_row();

    memcpy(&p2p_buf[pos],&t,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.vr(0,0),2*t*double_size); pos+=2*t*double_size;
    memcpy(&p2p_buf[pos],&atom.rhotot(0,0),2*t*double_size); pos+=2*t*double_size;
    memcpy(&p2p_buf[pos],&atom.corden(0,0),2*t*double_size); pos+=2*t*double_size;

    t=atom.ec.n_row();

    memcpy(&p2p_buf[pos],&t,int_size); pos+=int_size;
    memcpy(&p2p_buf[pos],&atom.ec(0,0),2*t*double_size); pos+=2*t*double_size;
    memcpy(&p2p_buf[pos],&atom.nc(0,0),2*t*int_size); pos+=2*t*int_size;
    memcpy(&p2p_buf[pos],&atom.lc(0,0),2*t*int_size); pos+=2*t*int_size;
    memcpy(&p2p_buf[pos],&atom.kc(0,0),2*t*int_size); pos+=2*t*int_size;

    shmem_int_wait_until((sync_send_flag+to),SHMEM_CMP_EQ,1);
    shmem_putmem(p2p_buf, p2p_buf, 1048576, to);
    shmem_int_add((sync_send_flag+to),-1,comm.comm.rank);
    shmem_int_add((sync_recv_flag+comm.comm.rank),1,to);
    shmem_quiet();

   }// end of false for loop
    
  }
  if(comm.comm.rank==to)
  {
for(i=0;i<ITER_MAX;i++) {
    int pos=0;

    sync_recv_flag[from]=0;
    shmem_int_add((sync_send_flag+comm.comm.rank),1,from);
    shmem_quiet();
    shmem_int_wait_until((sync_recv_flag+from),SHMEM_CMP_EQ,1);
  
    memcpy(&local_id,&p2p_buf[pos],int_size); pos+=int_size;
    memcpy(&atom.jmt,&p2p_buf[pos],int_size); pos+=int_size;
    memcpy(&atom.jws,&p2p_buf[pos],int_size); pos+=int_size;
    memcpy(&atom.xstart,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(&atom.rmt,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(atom.header,&p2p_buf[pos],80*char_size); pos+=80*char_size;
    memcpy(&atom.alat,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(&atom.efermi,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(&atom.vdif,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(&atom.ztotss,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(&atom.zcorss,&p2p_buf[pos],double_size); pos+=double_size;
    memcpy(atom.evec,&p2p_buf[pos],3*double_size); pos+=3*double_size;
    memcpy(&atom.nspin,&p2p_buf[pos],int_size); pos+=int_size;
    memcpy(&atom.numc,&p2p_buf[pos],int_size); pos+=int_size;

    memcpy(&t,&p2p_buf[pos],int_size); pos+=int_size;

    if(t!=atom.vr.n_row()) atom.resizePotential(t);

    memcpy(&atom.vr(0,0),&p2p_buf[pos],2*t*double_size); pos+=2*t*double_size;
    memcpy(&atom.rhotot(0,0),&p2p_buf[pos],2*t*double_size); pos+=2*t*double_size;
    memcpy(&atom.corden(0,0),&p2p_buf[pos],2*t*double_size); pos+=2*t*double_size;
    memcpy(&t,&p2p_buf[pos],int_size); pos+=int_size;

    if(t!=atom.nc.n_row()) atom.resizeCore(t);

    memcpy(&atom.ec(0,0),&p2p_buf[pos],2*t*double_size); pos+=2*t*double_size;
    memcpy(&atom.nc(0,0),&p2p_buf[pos],2*t*int_size); pos+=2*t*int_size;
    memcpy(&atom.lc(0,0),&p2p_buf[pos],2*t*int_size); pos+=2*t*int_size;
    memcpy(&atom.kc(0,0),&p2p_buf[pos],2*t*int_size); pos+=2*t*int_size;
    shmem_int_add((sync_recv_flag+from),-1,comm.comm.rank);
    shmem_quiet();    
   }
  }

 if(comm.comm.rank != 0)
     INITFREQ("2901000");
 if(comm.comm.rank >= PE_ON_FPGA)
     END(sec_id);
 
}

void expectTmatCommunication(LSMSCommunication &comm, LocalTypeInfo &local)
{

// prepost all recieves for tmats from remote nodes
/*
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    int from=comm.tmatFrom[i].remoteNode;
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
    {
      // printf("Node %d: expect tmat %d from %d\n",comm.rank,comm.tmatFrom[i].globalIdx[j],from);
      MPI_Irecv(&local.tmatStore(0,comm.tmatFrom[i].tmatStoreIdx[j]),2*local.lDimTmatStore,
                MPI_DOUBLE,from,comm.tmatFrom[i].globalIdx[j],comm.comm,
                &comm.tmatFrom[i].communicationRequest[j]);

    }
  }
*/
}

void sendTmats(LSMSCommunication &comm, LocalTypeInfo &local)
{

  void * temp_buff=(void*)shmalloc(2*local.lDimTmatStore*double_size);
  for(int i=0; i<comm.numTmatTo; i++)
  {
    int to=comm.tmatTo[i].remoteNode;
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
    {
      // printf("Node %d: send tmat %d to %d\n",comm.rank,comm.tmatTo[i].globalIdx[j],to);
      /*
      MPI_Isend(&local.tmatStore(0,comm.tmatTo[i].tmatStoreIdx[j]),2*local.lDimTmatStore,
                MPI_DOUBLE,to,comm.tmatTo[i].globalIdx[j],comm.comm,
                &comm.tmatTo[i].communicationRequest[j]);
      */
      // Assuming comm.tmatTo[i].numTmats == comm.tmatFrom[i].numTmats
      // If not ... check the iteration space mapping between 
      //    the receivers (above) and senders (here)
      memcpy(temp_buff,&local.tmatStore(0,comm.tmatFrom[i].tmatStoreIdx[j]),2*local.lDimTmatStore*double_size);
      // Note: there is a trade-off here: either reuse the same 
      // buffer and enjoy reduced cache misses (and potential 
      // power savings due to less data xfers between memories) 
      // or enjoy non-blocking communication with different buffers
      // but at the cost of multiple cache misses (and potentially
      // high power consumption)
      shmem_putmem(temp_buff, &local.tmatStore(0,comm.tmatTo[i].tmatStoreIdx[j]), 2*local.lDimTmatStore*double_size, to);
    }
  }
}
void finalizeTmatCommunication(LSMSCommunication &comm)
{
      // TODO: make it non-blocking ... beware of the memcpy & temp_buff
      shmem_quiet();
  /*
  MPI_Status status;
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    int from=comm.tmatFrom[i].remoteNode;
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
    {
      // printf("Finalize recieve request %d from node %d\n",j,from);
      MPI_Wait(&comm.tmatFrom[i].communicationRequest[j],&status);
    }
  }
  for(int i=0; i<comm.numTmatTo; i++)
  {
    int to=comm.tmatTo[i].remoteNode;
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
    {
      // printf("Finalize send request %d to node %d\n",j,to);
      MPI_Wait(&comm.tmatTo[i].communicationRequest[j],&status);
    }
  }
  */
}

void printCommunicationInfo(FILE *f, LSMSCommunication &comm)
{
  fprintf(f,"Communication: rank no. %d of %d\n",comm.comm.rank,comm.comm.size);
  fprintf(f,"Sending tmats to %d remote nodes:\n",comm.numTmatTo);
  for(int i=0; i<comm.numTmatTo; i++)
  {
    fprintf(f,"Node %d :",comm.tmatTo[i].remoteNode);
    for(int j=0; j<comm.tmatTo[i].numTmats; j++)
      fprintf(f," %d[%d]",comm.tmatTo[i].globalIdx[j],comm.tmatTo[i].tmatStoreIdx[j]);
    fprintf(f,"\n");
  }
  fprintf(f,"Recieving tmats from %d remote nodes:\n",comm.numTmatFrom);
  for(int i=0; i<comm.numTmatFrom; i++)
  {
    fprintf(f,"Node %d :",comm.tmatFrom[i].remoteNode);
    for(int j=0; j<comm.tmatFrom[i].numTmats; j++)
      fprintf(f," %d[%d]",comm.tmatFrom[i].globalIdx[j],comm.tmatFrom[i].tmatStoreIdx[j]);
    fprintf(f,"\n");
  }
}

/*inline*/ void globalMax_int(LSMSCommunication &comm,int &a)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  static int r_i;  
  r_i=a;
  shmem_int_max_to_all(&(a), &r_i, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_i, pSync1);
}

/*inline*/ void globalMax_double(LSMSCommunication &comm,double &a)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  static double r_d;  
  r_d=a;
  shmem_double_max_to_all(&(a), &r_d, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync1);
}

/*inline*/ void globalSum_int(LSMSCommunication &comm,int &a)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  static int r_i;  
  r_i=a;
  shmem_int_sum_to_all(&a, &r_i, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_i, pSync2);
}

/*inline*/ void globalSum_double(LSMSCommunication &comm,double &a)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  static double r_d;  
  r_d=a;
  shmem_double_sum_to_all(&a, &r_d, 1,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
}

/*inline*/ void globalSum_double(LSMSCommunication &comm,double &a, int n)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  static double r_d;  
  r_d=a;
  shmem_double_sum_to_all(&a, &r_d, n,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
}

/*inline*/ void globalSum_real(LSMSCommunication &comm,double *a, int n)
{
  shmem_barrier(comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size,pSync1);
  double* r_d = (double*)shmalloc(n*sizeof(double));  
  memcpy(r_d,a,n*sizeof(double));
  shmem_double_sum_to_all(a, r_d, n,comm.comm.start_pe, comm.comm.logPE_stride, comm.comm.size, pWrk_d, pSync2);
  shfree(r_d);
}

