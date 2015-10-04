// main driver for LSMS_3 class

#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL/VORPOL.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "lsmsClass.hpp"
#include "EvecGenerator.h"
#include "WangLandau.h"
// #include "WangLandau_withoutKernel.h"
#include "ExhaustiveIsing.h"
#include "WangLandau2d.h"

// #define USE_PAPI 1

#ifdef USE_PAPI
#include <papi.h>
#endif

#define R_VALUE_OFFSET 2

LocalTypeInfo LSMS::local;

int main(int argc, char *argv[])
{
  int size, rank, world_rank, my_group;
  int num_lsms; // number of parallel LSMS instances
  int size_lsms; // number of atoms in a lsms instance
  int num_steps; // number of energy calculations
  int initial_steps; // number of steps before sampling starts
  int stepCount=0; // count the Monte Carlo steps executed
  double max_time; // maximum walltime for this run in seconds
  bool restrict_time = false;       // was the maximum time specified?
  bool restrict_steps = false; // or the max. numer of steps?
  int align; // alignment of lsms_instances
  
  double magnetization;
  double energy_accumulator; // accumulates the enegy to calculate the mean
  int energies_accumulated;


  int new_peid,new_root;
  static int op,flag;
  double *evec,*r_values;
  evec=(double *)shmalloc(sizeof(double)*3*size_lsms);
  r_values=(double *)shmalloc(sizeof(double)*(R_VALUE_OFFSET+3*(size_lsms+1)));




  energy_accumulator=0.0;
  energies_accumulated=0;

  double walltime_0,walltime;

  double restartWriteFrequency=30.0*60.0;
  double nextWriteTime=restartWriteFrequency;

  MPI_Comm local_comm;
  int *lsms_rank0;
  MPI_Status status;

  char prefix[40];
  char i_lsms_name[64];
  char gWL_in_name[64], gWL_out_name[64];
  char mode_name[64];
  char energy_calculation_name[64];
  char stupid[37];

  char step_out_name[64];
  char wl_step_out_name[128];
  char *wl_stepf=NULL;
  bool step_out_flag=false;
  std::ofstream step_out_file;
  typedef enum {Constant, Random, WangLandau_1d, ExhaustiveIsing, WangLandau_2d} EvecGenerationMode;
  typedef enum {MagneticMoment, MagneticMomentZ, MagneticMomentX, MagneticMomentY} SecondDimension;

  EvecGenerationMode evec_generation_mode = Constant;
  SecondDimension second_dimension = MagneticMoment;
  double ev0[3];

  bool return_moments_flag=true; // true-> return all magnetic moments from lsms run at each step.
  bool generator_needs_moment=false;

  typedef enum {OneStepEnergy, MultiStepEnergy, ScfEnergy} EnergyCalculationMode;
  EnergyCalculationMode energyCalculationMode = OneStepEnergy;
  int energyIndex=1; // index for the return value to use for the MC step (0: total energy, 1: band energy)

  ev0[0]=ev0[1]=0.0; ev0[2]=1.0;
  // size has to be align + size_lsms*num_lsms
  align=1;
  num_lsms=1;
  size_lsms=-1;
  my_group=-1;
  num_steps=1;
  initial_steps=0;

  sprintf(i_lsms_name,"i_lsms");
  gWL_in_name[0]=gWL_out_name[0]=0;
  mode_name[0]=0;
  energy_calculation_name[0]=0;

  // check command line arguments
  for(int i=0; i<argc; i++)
  {
    if(!strcmp("-num_lsms",argv[i])) num_lsms=atoi(argv[++i]);
    if(!strcmp("-size_lsms",argv[i])) size_lsms=atoi(argv[++i]);
    if(!strcmp("-align",argv[i])) align=atoi(argv[++i]);
    if(!strcmp("-num_steps",argv[i])) {num_steps=atoi(argv[++i]); restrict_steps=true;}
    if(!strcmp("-initial_steps",argv[i])) initial_steps=atoi(argv[++i]); 
    if(!strcmp("-walltime",argv[i])) {max_time=60.0*atof(argv[++i]); restrict_time=true;}
    if(!strcmp("-i",argv[i])) strncpy(i_lsms_name,argv[++i],64);
    if(!strcmp("-random_dir",argv[i])) {evec_generation_mode = Random;}
    if(!strcmp("-step_out",argv[i]))
    {strncpy(step_out_name,argv[++i],64); step_out_flag=true;
      return_moments_flag=true;}
    if(!strcmp("-wl_out", argv[i])) strncpy(gWL_out_name,argv[++i],64);
    if(!strcmp("-wl_in", argv[i])) strncpy(gWL_in_name,argv[++i],64);
    if(!strcmp("-mode", argv[i])) strncpy(mode_name,argv[++i],64);
    if(!strcmp("-energy_calculation",argv[i])) strncpy(energy_calculation_name,argv[++i],64);
  }

  if(!(restrict_steps || restrict_time)) restrict_steps=true;

  if(mode_name[0]!=0)
  {
    if(!strcmp("constant",mode_name)) evec_generation_mode = Constant;
    if(!strcmp("random",mode_name)) evec_generation_mode = Random;
    if(!strcmp("1d",mode_name)) evec_generation_mode = WangLandau_1d;
    if(!strcmp("ising",mode_name)) evec_generation_mode = ExhaustiveIsing;
    if(!strcmp("2d",mode_name)) evec_generation_mode = WangLandau_2d;
    if(!strcmp("2d-m",mode_name)) {evec_generation_mode = WangLandau_2d; second_dimension=MagneticMoment;}
    if(!strcmp("2d-x",mode_name)) {evec_generation_mode = WangLandau_2d; second_dimension=MagneticMomentX;}
    if(!strcmp("2d-y",mode_name)) {evec_generation_mode = WangLandau_2d; second_dimension=MagneticMomentY;}
    if(!strcmp("2d-z",mode_name)) {evec_generation_mode = WangLandau_2d; second_dimension=MagneticMomentZ;}
  }

  if(energy_calculation_name[0]!=0)
  {
    if(energy_calculation_name[0]=='o') { energyCalculationMode = OneStepEnergy; energyIndex=1; }
    if(energy_calculation_name[0]=='m') { energyCalculationMode = MultiStepEnergy; energyIndex=1; }
    if(energy_calculation_name[0]=='s') { energyCalculationMode = ScfEnergy; energyIndex=0; }
  }

#ifdef USE_PAPI
#define NUM_PAPI_EVENTS 4
  int hw_counters = PAPI_num_counters();
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS]; // = {PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_FP_OPS,PAPI_VEC_INS};
  char *papi_event_name[] = {"PAPI_TOT_INS","PAPI_FP_OPS",
                             "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:OP_TYPE",
                             "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:OP_TYPE"};
  // "RETIRED_INSTRUCTIONS",
  // "RETIRED_MMX_AND_FP_INSTRUCTIONS:PACKED_SSE_AND_SSE2",
  // "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:1",
  // "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:1"
  // get events from names:
  for(int i=0; i<NUM_PAPI_EVENTS; i++)
  {
    if(PAPI_event_name_to_code(papi_event_name[i],&papi_events[i]) != PAPI_OK)
    {
      // printline("Error in obtaining PAPI event code for: "+ttos(papi_event_name[i]),
      //           std::cerr,parameters.myrankWorld);
      // printline("Skipping all following events",
      //           std::cerr,parameters.myrankWorld);
      if(hw_counters>i) hw_counters=i;
    }
  }
  long long papi_values[NUM_PAPI_EVENTS+4];
  // printline("PAPI: "+ttos(hw_counters)+" counters available",std::cout,parameters.myrankWorld);
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  long long papi_real_cyc_0 = PAPI_get_real_cyc();
  long long papi_real_usec_0 = PAPI_get_real_usec();
  long long papi_virt_cyc_0 = PAPI_get_virt_cyc();
  long long papi_virt_usec_0 = PAPI_get_virt_usec();
  PAPI_start_counters(papi_events,hw_counters);
#endif


  lsms_rank0=(int *)malloc(sizeof(int)*(num_lsms+1));

  // initialize MPI:
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  world_rank=rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  walltime_0 = get_rtc();

#ifndef SVN_REV
#define SVN_REV "unknown"
#endif

// make sure 'return_moments_flag' is set correctly
  switch(evec_generation_mode)
  {
  case Constant : break;
  case Random : break;
  case WangLandau_1d :
    return_moments_flag = true;
    generator_needs_moment = true;
    break;
  case ExhaustiveIsing : break;
  case WangLandau_2d :
    return_moments_flag = true;
    generator_needs_moment = true;
    break;
  default: std::cout<<" ERROR: UNKNOWN EVEC GENERATION MODE\n"; exit(1);
  }

  if(rank==0)
  {
    std::cout<<"LSMS_3"<<std::endl;
    std::cout<<" SVN revision "<<SVN_REV<<std::endl<<std::endl;
#ifdef USE_PAPI
    std::cout<<" Using Papi counters"<<std::endl<<std::endl; 
#endif
    std::cout<<" Size of LSMS instances = "<<size_lsms<<" atoms\n";
    std::cout<<" Number of LSMS instances = "<<num_lsms<<std::endl;
    std::cout<<" LSMS Energy calculated using ";
    switch(energyCalculationMode)
    {
    case OneStepEnergy: std::cout<<"oneStepEnergy [frozen potential band energy]"<<std::endl; break;
    case MultiStepEnergy: std::cout<<"multiStepEnergy [frozen potential band energy with converged Fermi energy]"<<std::endl; break;
    case ScfEnergy: std::cout<<"scfEnergy [self-consistent total energy]"<<std::endl; break;
    default: std::cout<<"UNKNOWN ENERGY CALCULATION METHOD"<<std::endl; exit(1);
    }
    if(restrict_steps) std::cout<<" Number of gWL steps = "<<num_steps<<std::endl;
    if(restrict_time) std::cout<<" Maximum walltime = "<<max_time<<"s\n";
    std::cout<<" Processor alignment (process allocation quantization) = "<<align<<std::endl;
    switch(evec_generation_mode)
    {
    case Constant : std::cout<<" Constant moments direction along "
                             <<ev0[0]<<" "<<ev0[1]<<" "<<ev0[2]<<std::endl;
      break;
    case Random : std::cout<<" Random distribution of moments (no Wang-Landau)"<<std::endl;
      break;
    case WangLandau_1d : std::cout<<" Wang-Landau for one continuous variable (energy)"<<std::endl;
//      return_moments_flag = true;
//      generator_needs_moment = true;
      break;
    case ExhaustiveIsing : std::cout<<" Exhaustive Ising sampling"<<std::endl; break;
    case WangLandau_2d : std::cout<<" Wang-Landau for two continuous variable (energy, ";
      switch(second_dimension)
      {
      case MagneticMoment  : std::cout<<"magnitude of magnetization)"; break;
      case MagneticMomentX : std::cout<<"x component of magnetization)"; break;
      case MagneticMomentY : std::cout<<"y component of magnetization)"; break;
      case MagneticMomentZ : std::cout<<"z component of magnetization)"; break;
      }
      std::cout<<std::endl;
//      return_moments_flag = true;
//      generator_needs_moment = true;
      break;
    default: std::cout<<" ERROR: UNKNOWN EVEC GENERATION MODE\n"; exit(1);
    }
    if(step_out_flag) std::cout<<" Step output written to: "<<step_out_name<<std::endl;
    std::cout<<std::endl;

    if(step_out_flag && (evec_generation_mode==WangLandau_1d))
    {
      // step_out_flag=false;
      snprintf(wl_step_out_name,127,"wl1d_%s",step_out_name);
      wl_stepf=wl_step_out_name;
    }

    if(step_out_flag)
    {
      step_out_file.open(step_out_name);
      step_out_file<<"#";
      for(int i=0; i<argc; i++) step_out_file<<" "<<argv[i];
      step_out_file<<std::endl<<size_lsms<<std::endl;
    }
  }

  if(generator_needs_moment) return_moments_flag=true;

  if(num_lsms==1)
  {
    SHMEM_activeset local_comm;
    local_comm.rank=_my_pe();
    local_comm.size=_num_pes();
    local_comm.start_pe=0;
    local_comm.logPE_stride=0;
    LSMS lsms_calc(local_comm,i_lsms_name,"1_");
      
    if(rank==0)
    {
      std::cout<<"executing LSMS(C++) for "<<lsms_calc.numSpins()<<" atoms\n";
      std::cout<<"  LSMS version = "<<lsms_calc.version()<<std::endl;
    }

    if(energyCalculationMode==OneStepEnergy)
      std::cout<<"one step Energy = "<<lsms_calc.oneStepEnergy()<<std::endl;
    else if(energyCalculationMode==MultiStepEnergy)
      std::cout<<"multi-step Energy = "<<lsms_calc.multiStepEnergy()<<std::endl;
    else if(energyCalculationMode==ScfEnergy)
      std::cout<<"self-consistent Energy = "<<lsms_calc.scfEnergy()<<std::endl;
    else
    {
      printf("ERROR: Unknown energy calculation mode for lsms_calc in wl-lsms main!\n");
     // MPI_Abort(MPI_COMM_WORLD,5);
      exit(5);
    }
  }
  else
  {
    // build the communicators
    //int color=MPI_UNDEFINED;
    //Assuming user passes a power of two while using "-align"
    int s = align;
    int comm_size=(size-align)/num_lsms;
    int world_rank;
    for(int i=0; i<num_lsms; i++)
    {
      if((world_rank>=s) && (world_rank<s+comm_size)) 
      { 
        my_group=i; 
        //color=i; 
        new_peid=world_rank-s;
        new_root=s;
      }
      lsms_rank0[i]=s;
      s+=comm_size;
    }
    if(world_rank==0){ 
      //color=num_lsms;
      new_peid=0;
      comm_size=1;
      new_root=0;
    }

    //MPI_Comm_split(MPI_COMM_WORLD, color, 0, &local_comm);
    SHMEM_activeset local_comm;
    local_comm.rank=new_peid;
    local_comm.size=comm_size;
    local_comm.start_pe=new_root;
    local_comm.logPE_stride=0;

    std::cout<<"world_rank="<<world_rank<<" -> group="<<my_group<<std::endl;

      
    snprintf(prefix,38,"Group %4d: ",my_group);

    // now we get ready to do some calculations...

    if(my_group>=0)
    {
      double energy;
      double band_energy;
      int static i_values[10];
      double static r_values[10];
      static int op;


      //MPI_Comm_rank(local_comm, &rank);
      rank = local_comm.rank;
      snprintf(prefix,38,"%d_",my_group);
      // to use the ramdisk on jaguarpf:
      // snprintf(prefix,38,"/tmp/ompi/%d_",my_group);
      LSMS lsms_calc(local_comm,i_lsms_name,prefix);
      snprintf(prefix,38,"Group %4d: ",my_group);

      if(rank==0 && my_group==0)
      {
        std::cout<<prefix<<"executing LSMS(C++) for "<<lsms_calc.numSpins()<<" atoms\n";
        std::cout<<prefix<<"  LSMS version = "<<lsms_calc.version()<<std::endl;
      }

      // wait for commands from master
      bool finished=false;
      while(!finished)
      {
        if(rank==0)
        {
          //MPI_Recv(evec,3*size_lsms,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          //op =status.MPI_TAG;
          if (lsms_rank0[0]==world_rank)
                shmem_barrier(0, lsms_rank0[0], 2, pSync1);

        }
        //MPI_Bcast(&op,1,MPI_INT,0,local_comm);
        shmem_broadcast32(&op, &op, 1, local_comm.start_pe, local_comm.start_pe, local_comm.logPE_stride, local_comm.size, pSync2); 

/* recognized opcodes:
   5: calculate energy

   recognized energy calculation modes:
   OneStepEnergy : calclulate frozen potential band energy in one step (don't converge Ef)
   use only if the Fermi energy will not change due to MC steps!
   The only method available in LSMS_1.9
   MultiStepEnergy : calculate frozen potential band energy after converging Fermi energy
   This should be the new default method. If the Fermi energy doesn't change
   multiStepEnergy only performs one step and should be equivalent to oneStepEnergy
   The tolerance for Ef convergence can be set with LSMS::setEfTol(Real).
   The default tolerance is set in the LSMS::LSMS constructor (currently 1.0e-6).
   The maximum number of steps is read from the LSMS input file 'nscf' parameter.
   ScfEnergy : this will calculate the selfconsistent total energy.
   The maximum number of steps is read from the LSMS input file 'nscf' parameter.
   NOT IMPLEMENTED YET!!!

   10: get number of sites
*/

        if(op==5)
        {
          lsms_calc.setEvec(evec);
          if(energyCalculationMode==OneStepEnergy)
            energy=lsms_calc.oneStepEnergy(&band_energy);
          else if(energyCalculationMode==MultiStepEnergy)
            band_energy=energy=lsms_calc.multiStepEnergy();
          else if(energyCalculationMode==ScfEnergy)
            energy=lsms_calc.scfEnergy(&band_energy);
          else
          {
            printf("ERROR: Unknown energy calculation mode for lsms_calc in wl-lsms main!\n");
            //MPI_Abort(MPI_COMM_WORLD,5);
            exit(5);
          }
          r_values[0]=energy;
          r_values[1]=band_energy;
          if(return_moments_flag)
          {
            lsms_calc.getMag(&r_values[R_VALUE_OFFSET]);
          }
          if(rank==0)
          {
            if(return_moments_flag)
            {
              //MPI_Send(r_values,R_VALUE_OFFSET+3*size_lsms,MPI_DOUBLE,0,1005,MPI_COMM_WORLD);
              shmem_double_put(r_values, r_values, R_VALUE_OFFSET+3*size_lsms, 0);

            } else {
              //MPI_Send(r_values,R_VALUE_OFFSET,MPI_DOUBLE,0,1005,MPI_COMM_WORLD);
              shmem_double_put(r_values, r_values, R_VALUE_OFFSET, 0);
            }
            shmem_fence();
            shmem_int_swap(&flag, world_rank, 0);

          }
              
        } else if(op==10) {
          i_values[0]=lsms_calc.numSpins();
          //MPI_Send(i_values,10,MPI_INT,0,1010,MPI_COMM_WORLD);
          shmem_int_put(i_values, i_values, 10, 0);
        } else {
          // printf("world rank %d: recieved exit\n",world_rank); 
          finished=true;
        }
      }

      shfree(evec);
      //shfree(r_values);
    }
    else if(world_rank==0)
    {
      int running;
      double **evecs;
      //double *r_values;
      //int i_values[10];
      int *init_steps;
      int total_init_steps;
      bool accepted;
        
      char *wl_inf=NULL;
      char *wl_outf=NULL;
      if(gWL_in_name) wl_inf=gWL_in_name;
      if(gWL_out_name) wl_outf=gWL_out_name;
        
      EvecGenerator *generator;

/*
      // get number of spins from first LSMS instance
      // temp r_values:
      r_values=(double *)malloc(sizeof(double)*10);
      MPI_Send(r_values,1,MPI_DOUBLE, lsms_rank0[0], 10, MPI_COMM_WORLD);
      free(r_values);
      MPI_Recv(i_values,10,MPI_INT,lsms_rank0[0],1010,MPI_COMM_WORLD,&status);
      if(i_values[0]!=size_lsms)
      {
        printf("Size specified for Wang-Landau and in LSMS input file don't match!\n");
        size_lsms=i_values[0];
      }
*/

      evecs=(double **)shmalloc(sizeof(double *)*num_lsms);
      init_steps=(int *)shmalloc(sizeof(int)*num_lsms);
      for(int i=0; i<num_lsms; i++)
      {
        evecs[i]=(double *)shmalloc(sizeof(double)*3*size_lsms);
        init_steps[i]=initial_steps;
      }
      total_init_steps=num_lsms*initial_steps;
        

      // Initialize the correct evec generator
      switch(evec_generation_mode)
      {
      case Random :  generator = new RandomEvecGenerator(size_lsms);
        break;
      case Constant: generator = new ConstantEvecGenerator(size_lsms, ev0, num_lsms);
        break;
     //case WangLandau_1d : generator = new WL1dEvecGenerator<std::mt19937>(size_lsms, num_lsms,
     //                                                                      evecs, wl_inf, wl_outf, wl_stepf);
     case WangLandau_1d : generator = new WL1dEvecGenerator<boost::mt19937>(size_lsms, num_lsms,
                                                                           evecs, wl_inf, wl_outf, wl_stepf);
        break;
      case ExhaustiveIsing : generator = new ExhaustiveIsing1dEvecGenerator(size_lsms, num_lsms,
                                                                            evecs, wl_inf, wl_outf);
        break;
      //case WangLandau_2d : generator = new WL2dEvecGenerator<std::mt19937>(size_lsms, num_lsms,
      //                                                                     evecs, wl_inf, wl_outf, wl_stepf);
      case WangLandau_2d : generator = new WL2dEvecGenerator<boost::mt19937>(size_lsms, num_lsms,
                                                                           evecs, wl_inf, wl_outf, wl_stepf);
        break;
      default: std::cerr<<"The code should never arrive here: UNKNOWN EVEC GENERATION MODE\n";
        exit(1);
      }

      for(int i=0; i<num_lsms; i++)
      {
        generator->initializeEvec(i,evecs[i]);
      }
      std::cout<<"This is the master node\n";
      // issue initial commands to all LSMS instances
      running=0;
      bool more_work=true;
      if(total_init_steps>0)
      {
        for(int i=0; i<num_lsms; i++)
        {
          std::cout<<"starting initial calculation in group "<<i<<std::endl;
          //MPI_Send(evecs[i], 3*size_lsms, MPI_DOUBLE, lsms_rank0[i], 5, MPI_COMM_WORLD);
          shmem_double_put(evec, evecs[i], 3*size_lsms, lsms_rank0[i]);
          shmem_int_p(&op, 5, lsms_rank0[i]);
          shmem_fence();


          num_steps--; running++; stepCount++;
          if(restrict_steps) std::cout<<"      "<<num_steps<<" steps remaining\n";
        }
        shmem_barrier(0, lsms_rank0[0], 2, pSync1);
        // first deal with the initial steps:
        while(running>0)
        {
          //if(return_moments_flag)
          //  MPI_Recv(r_values,R_VALUE_OFFSET+3*size_lsms,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          //else
          //  MPI_Recv(r_values,R_VALUE_OFFSET,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          
          shmem_int_wait(&flag,-1);

          running--;
          // std::cout<<"received energy E_tot ="<<r_values[0]<<std::endl;
          // std::cout<<"    band energy E_band="<<r_values[1]<<std::endl;
          if(total_init_steps>0)
          {
            //int r_group=(status.MPI_SOURCE-align)/comm_size;
            int r_group=(flag-align)/comm_size;
            std::cout<<"starting additional calculation in group "<<r_group<<std::endl;

            if(init_steps[r_group]>0)
            {
              more_work = !(generator->generateUnsampledEvec(r_group,evecs[r_group],r_values[energyIndex]));
              init_steps[r_group]--; total_init_steps--;
            }
                
            //MPI_Send(evecs[r_group], 3*size_lsms, MPI_DOUBLE, lsms_rank0[r_group], 5, MPI_COMM_WORLD);
            shmem_double_put(r_values, evecs[r_group],  3*size_lsms, lsms_rank0[r_group]); //TODO check this
            shmem_fence();
                
            num_steps--; running++; stepCount++;
            if(restrict_steps && num_steps<=0) more_work=false;
            if(restrict_steps) std::cout<<"      "<<num_steps<<" steps remaining\n";
            walltime = get_rtc() - walltime_0;
            if(restrict_time && walltime>=max_time) more_work=false;
            if(restrict_time) std::cout<<"      "<<max_time-walltime<<" seconds remaining\n";
          }
              
        }
      }
      more_work=true;
      running=0;
      for(int i=0; i<num_lsms; i++)
      {
        std::cout<<"starting main calculation in group "<<i<<std::endl;
        //MPI_Send(evecs[i], 3*size_lsms, MPI_DOUBLE, lsms_rank0[i], 5, MPI_COMM_WORLD);
        shmem_double_put(evec, evecs[i], 3*size_lsms, lsms_rank0[i]);
        shmem_int_p(&op, 5, lsms_rank0[i]);
        shmem_fence();
        num_steps--; running++; stepCount++;
        if(restrict_steps) std::cout<<"      "<<num_steps<<" steps remaining\n";
      }
      shmem_barrier(0, lsms_rank0[0], 2, pSync1);
        
      generator->startSampling();
      // wait for results and issue new commands or wind down
      while(running>0)
      {
        //MPI_Recv(r_values,R_VALUE_OFFSET+3*size_lsms,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        shmem_int_wait(&flag,-1);

        running--;
        std::cout<<"received energy E_tot ="<<r_values[0]<<std::endl;
        std::cout<<"    band energy E_band="<<r_values[1]<<std::endl;
        // printf("from status.MPI_SOURCE=%d\n",status.MPI_SOURCE);
        energy_accumulator+=r_values[0]; energies_accumulated++;
        if(more_work)
        {
          int r_group=(status.MPI_SOURCE-align)/comm_size;
          std::cout<<"starting additional calculation in group "<<r_group<<std::endl;
              
          if(generator_needs_moment)
          {
            double m0,m1,m2;
            m0=0.0; m1=0.0; m2=0.0;
            for(int i=0; i<3*size_lsms; i+=3)
            {
              m0+=r_values[R_VALUE_OFFSET+i];
              m1+=r_values[R_VALUE_OFFSET+i+1];
              m2+=r_values[R_VALUE_OFFSET+i+2];
            }
            switch(second_dimension)
            {
            case  MagneticMoment : magnetization=std::sqrt(m0*m0+m1*m1+m2*m2); break;
            case  MagneticMomentX : magnetization=m0; break;
            case  MagneticMomentY : magnetization=m1; break;
            case  MagneticMomentZ : magnetization=m2; break;
            }
            if(generator->generateEvec(r_group,evecs[r_group],r_values[energyIndex],magnetization, &accepted))
              more_work=false;
          } else {
            if(generator->generateEvec(r_group,evecs[r_group],r_values[energyIndex], &accepted)) more_work=false;
          }

          //MPI_Send(evecs[r_group], 3*size_lsms, MPI_DOUBLE, lsms_rank0[r_group], 5, MPI_COMM_WORLD);
          shmem_double_put(r_values, evecs[r_group],  3*size_lsms, lsms_rank0[r_group]); //TODO check this
          shmem_fence();

          num_steps--; running++; stepCount++;
          if(restrict_steps && num_steps<=0) more_work=false;
          if(restrict_steps) std::cout<<"      "<<num_steps<<" steps remaining\n";
          walltime = get_rtc() - walltime_0;
          if(restrict_time && walltime>=max_time) more_work=false;
          if(restrict_time) std::cout<<"      "<<max_time-walltime<<" seconds remaining\n";
        }
        else
        {
          // send an exit message to this instance of LSMS
          int r_group=(status.MPI_SOURCE-align)/comm_size;

          MPI_Send(evecs[r_group], 3*size_lsms, MPI_DOUBLE, lsms_rank0[r_group], 2, MPI_COMM_WORLD);
        }

        if(step_out_flag && accepted)
        {
          step_out_file<<"# iteration "<<energies_accumulated<<std::endl;
          step_out_file.precision(15);
          step_out_file<<energies_accumulated<<std::endl;
          step_out_file<<r_values[0]<<"  "<<r_values[1]<<std::endl;
          for(int j=0; j<3*size_lsms; j+=3)
          {
            step_out_file<<r_values[j+R_VALUE_OFFSET]<<"  "<<r_values[j+R_VALUE_OFFSET+1]
                         <<"  "<<r_values[j+R_VALUE_OFFSET+2]<<std::endl;
          }
        }
        // write restart file every restartWriteFrequency seconds
        if(walltime>nextWriteTime)
        {
          generator->writeState("WLrestart.jsn");
          nextWriteTime+=restartWriteFrequency;
        }

      }
      generator->writeState("WLrestart.jsn");
/*
  if(evec_generation_mode==WangLandau_1d)
  (static_cast<WL1dEvecGenerator<std::mt19937> *>(generator))->writeState("WLrestart.state");
  if(evec_generation_mode==ExhaustiveIsing)
  (static_cast<ExhaustiveIsing1dEvecGenerator *>(generator))->writeState("WLrestart.state");
*/
      for(int i=0; i<num_lsms; i++) free(evecs[i]);
      shfree(evecs);
      //shfree(r_values);
    }
  }

  if(world_rank==0)
  {
    if(step_out_flag)
    {
      step_out_file<<"# end\n-1\n"
                   <<energy_accumulator/double(energies_accumulated)<<std::endl;
      step_out_file.close();
    }
    std::cout<<"Finished all scheduled calculations. Freeing resources.\n";
    std::cout<<"Energy mean = "<<energy_accumulator/double(energies_accumulated)<<"Ry\n";
  }


  if(num_lsms>1)
  {
    // make sure averyone arrives here:
    MPI_Bcast(stupid,37,MPI_CHAR,0,MPI_COMM_WORLD);

    if(world_rank==0)
    {
      MPI_Comm_free(&local_comm);
    }
    else if(my_group>=0)
    {
      MPI_Comm_free(&local_comm);
    }
  }



  if(world_rank==0)
  {
    double walltime = get_rtc() - walltime_0;
    std::cout<<" WL-LSMS finished in "<<walltime<<" seconds.\n";
    std::cout<<" Monte-Carlo steps / walltime = "
             <<double(stepCount)/walltime<<"/sec\n";
  }

#ifdef USE_PAPI
  PAPI_stop_counters(papi_values,hw_counters);
  papi_values[hw_counters  ] = PAPI_get_real_cyc()-papi_real_cyc_0;
  papi_values[hw_counters+1] = PAPI_get_real_usec()-papi_real_usec_0;
  papi_values[hw_counters+2] = PAPI_get_virt_cyc()-papi_virt_cyc_0;
  papi_values[hw_counters+3] = PAPI_get_virt_usec()-papi_virt_usec_0;
  long long accumulated_counters[NUM_PAPI_EVENTS+4];
/*
  for(int i=0; i<hw_counters; i++)
  {
  printline(ttos(papi_event_name[i])+" = "+ttos(papi_values[i]),
  std::cout,parameters.myrankWorld);
  }
  printline("PAPI real cycles : "+ttos(papi_values[hw_counters]),
  std::cout,parameters.myrankWorld);
  printline("PAPI real usecs : "+ttos(papi_values[hw_counters+1]),
  std::cout,parameters.myrankWorld);
  printline("PAPI user cycles : "+ttos(papi_values[hw_counters+2]),
  std::cout,parameters.myrankWorld);
  printline("PAPI user usecs : "+ttos(papi_values[hw_counters+3]),
  std::cout,parameters.myrankWorld);
*/
  
  //MPI_Reduce(papi_values,accumulated_counters,hw_counters+4,
  //           MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);

  shmem_long_sum_to_all(accumulated_counters, papi_values, hw_counters+4,
      comm.pestart, comm.logPE_stride, comm.size, pWrk_i, pSync2);



  if(world_rank==0)
  {
    for(int i=0; i<hw_counters; i++)
    {
      std::cout<<"Accumulated: "<<(papi_event_name[i])<<" = "<<(accumulated_counters[i])<<"\n";
    }
    std::cout<<"PAPI accumulated real cycles : "<<(accumulated_counters[hw_counters])<<"\n";
    std::cout<<"PAPI accumulated user cycles : "<<(accumulated_counters[hw_counters+2])<<"\n";
    double gflops_papi = ((double)accumulated_counters[1])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_double = ((double)accumulated_counters[2])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_single = ((double)accumulated_counters[3])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gips = ((double)accumulated_counters[0])/(1000.0*(double)papi_values[hw_counters+1]);
    std::cout<<"PAPI_FP_OPS real GFLOP/s : "<<(gflops_papi)<<"\n";
    std::cout<<"PAPI hw double real GFLOP/s : "<<(gflops_hw_double)<<"\n";
    std::cout<<"PAPI hw single real GFLOP/s : "<<(gflops_hw_single)<<"\n";
    std::cout<<"PAPI real GINST/s : "<<(gips)<<"\n";
  }
#endif


  //MPI_Finalize();
  return 0;
}
