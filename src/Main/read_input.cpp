#include <stdlib.h>

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "PhysicalConstants.hpp"

#include "SystemParameters.hpp"
#include "mixing.hpp"
#include "LuaInterface/LuaSupport.hpp"

void repeatBasisCell(LSMSSystemParameters &lsms, CrystalParameters &crystal, int nx, int ny, int nz, int unique)
{
  int numBasis=crystal.num_atoms;
  int numSites=numBasis*nx*ny*nz;

  Matrix<Real> basis,basis_evecs;
  std::vector<int> basis_type;
  basis=crystal.position;
  basis_evecs=crystal.evecs;
  basis_type=crystal.type;

  crystal.position.resize(3,numSites);
  crystal.evecs.resize(3,numSites);
  crystal.type.resize(numSites);

  if(unique)
  {
    crystal.types.resize(numSites);
    crystal.num_types=numSites;
    for(int i=0; i<numBasis; i++) if(crystal.types[i].pot_in_idx<0) crystal.types[i].pot_in_idx=i;
  }

  int i=0;
  for(int ix=0; ix<nx; ix++)
    for(int iy=0; iy<ny; iy++)
      for(int iz=0; iz<nz; iz++)
      {
        for(int ib=0; ib<numBasis; ib++)
        {
          crystal.position(0,i)=basis(0,ib)+ix*crystal.bravais(0,0)+iy*crystal.bravais(0,1)+iz*crystal.bravais(0,2);
          crystal.position(1,i)=basis(1,ib)+ix*crystal.bravais(1,0)+iy*crystal.bravais(1,1)+iz*crystal.bravais(1,2);
          crystal.position(2,i)=basis(2,ib)+ix*crystal.bravais(2,0)+iy*crystal.bravais(2,1)+iz*crystal.bravais(2,2);
          if(unique)
          {
            crystal.type[i]=i;
            if(i>=numBasis)
            {
              crystal.types[i]=crystal.types[ib];
            }
          } else {
            crystal.type[i]=basis_type[ib];
          }
          crystal.evecs(0,i)=basis_evecs(0,ib);
          crystal.evecs(1,i)=basis_evecs(1,ib);
          crystal.evecs(2,i)=basis_evecs(2,ib);
          i++;
        }
      }
  crystal.num_atoms=numSites;
  crystal.bravais(0,0)*=nx;
  crystal.bravais(1,0)*=nx;
  crystal.bravais(2,0)*=nx;
  crystal.bravais(0,1)*=ny;
  crystal.bravais(1,1)*=ny;
  crystal.bravais(2,1)*=ny;
  crystal.bravais(0,2)*=nz;
  crystal.bravais(1,2)*=nz;
  crystal.bravais(2,2)*=nz;

  lsms.num_atoms=numSites;
}

int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix)
{

// c     read in the structure  identity.................................
//       read(10,'(a)') systemid
//                      CHARACTER*50

  luaGetStrN(L,"systemid",lsms.systemid,50);
  snprintf(lsms.potential_file_in,128,"v_%s",lsms.systemid);
  snprintf(lsms.potential_file_out,128,"w_%s",lsms.systemid);
  luaGetStrN(L,"potential_file_in",lsms.potential_file_in,128);
  luaGetStrN(L,"potential_file_out",lsms.potential_file_out,128);
  lsms.pot_in_type=0; // HDF5
  luaGetInteger(L,"pot_in_type",&lsms.pot_in_type);
  lsms.pot_out_type=-1; // don't write potential
  luaGetInteger(L,"pot_out_type",&lsms.pot_out_type);

// c     read in the standard output target switch.......................
//       read(10,'(a)') output_to_screen  

// c     ================================================================
// c     read in subroutine stop level
//       read(10,'(a)') istop
// c     write(6,'(a)') istop

  char ctmp[32]; strncpy(ctmp,"main",32);
  luaGetStrN(L,"istop",ctmp,32); lsms.global.setIstop(ctmp);

// c     ================================================================
// c     read in print level for a particular node and rest of the nodes.
//       read(10,*    ) node_print,print_instr,nprint
// c     write(6,'(3i5)') node_print,print_instr,nprint
  luaGetInteger(L,"print_node",&lsms.global.print_node);
  luaGetInteger(L,"default_iprint",&lsms.global.default_iprint);
  luaGetInteger(L,"iprint",&lsms.global.iprint);
#ifdef _OPENMP
  lsms.global.GPUThreads=std::min(16,omp_get_max_threads());
#else
  lsms.global.GPUThreads=1;
#endif
  luaGetInteger(L,"gpu_threads",&lsms.global.GPUThreads);
// c     ================================================================
// c     read in the number of atoms in the system.......................
//       read(10,*    ) num_atoms
  lsms.num_atoms=0;
  luaGetInteger(L,"num_atoms", &lsms.num_atoms);
//       if(num_atoms.lt.1 .or. num_atoms.gt.max_atoms) then
//          write(6,'(/,'' RDIN_ATOM_FT::'',
//      >               '' num_atoms exceeds the upper limit'')')
//          write(6,'(  ''                num_atoms:'',i5)')num_atoms
//          write(6,'(  ''                max_atoms:'',i5)')max_atoms
//          call fstop(sname)
//       endif

  lsms.nrelc=lsms.nrelv=0;
  luaGetInteger(L,"nrelc", &lsms.nrelc);
  luaGetInteger(L,"nrelv", &lsms.nrelv);
  lsms.clight=cphot*std::pow(10.0,lsms.nrelv);
// c     ================================================================
// c     read in indices that cotrol relativity (core & valence) and
// c     muffin tin/ASA switch
//       read(10,'(a)') text
//       read(10,*,err=1001) nrelc,nrelv,mtasa,rmt0,vshift,vshiftl
//       if(nrelc.ne.0) then
//          nrelc=10
//       endif
//       if(nrelv.ne.0) then
//          nrelv=10
//       endif
//       if( mtasa .ge. 3 ) then
  lsms.mtasa = 0;
  luaGetInteger(L,"mtasa",&lsms.mtasa);
  lsms.fixRMT = 0;
  luaGetInteger(L,"fixRMT",&lsms.fixRMT);
//       else if( mtasa .lt. -3 ) then
//            write(6,'(/,'' RDIN_ATOM_FT:: mtasa andor rmt0'',i5,d13.4)')
//      >     mtasa,rmt0
// 	   call fstop(sname)
//       endif

// c     read spin polarization index....................................
//       read(10,*    ) nspin,i_vdif,iexch
 
/* nspin = 1 : non spin polarized
           2 : spin polarized, collinear
           3 : spin polarized, non-collinear
           4 : fully relativistic
*/
  lsms.nspin=2;
  lsms.nrel_rel=false;
  luaGetInteger(L,"nspin",&lsms.nspin);
  if(lsms.nspin>3)
  {
    lsms.nrel_rel=true;
    lsms.nspin=3;
  }

  if(lsms.nspin>1) lsms.n_spin_pola=2; else lsms.n_spin_pola=1;
  if(lsms.nspin>2) lsms.n_spin_cant=2; else lsms.n_spin_cant=1;

//       if (nspin.lt.0 .or.(nspin.le.3.and. nspin.gt.ipspin+1)) then
//          write(6,'('' RDIN_ATOM_FT:: Wrong input for nspin'')')
//          call fstop(sname)
//       endif

// Read Bravais Vectors:
  for(int i=0; i<3; i++)
  {
    luaGetPositionInTable(L,"bravais",i+1);
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.bravais(j,i));
    lua_pop(L,1);
  }

// Read atom positions and evecs
  crystal.resize(lsms.num_atoms);
  crystal.resizeTypes(lsms.num_atoms);
  crystal.num_atoms=lsms.num_atoms;
  crystal.num_types=0;
  for(int i=0; i<crystal.num_atoms; i++)
  {
    int t; // the type to be assigned

    luaGetPositionInTable(L,"site",i+1);
    luaGetFieldFromStack(L,"pos");
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.position(j,i));
    lua_pop(L,1);

    if(!luaGetIntegerFieldFromStack(L,"type",&t) || (t-1)==i)
    {
      luaGetStrNFromStack(L,"atom",crystal.types[crystal.num_types].name,4);
      luaGetIntegerFieldFromStack(L,"lmax",&crystal.types[crystal.num_types].lmax);
      luaGetIntegerFieldFromStack(L,"Z",&crystal.types[crystal.num_types].Z);
      luaGetIntegerFieldFromStack(L,"Zc",&crystal.types[crystal.num_types].Zc);
      luaGetIntegerFieldFromStack(L,"Zs",&crystal.types[crystal.num_types].Zs);
      luaGetIntegerFieldFromStack(L,"Zv",&crystal.types[crystal.num_types].Zv);
      luaGetRealFieldFromStack(L,"rLIZ",&crystal.types[crystal.num_types].rLIZ);
printf("lmax=%d,Z=%d,Zc=%d,Zs=%d,Zv=%d,lmax=%f,",crystal.types[crystal.num_types].lmax,crystal.types[crystal.num_types].Z,crystal.types[crystal.num_types].Zc,crystal.types[crystal.num_types].Zs,crystal.types[crystal.num_types].Zv,crystal.types[crystal.num_types].rLIZ);
      luaGetFieldFromStack(L,"rsteps");
      for(int j=0; j<4; j++) luaGetRealPositionFromStack(L,j+1,&crystal.types[crystal.num_types].rsteps[j]);
      lua_pop(L,1);
      crystal.types[crystal.num_types].first_instance=i;
      crystal.types[crystal.num_types].number_of_instances=1;
      crystal.type[i]=crystal.num_types;
      crystal.num_types++;
    } else if(t<i && t>=0) {
      crystal.type[i]=crystal.type[t-1];
      crystal.types[crystal.type[i]].number_of_instances++;
    } else {
      fprintf(stderr,"Illegal type reference for atom %d : %d!\n",i+1,t);
      exit(1);
    }

    luaGetFieldFromStack(L,"evec");
    for(int j=0; j<3; j++) luaGetRealPositionFromStack(L,j+1,&crystal.evecs(j,i));
    lua_pop(L,1);
    lua_pop(L,1);
  }

  int xRepeat=1;
  luaGetInteger(L,"xRepeat",&xRepeat);
  int yRepeat=1;
  luaGetInteger(L,"yRepeat",&yRepeat);
  int zRepeat=1;
  luaGetInteger(L,"zRepeat",&zRepeat);
  int makeTypesUnique=0;
  luaGetInteger(L,"makeTypesUnique",&makeTypesUnique);

  repeatBasisCell(lsms, crystal, xRepeat, yRepeat, zRepeat,makeTypesUnique);

// c     ================================================================
// c     read in a title to identify the system .........................
//       read(10,'(a)') system_title
  luaGetStrN(L,"system_title",lsms.title,80);
// c     ================================================================
// c     Read number of Gaussian points for r and theta integrations.....
//       read(10,*    ) ngaussr,ngaussq
// c     ================================================================
// c     Read in name of the atom........................................
// c     Read in the atom position vector................................
// c     Read in cut off radius for the LIZ of the atom..................
// c     Read in radius steps for lmax, lmax-1, lmax-2, lmax-3 ..........
// c     ================================================================
// c
// c     ================================================================
// c     read in names of info_table & info_evec files:..................
// c     ================================================================
//       read(10,'(a)') text
//       read(10,'(2a30)')info_table,info_evec

// c     ================================================================
// c     read in parameters that control energy inregration:.............
// c     ================================================================
// c     igrid   : specifies energy contour  1=slow......................
// c     igrid   : specifies energy contour  2=gaussian..................
// c     igrid   : specifies energy contour  3=Don's Fermi function poles
// c
// c     for igrid =1...[Zero Temperature Calculations Only].............
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c     etop    : top of contour: real axis [usually reset to chempot]..
// c     eitop   : top    of contour on imaginary axis...................
// c     eibot   : bottom of contour on imaginary axis...................
// c     npts    : number of energy points per 0.1 ry....................
// c     kelvin  : not used..............................................
// c
// c     for igrid =2...[Zero Temperature Calculations Only].............
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c     etop    : top of contour: real axis [usually reset to chempot]..
// c     eitop   : not used..............................................
// c     eibot   : not used..............................................
// c     npts    : number of Gaussian distributed energy points..........
// c     kelvin  : not used..............................................
// c
// c     for igrid =3...[Finite Temperature Calculations Only]...........
// c     ebot    : bottom of contour: real axis [may be mod. by semcor]..
// c                                            [then reset in congauss].
// c     etop    : top of contour on real axis [usually reset to chempot]
// c     eitop   : not used..............................................
// c     eibot   : not used..............................................
// c     npts    : not used..............................................
// c     kelvin  : temperature in kelvin..................................
// c     nument  : # of gaussian points on elliptical contour for Entropy
// c     ================================================================

  luaGetIntegerFieldInTable(L,"energyContour","grid",&lsms.energyContour.grid);
  luaGetRealFieldInTable(L,"energyContour","ebot",&lsms.energyContour.ebot);
  luaGetRealFieldInTable(L,"energyContour","etop",&lsms.energyContour.etop);
  luaGetRealFieldInTable(L,"energyContour","eibot",&lsms.energyContour.eibot);
  luaGetRealFieldInTable(L,"energyContour","eitop",&lsms.energyContour.eitop);
  luaGetIntegerFieldInTable(L,"energyContour","npts",&lsms.energyContour.npts);
  lsms.energyContour.maxGroupSize=50;
  luaGetIntegerFieldInTable(L,"energyContour","maxGroupSize",&lsms.energyContour.maxGroupSize);

// c
// c     ================================================================
// c     read in controls for performing SCF calculation:................
// c     ================================================================
// c     nscf        : maximum number of scf iterations requested........
  lsms.nscf = 1;
  luaGetInteger(L, "nscf", &lsms.nscf);

// read the frequency of writing the potential during an scf calculation (writeSteps)
  lsms.writeSteps=30000;
  luaGetInteger(L, "writeSteps", &lsms.writeSteps);

// c     alpdv       : mixing parameter for chg. den. or potential.......
// c     alpma       : mixing parameter for moment density...............
// c     alpev       : mixing parameter for moment orientation...........
// c     mix_quant   : mixing charge density[potential] => 0[1]..........
// c     mix_algor   : mixing simple[DGAnderson] => 0[1].................
// lsms.mixing = 4*mix_algor+mix_quant; -1: no mixing
  lsms.mixing = -1;

  char quantity[80], algorithm[80];
  int numberOfMixQuantities = 0;

  // No mixing is set by default
  for (int i = 0; i < 5; i++)
  {
    if (i) mix.quantity[i] = false;
    else mix.quantity[i] = true;
    mix.algorithm[i] = MixingParameters::noAlgorithm;
    mix.mixingParameter[i] = 0.0;
  }

  luaGetInteger(L, "numberOfMixQuantities", &numberOfMixQuantities);

  for (int i = 0; i < numberOfMixQuantities; i++)
  {
    luaGetPositionInTable(L, "mixing", i+1);
    luaGetStrNFromStack(L, "quantity", quantity, 50);

    mix.quantity[0] = false;
    int quantityIdx = 0;
    if (strcmp("charge", quantity) == 0) 
      quantityIdx = MixingParameters::charge;
    else if (strcmp("potential", quantity) == 0) 
      quantityIdx = MixingParameters::potential;
    else if (strcmp("moment_magnitude", quantity) == 0) 
      quantityIdx = MixingParameters::moment_magnitude;
    else if (strcmp("moment_direction", quantity) == 0)
      quantityIdx = MixingParameters::moment_direction;

    mix.quantity[quantityIdx] = true;

    luaGetStrNFromStack(L, "algorithm", algorithm, 50);

    if (strcmp("simple", algorithm) == 0)
      mix.algorithm[quantityIdx] = MixingParameters::simple;
    else if (strcmp("broyden", algorithm) == 0)
      mix.algorithm[quantityIdx] = MixingParameters::broyden;

    luaGetRealFieldFromStack(L, "mixingParameter", &mix.mixingParameter[quantityIdx]);

    lua_pop(L,2);

   }

// c     iharris = 0 : do not calculate harris energy....................
// c     iharris = 1 : calculate harris energy using updated chem. potl..
// c     iharris >=2 : calculate harris energy at fixed chem. potl.......
// c     i_potwrite  : the number of iterations between potential writes.
// c     movie   = 0 : no movie data will be written.....................
// c             = 1 : movie data will be written........................
// c     ctq         : coefficient of torque ............................
// c     ================================================================
//       read(10,'(a)') text
//       read(10, *   ) nscf,alpdv,alpma,alpev,mix_quant,mix_algor,
//      >               iharris,i_potwrite,movie
// c     ================================================================
// c     check consistencey of these parameters..........................
// c     ================================================================
//       if(i_potwrite.gt.nscf) then
//          i_potwrite=nscf
//       endif
//       if(i_potwrite.lt.-nscf) then
//          i_potwrite=-nscf
//       end if
//       if(mix_quant.ne.0 .and. mix_quant.ne.1) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for mix_quant ='',
//      >             1i3)')mix_quant
//          call fstop(sname)
//       else if(mix_algor.gt.2) then 
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for mix_algor ='',
//      >             1i3)')mix_algor
//          call fstop(sname)
//       endif
// c     ================================================================
// c     if calculating the harris energy make sure that mix_quant switch
// c     is set mix the charge density...................................
// c     ================================================================
//       if(iharris.ne.0) then
//          mix_quant=0
//       endif
// c
// c     ================================================================
// c     read in quantities that control Spin Dynamics :.................
// c     ================================================================
// c     nstep        : number of time steps [for SCF only : nstep=1 ]...
// c     tstep        : time step........................................
// c     etol         : over-ride scf convergence tolerence : energy ....
// c     ptol         : over-ride scf convergence tolerence : pressure ..
// c     eftol        : over-ride scf convergence tolerence : Fermi engy.
// c     rmstol       : over-ride scf convergence tolerence : rmstol ....
// c     ----------------------------------------------------------------
// c     etol,ptol,eftol, & rmstol can be relaxed for SD.................
// c     ================================================================
//       read(10,'(a)') text
// c     write(6,'(a70)')text
//       read(10,*) ntstep,tstep,etol,ptol,eftol,rmstol
// c     write(6,'('' RDIN_ATOM_FT: ntstep,tstep,etol,ptol,eftol,rmstol'',
// c    >i4,f10.4,3x,4d10.5)') ntstep,tstep,etol,ptol,eftol,rmstol
// c
// c     ================================================================
// c     read controls for calculation of torque & exchange interactions.
// c     ================================================================
// c     ctq        : 
// c     j_ij       : 
// c     ================================================================
//       read(10,'(a)') text
//       read(10, *   ) ctq,j_ij
// c     ================================================================
// c     check consistencey of these parameters..........................
// c     ================================================================
//       if(ctq.lt.0.0d0) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for ctq < 0.0'',
//      >             1f10.5)')ctq
//          call fstop(sname)
//       else if(j_ij.ne.0 .and. j_ij.ne.1) then
//          write(6,'('' RDIN_ATOM_FT::'',
//      >             '' Incorrect input data for j_ij <> 0, and <> 1'',
//      >             1i5)')j_ij
//          call fstop(sname)
//       endif
  return 0;
}
