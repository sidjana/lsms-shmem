#include "Real.hpp"
#include "Complex.hpp"
#include <vector>
#include <cmath>
#include <cblas.h>
#include "PhysicalConstants.hpp"

#include "lapack.h"

// #include "Communication/LSMSCommunication.hpp"
#include "SingleSite/SingleSiteScattering.hpp"
// #include "MultipleScattering.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"

void solveSingleScatterers(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                           std::vector<Matrix<Real> > &vr, Complex energy,
                           std::vector<NonRelativisticSingleScattererSolution> &solution,int iie)
{
// ========================== SINGLE SCATTERER STUFF
  Complex prel=std::sqrt(energy*(1.0+energy*c2inv));
  Complex pnrel=std::sqrt(energy);


  for(int i=0; i<local.num_local; i++)
    solution[i].init(lsms,local.atom[i],&local.tmatStore(iie*local.blkSizeTmatStore,i));

  if(lsms.nrelv>0) prel=pnrel;

  if(local.atom.size()>solution.size()) solution.resize(local.atom.size());

// this is not ready for multithreading yet
// else: pragma omp parallel for
  // if(lsms.global.iprint>=0) printf("calculate single scatterer solutions.\n");
  for(int i=0; i<local.atom.size(); i++)
  {
    // if(lsms.global.iprint>=1) printf("calc single scatterer %d.%d\n",comm.rank,i);
    //printf("calc single scatterer atom no. = %d\n",i);
    calculateSingleScattererSolution(lsms,local.atom[i],vr[i],energy,prel,pnrel,solution[i]);
// calculate pmat_m (needed for tr_pxtau)
    int kkrsz=local.atom[i].kkrsz;
    int kkrszsqr=kkrsz*kkrsz;
    int info;
    int ipvt[kkrsz];
    if(lsms.n_spin_cant==2 && lsms.nrel_rel==0)
    {
      local.atom[i].pmat_m[iie].resize(kkrsz,kkrsz);
      Complex *pmat = new Complex[kkrszsqr];
      Complex *wbig = new Complex[kkrszsqr];
      Complex *pmat_m_ptr=&local.atom[i].pmat_m[iie](0,0);
      cblas_zcopy(kkrszsqr,&solution[i].tmat_l(0,0,0),1,pmat,1);
      cblas_zcopy(kkrszsqr,&solution[i].tmat_l(0,0,1),1,pmat_m_ptr,1);
      zgetrf_(&kkrsz,&kkrsz,pmat_m_ptr,&kkrsz,ipvt,&info);
      zgetri_(&kkrsz,pmat_m_ptr,&kkrsz,ipvt,wbig,&kkrszsqr,&info);
//    -------------------------------------------------------------
      zgetrf_(&kkrsz,&kkrsz,pmat,&kkrsz,ipvt,&info);
      zgetri_(&kkrsz,pmat,&kkrsz,ipvt,wbig,&kkrszsqr,&info);

      for(int j=0; j<kkrszsqr; j++) pmat_m_ptr[j]-=pmat[j];

      delete [] wbig;
      delete [] pmat;
    }
  }



// ========= END SINGLE SCATTERER STUFF ======================
}
