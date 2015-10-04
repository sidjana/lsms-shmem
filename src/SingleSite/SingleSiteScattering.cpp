#include <cblas.h>
#include "SingleSiteScattering.hpp"

extern "C"
{
void trltog_(int *, int *, Complex *, Complex *, Complex *, Complex *, Complex *);
}

void calculateSingleScattererSolution(LSMSSystemParameters &lsms, AtomData &atom,
                                      Matrix<Real> &vr,
                                      Complex energy, Complex prel, Complex pnrel,
                                      NonRelativisticSingleScattererSolution &solution)
{
  int iprpts=atom.r_mesh.size();
  solution.energy=energy;
  // Real r_sph=atom.r_mesh[atom.jws];
  // if(lsms.mtasa==0) r_sph=atom.r_mesh[atom.jmt];
  Real r_sph=atom.rInscribed;
  if(lsms.mtasa>0) r_sph=atom.rws;

  if(lsms.n_spin_pola==1) // non spin polarized
  {
    single_scatterer_nonrel_(&lsms.nrelv, &lsms.clight, &atom.lmax, &atom.kkrsz,
                             &energy,&prel,&pnrel,
                             &vr(0,0),&atom.r_mesh[0],&atom.h,&atom.jmt,&atom.jws,
                             &solution.tmat_l(0,0,0),&solution.matom(0,0),
                             &solution.zlr(0,0,0),&solution.jlr(0,0,0),
                             &r_sph,
                             &iprpts,
                             &lsms.global.iprint,lsms.global.istop,32);
    int kkrszsqr=atom.kkrsz*atom.kkrsz;
    int one=1;
    cblas_zcopy(kkrszsqr,&solution.tmat_l(0,0,0),1,&solution.tmat_g(0,0),1);
  } else {
    for(int is=0; is<lsms.n_spin_pola; is++)
    {
      single_scatterer_nonrel_(&lsms.nrelv, &lsms.clight, &atom.lmax, &atom.kkrsz,
                               &energy,&prel,&pnrel,
                               &vr(0,is),&atom.r_mesh[0],&atom.h,&atom.jmt,&atom.jws,
                               &solution.tmat_l(0,0,is),&solution.matom(0,is),
                               &solution.zlr(0,0,is),&solution.jlr(0,0,is),
                               &r_sph,
                               &iprpts,
                               &lsms.global.iprint,lsms.global.istop,32);
    }
    if(lsms.n_spin_cant>1)
    {
      trltog_(&atom.kkrsz,&atom.kkrsz,&atom.ubr[0],&atom.ubrd[0],
              &solution.tmat_l(0,0,0), &solution.tmat_l(0,0,1),&solution.tmat_g(0,0));
    } else {
      int kkrszsqr=atom.kkrsz*atom.kkrsz;
      int one=1;
      cblas_zcopy(kkrszsqr,&solution.tmat_l(0,0,0),1,&solution.tmat_g(0,0),1);
      cblas_zcopy(kkrszsqr,&solution.tmat_l(0,0,1),1,&solution.tmat_g(0,atom.kkrsz),1);
    }
  }
}

void calculateScatteringSolutions(LSMSSystemParameters &lsms, std::vector<AtomData> &atom,
                                  Complex energy, Complex prel, Complex pnrel,
                                  std::vector<NonRelativisticSingleScattererSolution> &solution)
{
  // if(atom.size()>solution.size()) solution.resize(atom.size());
  for(int i=0; i<atom.size(); i++)
    calculateSingleScattererSolution(lsms,atom[i],atom[i].vr,energy,prel,pnrel,solution[i]);
}
