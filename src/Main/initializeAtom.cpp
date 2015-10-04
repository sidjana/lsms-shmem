#include <vector>
#include <algorithm>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "SingleSite/AtomData.hpp"
#include "initializeAtom.hpp"
#include "Core/CoreStates.hpp"

/* initializeAtom(AtomData &a)
   initialize an atom using the following information in AtomData

    Mesh information: xstart, rmt, jmt, jws
    Atom information: ztotss, zcorss, zsemss, zvalss
    (Note that ztotss=zcorss+zsemss+zvalss)
    lmax to limit core search
*/
class InitialAtomLevels {
public:
  Real energy;
  int n,l;
  // bool operator<()(InitialAtomLevels const &a, InitialAtomLevels const &b) { return a.energy<b.energy; }
} ;

struct compareInitialAtomLevels {
  bool operator()(InitialAtomLevels const &a, InitialAtomLevels const &b) { return a.energy<b.energy; }
};

/*
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deepst(nqn,lqn,kqn,en,rv,r,rf,h,z,c,
     >                  nitmax,tol,nws,nlast,iter,iprpts,ipdeq)

c               nqn:     principal quantum number; 
c               lqn:     orbital quantum number;
c               kqn:     kappa quantum number; 
c               en:      energy; 
c               rv:      potential in rydbergs times r;
c               r:       radial log. grid; 
c               rg:      big component; 
c               rf:      small component;
c               h:       exp. step; 
c               z:       atomic number; 
c               nitmax:  number of iterations;
c               tol:     energy tolerance; 
c               nws:     bounding sphere radius index;
c               nlast:   last tabulation point; 
c               c:       speed of light in rydbergs;
c               drg,drf: wavefunctions derivatives times r;
c               gam:     first power for small r expansion; 
c               slp:     slope at the origin;
c               dm:      h/720;
*/
extern "C"
{
  void deepst_(int *nqn, int *lqn, int *kqn, Real *en, Real *rv, Real*r,
               Real *rf, Real*h, Real *z, Real *c,
               int *nitmax, Real *tol, int *nws, int *nlast, int *iter,
               int *iprpts, int *ipdeq);
}

// approximation for the error function
Real approxErfc(Real x)
{
  const Real p=0.47047;
  const Real a1=0.3480242;
  const Real a2=-0.0958798;
  const Real a3=0.7478556;

  Real t=1.0/(1.0+p*x);
  return 1.0-((a1+(a2+a3*t)*t)*t)*std::exp(-(x*x));
}

void initializeAtom(AtomData &a)
{
  a.generateRadialMesh();

  // inititalize potential to be V(r) = -2Z/r
  // add a potential corresponding to a gaussian charge distribution
  // for the core charge with \sigma=0.2*rmt 
  // (vr stores V(r) * r)
  Real sigmaSqrt2Inv=1.0/(0.2*std::sqrt(2.0)*a.rmt);
  Real q=a.zcorss+a.zsemss;
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.vr(ir,0)=a.vr(ir,1)= -2.0*a.ztotss+q*approxErfc(a.r_mesh[ir]*sigmaSqrt2Inv);
  }
  // add homogeneous charge density inside the sphere from valence electrons
  // do make the site neutral. (i.e. vr(r_mesh.size(),*)=0)
  // q=-a.vr(a.r_mesh.size()-1,0)/a.r_mesh[a.r_mesh.size()-1];
  // for(int  ir=0; ir<a.r_mesh.size(); ir++)
  // {
  //   Real r=a.r_mesh[ir];
  //   a.vr(ir,0)+=q*r*r;
  //   a.vr(ir,1)+=q*r*r;
  // }

  // find the lowest zcorss+zsemss levels
  // maximum number of states to test: (lmax+1)*(lmax+2)/2
  // assuming kappa degeneracy, iterating over principal quantum number and l
  // we treat core and semi-core states the same for the initialization,
  // using deepst for both

  int lmaxCore=std::min(a.lmax,3); // only consider s,p,d,f electrons (we should never need g or higher)

  int numAtomLevels=((lmaxCore+1)*(lmaxCore+2))/2;
  int kappa;
  Real energy;
  Real cLight=2.0*137.0359895;
  int nitmax=100; // maximum number of iterations
  Real tol=1.0e-8; // energy tolerance
  int last = a.r_mesh.size();
  int ipdeq=5;
  int iter;
  std::vector<InitialAtomLevels> atomLevels(numAtomLevels);
  std::vector<Real> rf(a.r_mesh.size()+1);
  int nl=0;
  for(int n=1; n<=lmaxCore+1; n++)
  {
    for(int l=0; l<n; l++)
    {
      kappa=-l-1; //only use kappa -l-1 this will work for all l (including l=0)
      energy=-(a.ztotss)*(a.ztotss)/Real((n)*(n));

      printf("calculating n=%2d  l=%2d  :  Energy guess=%12.6lf\n",n,l,energy);

      deepst_(&n, &l, &kappa, &energy, &a.vr(0,0), &a.r_mesh[0], &rf[0], &a.h,
             &a.ztotss, &cLight, &nitmax, &tol, &a.jws, &last, &iter,
             &last, &ipdeq);

      atomLevels[nl].n=n; atomLevels[nl].l=l; atomLevels[nl].energy=energy;
      printf("n=%2d  l=%2d  :  Energy=%12.6lf\n",n,l,energy);
      nl++;
    }
  }
  // sort
  std::sort(atomLevels.begin(),atomLevels.end(),compareInitialAtomLevels());
//            [](myclass const & a, myclass const &b){return a.energy < b.energy;});
// fill core states:
  int coreTarget=a.zcorss+a.zsemss;
  int coreElectrons=0;
  a.numc=0;
  for(int i=0; i<atomLevels.size(); i++)
  {
    if(coreElectrons<coreTarget)
    {
      coreElectrons+=2*(2*atomLevels[i].l+1); printf("C ");
      if(atomLevels[i].l==0) a.numc+=1; else a.numc+=2;
    } else printf("V ");
    printf("n=%2d  l=%2d  :  Energy=%12.6lf\n",atomLevels[i].n,atomLevels[i].l,atomLevels[i].energy);
  }
  if(coreElectrons!=coreTarget)
  {
    printf("Warning: initializeAtom can't satisfy the core electron requirement:\n  Target: %d (%lf + %lf)\n  Actual: %d\n",
           coreTarget,a.zcorss,a.zsemss,coreElectrons);
  }
  a.resizeCore(a.numc);
  int j=0;
  for(int i=0; i<a.numc; i++)
  {
    a.ec(i,0)=a.ec(i,1)=atomLevels[j].energy;
    a.nc(i,0)=a.nc(i,1)=atomLevels[j].n;
    a.lc(i,0)=a.lc(i,1)=atomLevels[j].l;
    a.kc(i,0)=a.kc(i,1)=-atomLevels[j].l-1;
    if(atomLevels[j].l!=0)
    {
      i++;
      a.ec(i,0)=a.ec(i,1)=atomLevels[j].energy;
      a.nc(i,0)=a.nc(i,1)=atomLevels[j].n;
      a.lc(i,0)=a.lc(i,1)=atomLevels[j].l;
      a.kc(i,0)=a.kc(i,1)=atomLevels[j].l;
    }
    j++;
    // printf("%d %d\n",i,j);
  }
  //exit(1);
}


int initializeNewPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  for(int i=0; i<local.num_local; i++)
  {
    printf("Initializing potential %d.\n",i);
    snprintf(local.atom[i].header,40,"Potential Initialized by LSMS_3");
    for(int j=31;j<80;j++) local.atom[i].header[j]=' ';
    local.atom[i].resizePotential(lsms.global.iprpts);
    local.atom[i].ztotss=(Real)crystal.types[local.global_id[i]].Z;
    local.atom[i].zcorss=(Real)crystal.types[local.global_id[i]].Zc;
    local.atom[i].zsemss=(Real)crystal.types[local.global_id[i]].Zs;
    local.atom[i].zvalss=(Real)crystal.types[local.global_id[i]].Zv;
    local.atom[i].vdif=0.0;
    local.atom[i].xvalws[0]=local.atom[i].xvalws[1]=0.5*local.atom[i].zvalss;
    local.atom[i].lmax=crystal.types[local.global_id[i]].lmax;
    local.atom[i].kkrsz=(local.atom[i].lmax+1)*
                        (local.atom[i].lmax+1);

    local.atom[i].xstart=-11.1309674;
    local.atom[i].jmt=1001;
    local.atom[i].jws=1001;
    local.atom[i].nspin=2;
    local.atom[i].evec[0]=local.atom[i].evec[1]=0.0; local.atom[i].evec[2]=1.0;

    initializeAtom(local.atom[i]);
    local.atom[i].alat=local.atom[i].rmt;
    local.atom[i].efermi=0.5;
  }
  lsms.chempot=0.5;
  calculateCoreStates(comm, lsms, local);
  return 0;
}
