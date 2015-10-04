#include "mixing.hpp"


class NoMixing : public Mixing {
public:
  void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a) {}
  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
  void updatePotential(LSMSSystemParameters &lsms, AtomData &a) {}
  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
  void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
};



class FrozenPotential : public Mixing {
 
public:

  void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  {
    a.rhotot = a.rhoNew;
    a.xvalws[0] = a.xvalwsNew[0];
    a.xvalws[1] = a.xvalwsNew[1];
  }

  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i = 0; i < as.size(); i++) {
      as[i].rhotot = as[i].rhoNew;
      as[i].xvalws[0] = as[i].xvalwsNew[0];
      as[i].xvalws[1] = as[i].xvalwsNew[1];
    }
  }

  void updatePotential(LSMSSystemParameters &lsms, AtomData &a) {}

  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

  void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};



class SimpleChargeDensityMixing : public Mixing {

public:
  Real alpha;

  SimpleChargeDensityMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  { 
    simpleMixing( &a.rhotot(0,0), &a.rhoNew(0,0), a.rhotot.size(), alpha);
    simpleMixing( &a.xvalws[0], &a.xvalwsNew[0], 2, alpha);
  }

  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++)
    { 
      simpleMixing( &as[i].rhotot(0,0), &as[i].rhoNew(0,0), as[i].rhotot.size(), alpha);
      simpleMixing( &as[i].xvalws[0], &as[i].xvalwsNew[0], 2, alpha);
    }
  }

  void updatePotential(LSMSSystemParameters &lsms, AtomData &a)
  { 
    a.vr = a.vrNew;
  }

  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++)
      as[i].vr = as[i].vrNew;
  }

  void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};



class SimplePotentialMixing : public Mixing {

public:
  Real alpha;
  
  SimplePotentialMixing(Real _alpha) : alpha(_alpha) {}
  
  void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  { 
    a.rhotot = a.rhoNew;
    a.xvalws[0] = a.xvalwsNew[0];
    a.xvalws[1] = a.xvalwsNew[1];
  }
  
  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++) {
      as[i].rhotot = as[i].rhoNew;
      as[i].xvalws[0] = as[i].xvalwsNew[0];
      as[i].xvalws[1] = as[i].xvalwsNew[1];
    }
  }   
  
  void updatePotential(LSMSSystemParameters &lsms, AtomData &a)
  { 
    simpleMixing( &a.vr(0,0), &a.vrNew(0,0), a.vr.size(), alpha);
  }
  
  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++)
      simpleMixing( &as[i].vr(0,0), &as[i].vrNew(0,0), as[i].vr.size(), alpha);
  }
  
  void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};



class EfMixing : public Mixing {
  Real efOld, alpha;

public:

  EfMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  {
    lsms.chempot = alpha * lsms.chempot + (1.0-alpha) * efOld;
  }

  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    lsms.chempot = alpha * lsms.chempot + (1.0-alpha) * efOld;
  }

  void updatePotential(LSMSSystemParameters &lsms, AtomData &a)
  {
    //a.vr = a.vrNew;
  }

  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    //for(int i=0; i<as.size(); i++)
    //  as[i].vr = as[i].vrNew;
  }

  void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    efOld = lsms.chempot;
  }

};



void setupMixing(MixingParameters &mix, Mixing* &mixing)
{

  printf("\n");
  mixing = NULL;

  // frozen potential by default
  if (!mix.quantity[MixingParameters::no_mixing] &&
      !mix.quantity[MixingParameters::charge] && 
      !mix.quantity[MixingParameters::potential] && 
      !mix.quantity[MixingParameters::moment_magnitude] && 
      !mix.quantity[MixingParameters::moment_direction])
  {
    mixing = new FrozenPotential;
    printf("Mixing method     : frozen potential (default)\n");
  }
  // no mixing
  else if (mix.quantity[MixingParameters::no_mixing])
  {
    mixing = new NoMixing;
    printf("Mixing method     : no mixing\n");
  }
  // charge mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
            mix.quantity[MixingParameters::charge] && 
           !mix.quantity[MixingParameters::potential] && 
           !mix.quantity[MixingParameters::moment_magnitude] && 
           !mix.quantity[MixingParameters::moment_direction])
  {
    switch (mix.algorithm[MixingParameters::charge]) {
      case 1 :
        mixing = new SimpleChargeDensityMixing(mix.mixingParameter[MixingParameters::charge]);
        printf("Mixing method     : simple\n");
        break;
      case 2 :
        // broyden (not yet implemented)
        printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new NoMixing;
        printf("Mixing method     : no mixing\n");
    }
    printf("Mixing quantity   : charge\n");
    printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::charge]);
  }
  // potential mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
           !mix.quantity[MixingParameters::charge] && 
            mix.quantity[MixingParameters::potential] && 
           !mix.quantity[MixingParameters::moment_magnitude] && 
           !mix.quantity[MixingParameters::moment_direction])
  {
    switch (mix.algorithm[MixingParameters::potential]) {
      case 1 :
        if (mix.mixingParameter[MixingParameters::potential] == 0.0) {
          mixing = new FrozenPotential;
          printf("Mixing method     : frozen potential\n");
        }
        else {
          mixing = new SimplePotentialMixing(mix.mixingParameter[MixingParameters::potential]);
          printf("Mixing method     : simple\n");
        }
        break;
      case 2 :
        // broyden (not yet implemented)
        printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new FrozenPotential;
        printf("Mixing method     : frozen potential\n");
    }
    printf("Mixing quantity   : potential\n");
    printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::potential]);
  }
  else
  {
    printf("Type of mixing is not supported.\n");
    for (int i = 0; i < mix.numQuantities; i++) {
      printf("quantity = %5d, algorithm = %5d, mixing parameter = %6.3f\n", 
             mix.quantity[i], mix.algorithm[i], mix.mixingParameter[i]);
    }
    exit(1);
  }

}

