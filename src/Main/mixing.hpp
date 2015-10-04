#ifndef LSMS_MIXING_H
#define LSMS_MIXING_H
#include "Real.hpp"
#include "SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
#include <vector>
#include <cmath>


struct MixingParameters {

  // Different mixing quantities and algorithms
  static const int numQuantities = 5;

  enum mixQuantity {no_mixing = 0, charge = 1, potential = 2, moment_magnitude = 3,
                    moment_direction = 4};
  enum mixAlgorithm {noAlgorithm = 0, simple = 1, broyden = 2};
  
  // These parameters specify the which quantity(ies) is (are) being mixed and which algorithm(s) to used.
  // The correspondances of the indices are specified in mixQuantity.
  // bool values:
  // 0 : quantity is not used for mixing
  // 1 : quantity is used for mixing
  bool quantity[numQuantities];
  mixAlgorithm algorithm[numQuantities];
  Real mixingParameter[numQuantities];

};


template <typename T>
void simpleMixing(T *fold, T* fnew, int n, Real alpha)
{
  if(alpha>1.0) alpha = 1.0;
  if(alpha<0.0) alpha = 0.0;
  Real beta = 1.0 - alpha;

  for(int i=0; i<n; i++)
    fold[i] = alpha * fnew[i] + beta * fold[i];
}


class Mixing {
public:
  virtual void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
  virtual void updatePotential(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
  virtual void prepare(LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
};

/*
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
  }

  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i=0; i<as.size(); i++)
      as[i].rhotot = as[i].rhoNew;
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
    for(int i=0; i<as.size(); i++)
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
    for(int i=0; i<as.size(); i++)
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
    a.xvalws[0]=a.xvalwsNew[0];
    a.xvalws[1]=a.xvalwsNew[1];
  }

  void updateChargeDensity(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i=0; i<as.size(); i++)
    {
      as[i].rhotot = as[i].rhoNew;
      as[i].xvalws[0]=as[i].xvalwsNew[0];
      as[i].xvalws[1]=as[i].xvalwsNew[1];
    }
  }

  void updatePotential(LSMSSystemParameters &lsms, AtomData &a)
  {
    simpleMixing( &a.vr(0,0), &a.vrNew(0,0), a.vr.size(), alpha);
  }

  void updatePotential(LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i=0; i<as.size(); i++)
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
*/

void setupMixing(MixingParameters &mix, Mixing* &mixing);


#endif
