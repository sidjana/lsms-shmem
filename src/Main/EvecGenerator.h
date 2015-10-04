// mode: -*- c++ -*-

// Header for the evec generator class to be used in lsms.cc
// also example implementations:
// ConstantEvecGenerator
// RandomEvecGenerator

#ifndef EVEC_GENERATOR_H
#define EVEC_GENERATOR_H

#include <stdio.h>
#include <vector>
#include "random_evec.h"

class EvecGenerator
{
 public:
  virtual bool generateEvec(int instance, double *evecs, double energy) { return false; }
  virtual bool generateEvec(int instance, double *evecs, double energy, bool *accepted)
    { *accepted=true; return generateEvec(instance,evecs,energy); }
  virtual bool generateEvec(int instance, double *evecs, double energy, double magnetization, bool *accepted)
    { *accepted=true; return generateEvec(instance,evecs,energy); }
  virtual void initializeEvec(int instance, double *evecs) { generateEvec(instance,evecs,0.0); }
  virtual bool generateUnsampledEvec(int instance, double *evecs, double energy) {return generateEvec(instance,evecs,energy); }
  virtual void startSampling(void) {;}
  virtual void writeState(const char *name) {;}
  void setVerbosity(int v) {verbosity=v;}

  int verbosity;
};

class ConstantEvecGenerator : public EvecGenerator
{
 public:
  ConstantEvecGenerator(size_t num_spins, double evec[3], int nw=0)
  {n_spins=num_spins; setEvec(evec); n_walker=nw;
   if(nw>0) {walker_step.resize(nw); for(int i=0; i<nw; i++) walker_step[i]=0;}}
  void setEvec(double evec[3]){ev[0]=evec[0]; ev[1]=evec[1]; ev[2]=evec[2];}
  bool generateEvec(int instance, double *evecs, double energy)
  {
    if(n_walker>0) printf("Walker %4d Step %5d     Energy %.6lf\n",instance,walker_step[instance],energy);
    for(size_t j=0; j<3*n_spins; j+=3)
      {
        evecs[j]=ev[0];
        evecs[j+1]=ev[1];
        evecs[j+2]=ev[2];
      }
    if(instance<n_walker) walker_step[instance]++;
    return false;
  }
  void initializeEvec(int instance, double *evecs)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
      {
        evecs[j]=ev[0];
        evecs[j+1]=ev[1];
        evecs[j+2]=ev[2];
      }
  }
 private:
  size_t n_spins;
  double ev[3];
  int n_walker;
  std::vector<int> walker_step;
};

class RandomEvecGenerator : public EvecGenerator
{
 public:
  RandomEvecGenerator(size_t num_spins) {n_spins=num_spins;}
  bool generateEvec(int inst, double *evecs, double energy)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
      random_evec(&evecs[j]);
    return false;
  }
  void initializeEvec(int instance, double *evecs) { generateEvec(instance,evecs,0.0); }
 private:
  size_t n_spins;
};

#endif
