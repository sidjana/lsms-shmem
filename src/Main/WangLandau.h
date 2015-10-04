// -*- mode: c++ -*-
#ifndef LSMS_WANG_LANDAU_H
#define LSMS_WANG_LANDAU_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// we use BOOST for the random number generator
 #include <boost/random.hpp>
//#include <random>
#include "../../mjson/json.h"
#include "EvecGenerator.h"
#include "Graph1dMoments.hpp"

void inline performGlobalUpdate(Graph1dMoments<double,double> &g, double kappa, double lambda, double omega)
{
  for(int i=0; i<g.getN(); i++)
    if(g[i]>omega) g[i]+=kappa*std::exp(-lambda/(g[i]-omega));
}

class StatesWriter
{
public:
  StatesWriter(const char *filename=NULL)
  {
    if(filename==NULL) writeFlag=false;
    else {writeFlag=true; of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(17);}
  }
  ~StatesWriter() {if(writeFlag) of.close();}
  void writeHeader(double lnF, int numWalk, int numSpin, double **spins)
  {
    if(writeFlag)
      {
	of<<lnF<<" "<<numWalk<<" "<<numSpin<<std::endl;
	for(int i=0; i<numWalk; i++)
	  {
	    of<<i;
	    for(int j=0; j<3*numSpin; j++) of<<" "<<spins[i][j];
	    of<<std::endl;
	  }
      }
  }
  void writeChange(int iWalk, int numRet, int ispin, double *ev, double E)
  {
    if(writeFlag)
      of<<iWalk<<" "<<numRet<<" "<<ispin<<" "<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<" "<<E<<std::endl;
  }
  void newFile(const char *filename=NULL)
  {
    if(writeFlag) of.close();
    writeFlag=false;
    if(filename!=NULL){
      writeFlag=true;
      of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);
    }
  }
private:
  bool writeFlag;
  std::ofstream of;
};

//template<class RNG = std::mt19937>
template<class RNG = boost::mt19937>
class WL1dEvecGenerator : public EvecGenerator
{
 public:
  WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                    const char *init_file_name=NULL, const char *out_file_name=NULL, const char *states_filename=NULL);
  void initializeEvec(int instance, double *evecs);
  bool generateUnsampledEvec(int instance, double *evecs, double energy) {initializeEvec(instance, evecs); return false;}
  void startSampling(void) { sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);}
  bool generateEvec(int instance, double *evecs, double energy, bool *accepted);
  bool generateEvec(int instance, double *evecs, double energy, double magnetization, bool *accepted);
  bool generateEvec(int instance, double *evecs, double energy) {bool h; return generateEvec(instance, evecs, energy, &h);}
  bool generateEvec(int instance, double *evecs) {std::cerr<<"Need energy for WL1dEvecGenerator\n"; exit(1);}
  void writeState(const char *name);
  void writeDos(const char *name);
 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
   boost::uniform_real<double> rnd; //rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  //std::uniform_real_distribution<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);

  /*
  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  */

  double hMinimum;
  Graph1dMoments<double,double> dos, histo;
  Kernel1d<double,double> dosKernel, histoKernel, nullKernel;
  KernelType kernelType;

  unsigned long accept, reject, acceptSinceLastChange;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  // instance specific:
  std::vector<long> ref0, ref1;
  std::vector<double> position, magnetizationAtPosition;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<int> lastAccepted;
  std::vector<double> lastAcceptedEnergy;
  std::vector<double> lastAcceptedEvec;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]

  std::vector<int> numRetentions;

  char *statesFile;
  StatesWriter sw;

  int changeMode;
  bool histogramUpdateMode;
  int updatesPerBin;

  struct {double kappa, lambda, omega; int frequency, changes;} globalUpdate;

#ifdef ISING
    void inline random_evec_1(double ev[3])
  {
    ev[0]=ev[1]=0.0;
    ev[2]=1.0;
    if(rng()%2 == 0) ev[2]=-ev[2];
  }
#else
  void inline random_evec_1(double ev[3])
  {
    double x,y,z;
    do {
      x = rnd(rng);
      y = rnd(rng);
    } while(x*x+y*y>1);
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r;
    y *= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z;
    r=1.0/sqrt(x*x+y*y+z*z);
    ev[0]=x*r; ev[1]=y*r; ev[2]=z*r;
  }
#endif

#ifdef ISING    
  void inline random_evec(double ev[3])
  {
    ev[2]=-ev[2];
  }
#else
  void inline random_evec(double ev[3])
  {
    double x, y, z;
    do {
      x = rnd(rng); y = rnd(rng);
    } while(x*x+y*y>1); 
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r; y*= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z; 
    // Project out the parallel component;
    r = x*ev[0] + y*ev[1] + z*ev[2];
    x -= r*ev[0]; y -= r*ev[1]; z -= r*ev[2];
    r = x*x + y*y + z*z;
    double t = 1-0.3*rnd(rng);
    ev[0] *= t; ev[1] *= t; ev[2] *= t;
    r = sqrt((1-t*t)/r);
    ev[0] += x*r; ev[1] += y*r; ev[2] += z*r;
    r=1.0/sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    ev[0]*=r; ev[1]*=r; ev[2]*=r;
    
    /*  
    ev[2]=1.0-2.0*rnd(rng);
    // ev[2]=rnd11(rng);
    double phi=2.0*M_PI*rnd(rng);
    // double phi=rnd0pi(rng);
    double cos_theta=sqrt(1-ev[2]*ev[2]);
    ev[0]=cos_theta*cos(phi);
    ev[1]=cos_theta*sin(phi);
    */  
  }
#endif
};

template<class RNG>
WL1dEvecGenerator<RNG>::WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
					  const char *init_file_name, const char *out_file_name, const char *states_filename)
  : sw(states_filename)
{
  verbosity=3;

  changeMode = 8+16;

  globalUpdate.frequency=0;
  globalUpdate.changes=0;
  globalUpdate.kappa=1.0;
  globalUpdate.lambda=1.0;
  globalUpdate.omega=0.5;

  long nX=-1;
  histogramUpdateMode=false;
  updatesPerBin=100;
  double interval=0.01;
  double kernelWidth=0.1;
  double xMin=-HUGE;
  double xMax=1.0;
  bool readPositions=false;
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  evecs_pointer = ev_p;
  ref0.resize(n_walkers);  for(int i=0; i<n_walkers; i++) ref0[i]=-1; // ref0[i]=HUGE;
  ref1.resize(n_walkers);
  position.resize(n_walkers);
  magnetizationAtPosition.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);
  for(int i=0; i<n_walkers; i++)
  {
    lastAccepted[i]=-2;
  }
  for(int i=0; i<3*n_walkers; i++)  lastAcceptedEvec[i]=0.0;

  oldSpin.resize(3*n_walkers);

  statesFile=NULL;
  if(states_filename!=NULL)
    {
      statesFile=(char *)malloc(sizeof(char)*(1+strlen(states_filename)));
      strcpy(statesFile,states_filename);
    }

  numRetentions.resize(n_walkers);
  for(int i=0; i<n_walkers; i++) numRetentions[i]=0;

  /*
  nX = -1;
  xMin = -HUGE; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  */
  dos_out_name="dos1d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  cycleCount=0;
  hMinimum= 1;     //   10
  acceptSinceLastChange=accept=reject=0;
  gammaFinal=1.e-6;
  flipPerUpdate=1; //  100
  updateCycle= 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.75;
  updatesPerBin =100;

  kernelType=None;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(interval);
  histo.setDeltaAndClear(interval);
  if(init_file_name!=NULL && init_file_name[0]!=0)
  {
    std::string label, value;

    std::cout<<"Reading Wang-Landau configuration from: "<<init_file_name<<std::endl;

    dos_out_name=std::string(init_file_name)+".out";
    if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
    std::ifstream inp(init_file_name);
    std::ostringstream buff;

    std::string line;
    while(std::getline(inp,line)) 
      buff << line << std::endl;

    inp.close();

    std::string fileString = buff.str();
    const char* fileChars  = fileString.c_str();
    json_t *json_root=NULL;

    json_parse_document(&json_root, (char *)fileChars);

    if(json_root == NULL || json_root->type != JSON_OBJECT)
    {
      std::ostringstream message;
      std::cerr << "In WL1dEvecGenerator(" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
  
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;
      if(label=="xMin") xMin=atof(it->child->text);
      else if(label=="xMax") xMax=atof(it->child->text);
      else if(label=="interval") interval=atof(it->child->text);
      else if(label=="kernelWidth") kernelWidth=atof(it->child->text);
      else if(label=="kernelType")
      {
        std::string strValue(it->child->text);
        kernelType=getKernelType(strValue);
      }
      else if(label=="gamma") gamma=atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal=atof(it->child->text);
      else if(label=="nX")
      {
        nX=atoi(it->child->text);
	dos.setRangeAndClear(xMin,xMax,nX);
	histo.setRangeAndClear(xMin,xMax,nX);
	/*
        if(dos!=NULL) free(dos);
        if(histo!=NULL) free(histo);
        dos=(double *)calloc(nX,sizeof(double));
        histo=(int *)calloc(nX,sizeof(int));
	*/
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="cycleCount") cycleCount=atoi(it->child->text);
      else if(label=="changeMode") changeMode=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
      else if(label=="histogramMinimum") hMinimum=atof(it->child->text);
      else if(label=="updatesPerBin") updatesPerBin=atoi(it->child->text);
      else if(label=="globalUpdate.frequency") globalUpdate.frequency=atoi(it->child->text);
      else if(label=="globalUpdate.changes") globalUpdate.changes=atoi(it->child->text);
      else if(label=="globalUpdate.kappa") globalUpdate.kappa=atof(it->child->text);
      else if(label=="globalUpdate.lambda") globalUpdate.lambda=atof(it->child->text);
      else if(label=="globalUpdate.omega") globalUpdate.omega=atof(it->child->text);
      else if(label=="seed") rng.seed(atoi(it->child->text));
      else if(label=="accept") accept=atol(it->child->text);
      else if(label=="acceptSinceLastChange") acceptSinceLastChange=atol(it->child->text);
      else if(label=="reject") reject=atol(it->child->text);
//*
      else if(label=="rngState")
      {
        std::string strValue(it->child->text);
        std::stringstream strStream(strValue, std::stringstream::in);
        strStream>>rng;
      }
//*/
      else if(label=="dos")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          dos[j++]=atof(i->text);
        }
        if(j!=dos.getN()) {std::cout<<"ERROR #(dos) "<<j<<" != nX "<<dos.getN()<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          histo[j++]=atof(i->text);
        }
        if(j!=histo.getN()) {std::cout<<"ERROR #(histo) "<<j<<" != nX "<<histo.getN()<<std::endl; exit(1);}
      }
      else if(label=="moments")
      {
        json_t *a = it->child->child;
        int k = atoi(a->text);
        dos.setNumberOfMoments(k);
        a=a->next;
        int jj=0;
        for(json_t *j=a->child; j!=NULL; j=j->next)
        {
          dos.setNumberOfSamplesAtIdx(jj,atoi(j->text));
          jj++;
        }
        json_t *i=a->next;
        for(int kk=0; kk<k; kk++)
        {
          int jj=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
          {
            dos.setMomentAtIdx(jj,kk,atof(j->text));
            jj++;
          }
          i=i->next;
        }
      }
      else if(label=="moments.k")
      {
        dos.setNumberOfMoments(atoi(it->child->text));
      }
      else if(label=="evecs")
      {
        json_t *a = it->child;
        n_initialized_from_file=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              evecs_pointer[n_initialized_from_file][k++]=atof(j->text);
            }
// initialize oldSpin and lastChange to point to site 0
            lastChange[n_initialized_from_file]=0;
            oldSpin[  3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][0];
            oldSpin[1+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][1];
            oldSpin[2+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][2];
//
            n_initialized_from_file++;
            if(k!=3*n_spins) {std::cout<<"ERROR #(evecs) "<<k<<" != 3*n_spins "<<3*n_spins<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="oldSpin")
      {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              oldSpin[n_initialized*3+k++]=atof(j->text);
            }
            n_initialized++;
            if(k!=3) {std::cout<<"ERROR #(oldSpin) "<<k<<" != 3"<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="lastChange")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers) lastChange[j++]=atoi(i->text);
          n_initialized++;
        }
      }
      else if(label=="position")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            position[j++]=atof(i->text);
            n_initialized++;
          }
        }
        readPositions=true;
      }
      else if(label=="magnetizationAtPosition")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            magnetizationAtPosition[j++]=atof(i->text);
            n_initialized++;
          }
        }
//        readMagnetizationAtPositions=true;
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else std::cout<<"WARNING: unknown label: "<<label<<std::endl;
    }

    // set xMax and xMin or interval depending on nX:
    if(nX>0)
    {
      interval=(xMax-xMin)/double(nX);
      dos.setRange(xMin,xMax); histo.setRange(xMin,xMax);
    } else {
      dos.setDeltaAndClear(interval); histo.setDeltaAndClear(interval);
    }

    json_free_value(&json_root);
  }
  dosKernel.setWidthAndClear(interval,kernelWidth);
  histoKernel.setWidthAndClear(interval,kernelWidth);
  nullKernel.setWidthAndClear(interval,kernelWidth);
  initKernel(kernelType,dosKernel);
  dosKernel.scale(gamma/dosKernel(0.0));
  initKernel(kernelType,histoKernel);
  histoKernel.scale(1.0/histoKernel(0.0));
  initKernel(kernelType,nullKernel);
  nullKernel.scale(0.0);
//  histoKernel.scale(interval);

  if(readPositions)
    for(int i=0; i<n_walkers; i++)
      ref0[i]=dos.idx(position[i]);

  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout<<"Wang-Landau output will be written to: "<<dos_out_name<<std::endl;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  bool firstFerromagnetic=false;
  if(inst>=n_initialized_from_file)
    {
      if(firstFerromagnetic)
        if(inst==0)
	  for(size_t j=0; j<3*n_spins; j+=3)
	    {evecs[j]=0.0; evecs[j+1]=0.0; evecs[j+2]=1.0;}
        else if(inst==1)
          for(size_t j=0; j<3*n_spins; j+=3)
            {evecs[j+2]= (j/3)%2 ? 1.0 : -1.0 ; evecs[j+1]=0.0; evecs[j]=0.0;}
        else
	  for(size_t j=0; j<3*n_spins; j+=3)
	    random_evec_1(&evecs[j]);
      else
        for(size_t j=0; j<3*n_spins; j+=3)
          random_evec_1(&evecs[j]);
    }

  out[inst]=false;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::writeState(const char* name)
{
  // if(!syncronizeGraphs(dos,histo))
  if(dos.getN()!=histo.getN())
  {
    std::cout<<"Histogramm size dosn't match DOS! Clearing histogramm!\n";
    histo.setRangeAndClear(dos.getMinX(),dos.getMaxX(),dos.getN());
  }
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(8);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    ofile<<"\"nX\" : " << dos.getN() << ",\n";
    if(kernelType!=None)
      ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma << ",\n";
    ofile<<"\"accept\" : " << accept <<",\n";
    ofile<<"\"acceptSinceLastChange\" : " << acceptSinceLastChange <<",\n";
    ofile<<"\"reject\" : " << reject <<",\n";
    ofile<<"\"gammaFinal\" : " << gammaFinal << ",\n";

    ofile<<"\"changeMode\" : "<< changeMode <<",\n";
    ofile<<"\"flatnessCriterion\" : "<< flatnessCriterion <<",\n";
    ofile<<"\"histogramMinimum\" : "<< hMinimum <<",\n";
    ofile<<"\"updatesPerBin\" : "<< updatesPerBin <<",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << ",\n";
    ofile<<"\"updateCycle\" : " << updateCycle << ",\n";

    ofile<<"\"globalUpdate.frequency\" : "<<  globalUpdate.frequency << ",\n";
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"globalUpdate.kappa\" : "<<  globalUpdate.kappa << ",\n";
    ofile<<"\"globalUpdate.lambda\" : "<<  globalUpdate.lambda << ",\n";
    ofile<<"\"globalUpdate.omega\" : "<<  globalUpdate.omega << ",\n";

    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"rngState\" : \""<<rng<<"\",\n";
    ofile<<"\"evecs\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[\n";
      for(int j=0; j<3*n_spins; j+=3)
        ofile<<evecs_pointer[i][j]<<", "<<evecs_pointer[i][j+1]
             <<", "<<evecs_pointer[i][j+2]<<((j==3*n_spins-3)?"\n":",\n");
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"oldSpin\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[ "<<oldSpin[3*i]<<", "<<oldSpin[3*i+1]
             <<", "<<oldSpin[3*i+2];
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"lastChange\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< lastChange[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"position\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< position[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"magnetizationAtPosition\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< magnetizationAtPosition[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"stepsSinceLastHistogramUpdate\" : " << stepsSinceLastHistogramUpdate << ",\n";
    ofile<<"\"numberOfUpdatesSinceLastBoost\" : " << numberOfUpdatesSinceLastBoost << ",\n";
    ofile<<"\"modificationFactorChanges\" : " << modificationFactorChanges << ",\n";
    ofile<<"\"cycleCount\" : " << cycleCount << "\n";
    ofile<<"}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
} 

template<class RNG>
void WL1dEvecGenerator<RNG>::writeDos(const char* name)
{
  std::ofstream ofile(name);
  // write dos;
  // we are using JSON as our file format
  if(ofile)
  {
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << "," <<std::endl;
    ofile<<"\"xMax\" : " << dos.getMaxX() << "," <<std::endl;
    ofile<<"\"nX\" : " << dos.getN() << "," <<std::endl;
    ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma <<"," <<std::endl;
    ofile<<"\"gammaFinal\" : " << gammaFinal << "," <<std::endl;
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << "," <<std::endl;
    ofile<<"\"updateCycle\" : " << updateCycle << "," <<std::endl;
    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"]\n}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
}

template<class RNG>
bool WL1dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, double energy, bool *accepted)
{
  return generateEvec(instance, evecs, energy, 0.0, accepted);
}

template<class RNG>
bool WL1dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, double energy, double magnetization, bool *accepted)
{
  bool accept_step;
  // +++++++++ Added by Odbadrakh and Don on Aug 11, 2010
  // +++ If a flip of spin results in energy in the same bin as previous accepted 
  // +++ state, it is accepted all the time. We change it by accepting if the
  // +++ new state has lower density of states according to a linear spline of DOS
  
  double dos_differ;
  double energy_differ;
  double to_go_or_not;

  // +++++++++ end of modification. Follow the new variables to see the actual changes 
  // energy between xMin and xMax ?
  int grid = dos.idx(energy);

  // dos.test();

  stepsSinceLastHistogramUpdate++; // counter[0]++
  *accepted=false;
  // record initial energy for step out
  if(lastAccepted[instance]==-2)
    {position[instance]=lastAcceptedEnergy[instance]=energy; lastAccepted[instance]=-1;}
  
  if(grid < 0 || grid > dos.getN()-1)
  {
    dos.extendTo(energy-dosKernel.getWidth()); histo.extendTo(energy-histoKernel.getWidth());
    dos.extendTo(energy+dosKernel.getWidth()); histo.extendTo(energy+histoKernel.getWidth());
    grid=dos.idx(energy);
    // we need to adjust the grid point of the walker positions (ref0) for all walkers
    for(long ii=0; ii<n_walkers; ii++)
    {
      if(ref0[ii]>=0) ref0[ii]=dos.idx(position[ii]);
    }
  }
  numRetentions[instance]++;
  out[instance]=false;
  ref1[instance] = grid;
  if(ref0[instance]<0)
    {
      accept_step = true;
      ref0[instance] = ref1[instance];
      position[instance] = energy;
      magnetizationAtPosition[instance]=magnetization;
    } else {
    if(abs(ref1[instance] - ref0[instance]) < 1 ) {
// Actual change made by Odbadrakh, Aug 30, 2010
      dos_differ = dos[ref0[instance]-1] - dos[ref0[instance]];
      energy_differ = energy-lastAcceptedEnergy[instance];
      to_go_or_not = dos_differ * energy_differ;
      if(to_go_or_not >= 0.0) { // accepts all downhill changes
        accept_step = true;
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance]=magnetization;
      } else { // uphill moves
        if(rnd(rng) < exp(to_go_or_not/dos.getDelta()))
        { //std::cout<<" dos "<<instance<<"  weight "<<dos.getDelta()<<std::endl;
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance]=magnetization;
        } else {
          accept_step = false;
        }
      }
    } else {

      if(dos[ref1[instance]] <= dos[ref0[instance]])
      {
        accept_step = true;
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance]=magnetization;
      } else {
        if(rnd(rng) < exp(dos[ref0[instance]] - dos[ref1[instance]]))
        {
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance]=magnetization;
        } else {
          accept_step = false;
        }
      }
      
    }

 }

// End of change made on Aug 30, 2010
  *accepted=accept_step;
  if(verbosity>2)
  std::cout<<"WangLandau 1d EvecGenerator step "
           <<modificationFactorChanges<<":"<<numberOfUpdatesSinceLastBoost<<":"
           <<stepsSinceLastHistogramUpdate<<" nX="<<dos.getN()
	   <<" ["<<dos.getMinX()<<", "<<dos.getMaxX()<<"] "
           <<(accept_step ? "accepted" : "rejected")
           <<" Energy = "<<energy<<", Instance "<<instance<<std::endl;

  if(stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate=0; // counter[0]
    numberOfUpdatesSinceLastBoost++; // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      // addKernel(dos,dosKernel,energy);
      if(!histogramUpdateMode) // standard Wang-Landau
      {
        addKernel(dos,dosKernel,position[instance]);
        dos.addMomentsAtIdx(dos.idx(position[instance]),magnetizationAtPosition[instance]);
        // if(accept_step) addKernel(histo,histoKernel,position[instance]);
        addKernel(histo,histoKernel,position[instance]);
        // addKernel(dos,dosKernel,energy);
        // addKernel(histo,histoKernel,energy);
      } else {
        // addKernel(dos,nullKernel,position[instance]);
        if(accept_step) addKernel(histo,histoKernel,position[instance]);
      }
    } else {
      std::cerr<<"ATTENTION: We should never reach this place in WL1dEvecGenerator!\n";
      exit(1);
    }
  }
  if(changeMode & (8+16+32)) // change gamma at every step
  {
    long n=accept+reject;
    double dn;
    if(changeMode &  (8+16) == 8) n=accept;
    else if(changeMode &  (8+16) == 16) n=reject;
    else if(changeMode &  (8+16) == 8+16) n=accept+reject;
    if(changeMode & 32) dn=double(n)/double(n_walkers);
    else dn=double(n);
    gamma=1.0/dn; if(gamma>1.0) gamma=1.0;
    initKernel(kernelType,dosKernel);
    dosKernel.scale(gamma);
  }
  if(cycleCount >= updateCycle)
  {
    cycleCount=0;
    // syncronizeGraphs(dos,histo);
    if(dos.getN()!=histo.getN())
    {
      std::cout<<"Histogramm size dosn't match DOS! Clearing histogramm!\n";
      histo.setRangeAndClear(dos.getMinX(),dos.getMaxX(),dos.getN());
    }
    writeState("WL1d.state");
    if(!histogramUpdateMode) // Normal Wang-Landau
    {
    // calculate minimum nonzero histogram
      double hMin, hMax, hMean;
    // we look only at the histogram inside the energy interval that was actually sampled if we use kernel updates
      if(kernelType==None)
      {
        hMean=histo.getMeanY();
        histo.getMinMaxY(hMin,hMax);
      } else {
        hMean=histo.getMeanYInInterval(dos.getMinX()+dosKernel.getWidth(),
                                       dos.getMaxX()-dosKernel.getWidth());
        histo.getMinMaxYInInterval(dos.getMinX()+dosKernel.getWidth(),
                                   dos.getMaxX()-dosKernel.getWidth(),
                                   hMin, hMax);
      }
    // double currentFlatnessCriterion = double(hMin-hMax)/hMean;
    double currentFlatnessCriterion = double(hMin)/hMean;
    
    std::cout <<"# accepence ratio = "<<double(accept)/double(accept+reject)<<"\n";
    std::cout <<"# current flatness = "<<currentFlatnessCriterion<<(changeMode & 4 ? " *":"")<<"\n";
    std::cout <<"# current histogram minimum = "<<hMin<<(changeMode & 2 ? " *":"")<<"\n";
    std::cout <<"# average accepted steps/bin since last gamma change = "
              <<double(acceptSinceLastChange)/double(histo.getN())
              <<(changeMode & 1 ? " *":"")<<"\n";

    if(changeMode != 0 && changeMode<8)
    {
      if(globalUpdate.frequency>0 && (globalUpdate.frequency*histo.getN())<acceptSinceLastChange) //global update
      {
        char fn[256];
        snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
        writeDos(fn);
        globalUpdate.changes++;
        std::cout<<"# global update "<<globalUpdate.changes<<std::endl;
        performGlobalUpdate(dos,globalUpdate.kappa,globalUpdate.lambda,globalUpdate.omega);
        histo.clear();
        acceptSinceLastChange=0;
      }
      else if(((acceptSinceLastChange>=histo.getN()*updatesPerBin || !(changeMode &1)) &&
          (hMin >= hMinimum || !(changeMode & 2)) &&
          (currentFlatnessCriterion>flatnessCriterion || !(changeMode & 4)) ))
      {
        std::cout <<"# level "<<modificationFactorChanges<<" with gamma="<<gamma<<" is finished.\n";
        modificationFactorChanges++; // counter[2]

        // write dos;
        // we are using JSON as our file format
        char fn[256];
        snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
        writeDos(fn);
        // writeDos(dos_out_name.data());
        // clear the Histogram
        histo.clear();
        acceptSinceLastChange=0;
        // change gamma
        dosKernel.scale(0.5);
        gamma = 0.5*gamma;
        std::cout<<"# level "<<modificationFactorChanges<<" with gamma="<<gamma<<" begins.\n";
        if(statesFile!=NULL)
	{
	  // char fn[256];
	  snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
	  sw.newFile(fn);
	  sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
	}
      }
    }
    } else { // histogram Update mode
      if(acceptSinceLastChange>=histo.getN()*updatesPerBin)
      {
        acceptSinceLastChange=0;
        for(int ix=0; ix<dos.getN(); ix++) dos[ix]+=gamma*histo[ix];
        histo.clear();
        gamma = 0.5*gamma;

        std::cout<<"# level "<<modificationFactorChanges<<" with gamma="<<gamma<<" begins.\n";
        if(statesFile!=NULL)
        {
          char fn[256];
          snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
          sw.newFile(fn);
          sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
        }
      }
    }
  }

  if(gamma<gammaFinal) return true;

  if(accept_step)
  {
    sw.writeChange(instance,numRetentions[instance],lastAccepted[instance], &lastAcceptedEvec[3*instance],
                   lastAcceptedEnergy[instance]);
    lastAccepted[instance]=lastChange[instance];
    lastAcceptedEnergy[instance]=energy;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]];
    accept++;
    acceptSinceLastChange++;
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
    numRetentions[instance]=0;
  } else {
    reject++;
    evecs[  3*lastChange[instance]] = oldSpin[  3*instance];
    evecs[1+3*lastChange[instance]] = oldSpin[1+3*instance];
    evecs[2+3*lastChange[instance]] = oldSpin[2+3*instance];
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
  }

  return false;
}


#endif

