#ifndef __SOLMRPCDIGIT_H__
#define __SOLMRPCDIGIT_H__

#include "utlSolidBase.h"

#include <fstream>
#include <vector>

using namespace std;

enum mode{BEAM, GENERATOR, COSMIC, SINGLE};

struct otree
{

  /*
  double x, y, z, lz, edep;
  int evt_num;
  int nions;
  int ion;
  int id;
  double time0, timeref;
  double nElectrons;
  double nmult[200], charge[200], current[200], time1[200];
  double ind_charge;
  double sum_charge;
  */
  int EvtNumber;
  double rate;
  int pid;
  int mpid;
  int tid;
  int mtid;

  double gen_px;
  double gen_py;
  double gen_pz;
  double gen_vx;
  double gen_vy;
  double gen_vz;

  double trackE;
  double timeLeading;
  double timeOverThr;
  double timeRef;
  double totcharge;
  double fastCharge[200];

  double vp_x, vp_y, vp_z;
  double vp_vx, vp_vy, vp_vz;
  double vp_px, vp_py, vp_pz;
  int vp_module, vp_strip;

  int pid_in_gas[10];
  int mpid_in_gas[10];
  int tid_in_gas[10];
  int mtid_in_gas[10];

  double E[10];
  double edep[10];
  double charge[10];
  int nion[10];
  double x[10], y[10], z[10];
  double px[10], py[10], pz[10];
  int nAvalanche;
  int module[10], strip[10];

};

struct avlc_pair
{

  int istep;
  double nelec;
  double time;

};

class SolMRPCdigit : public utlSolidBase
{

 public:

  SolMRPCdigit(string name, string fname);
  ~SolMRPCdigit();

  void Init();
  void End();
  int process(int nevt, int ngen);
  void SetConstant(int mode);
  void SetOutFileName(string foutname){ outfilename = foutname; }
  void SetSourceMode(int mode){ _mode = mode; }
  void SetQthreshold(double Qthr){ Qthreshold = Qthr; }

 private:

  string fin;
  string outfilename;
  otree nt;

  void avalanche(int gap, double lz, vector<avlc_pair>& v_sig, double& ind_chg);
  double getSigma(double x);
  double prob(double x, double s);
  double calEweight();  
  void InitStruct();
  void SetOutTree(TTree* tree);

  int _mode;
  int _nsteps;
  double _stepsize;
  double Qthreshold;

  //avalanche sim parameter
  double alpha;
  double eta;
  double drift_v;
  double emin;
  double E_weight;


  TH1F* hToT;
  TH1F* hTimeLeading;
  TH1F* htotCharge;
  TH1F* hT0;
  TH1F* hRate;
  TH2F* h2ChargeTime;
  TH2F* h2ChargeToT;
  TH1F* hchannel;

  TFile* fout;
  TTree* T;

};

#endif /* __SOLMRPCDIGIT_H__ */
  
