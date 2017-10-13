#ifndef __UTLSOLIDBASE_H__
#define __UTLSOLIDBASE_H__

#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

enum par{ALL, HEADER, GENERATED, FLUX, SPD, MRPC};
 
class utlSolidBase
{

 public:

  utlSolidBase(std::string name);
  ~utlSolidBase();

  int Init();
  int End();
  void SetOutputName(std::string fname){ outfname = fname; }
  void AddFile(std::string fname){ filename = fname; }

  int process_tree();
  //  int process_tree(TTree* T);

  void InitContainers();
  void SetHeaderTreeVars(TTree* T);
  void SetGenTreeVars(TTree* T);
  void SetFluxTreeVars(TTree* T);
  void SetSPDTreeVars(TTree* T);
  void SetMRPCTreeVars(TTree* T);
  void ClearContainers(int mode);

 private:

  std::string outfname;
  std::string filename;
  TFile* fout ;

  bool ChkConversion(std::vector<double> *v_hit, std::vector<double> *v_id);

  double get_totEdep(std::vector<double> *v_hit, std::vector<double> *v_id, std::vector<double> *v_totE);

  TH1F* h1_edep;
  TH1F* h_edep_rate;
  TH2F* h2_edep;
  TH1F* h_edep_conv;
  TH1F* h_edep_conv0;
  TH1F* h_edep_conv_solo;
  TH1F* h_edep_conv1;
  TH1F* h_edep_conv2;
  TH1F* h_edep_pho;

 public:

  //Header
  std::vector<double> *evn;
  std::vector<double> *evn_type;
  std::vector<double> *beamPol;
  std::vector<double> *var1;
  std::vector<double> *var2;
  std::vector<double> *var3;
  std::vector<double> *var4;
  std::vector<double> *var5;
  std::vector<double> *var6;
  std::vector<double> *var7;
  std::vector<double> *var8;

  //Generated
  std::vector<int> *gen_pid;
  std::vector<double> *gen_px;
  std::vector<double> *gen_py;
  std::vector<double> *gen_pz;
  std::vector<double> *gen_vx;
  std::vector<double> *gen_vy;
  std::vector<double> *gen_vz;

  //Flux
  std::vector<double> *flux_hitn;
  std::vector<double> *flux_id;
  std::vector<double> *flux_pid;
  std::vector<double> *flux_mpid;
  std::vector<double> *flux_tid;
  std::vector<double> *flux_mtid;
  std::vector<double> *flux_otid;
  std::vector<double> *flux_trackE;
  std::vector<double> *flux_totEdep;
  std::vector<double> *flux_avg_x;
  std::vector<double> *flux_avg_y;
  std::vector<double> *flux_avg_z;
  std::vector<double> *flux_avg_lx;
  std::vector<double> *flux_avg_ly;
  std::vector<double> *flux_avg_lz;
  std::vector<double> *flux_px;
  std::vector<double> *flux_py;
  std::vector<double> *flux_pz;
  std::vector<double> *flux_vx;
  std::vector<double> *flux_vy;
  std::vector<double> *flux_vz;
  std::vector<double> *flux_mvx;
  std::vector<double> *flux_mvy;
  std::vector<double> *flux_mvz;
  std::vector<double> *flux_avg_t;
  
  //SPD Tree
  std::vector<double> *spd_id;
  std::vector<double> *spd_hitn;  
  std::vector<double> *spd_pid;
  std::vector<double> *spd_mpid;
  std::vector<double> *spd_tid;
  std::vector<double> *spd_mtid;
  std::vector<double> *spd_otid;
  std::vector<double> *spd_trackE;
  std::vector<double> *spd_totEdep;
  std::vector<double> *spd_avg_x;
  std::vector<double> *spd_avg_y;
  std::vector<double> *spd_avg_z;
  std::vector<double> *spd_avg_lx;
  std::vector<double> *spd_avg_ly;
  std::vector<double> *spd_avg_lz;
  std::vector<double> *spd_px;
  std::vector<double> *spd_py;
  std::vector<double> *spd_pz;
  std::vector<double> *spd_vx;
  std::vector<double> *spd_vy;
  std::vector<double> *spd_vz;
  std::vector<double> *spd_mvx;
  std::vector<double> *spd_mvy;
  std::vector<double> *spd_mvz;
  std::vector<double> *spd_avg_t;

  //MRPC Tree
  std::vector<double> *mrpc_pid;
  std::vector<double> *mrpc_mpid;
  std::vector<double> *mrpc_tid;
  std::vector<double> *mrpc_mtid;
  std::vector<double> *mrpc_otid;
  std::vector<double> *mrpc_trackE;
  std::vector<double> *mrpc_totEdep;
  std::vector<double> *mrpc_avg_x;
  std::vector<double> *mrpc_avg_y;
  std::vector<double> *mrpc_avg_z;
  std::vector<double> *mrpc_avg_lx;
  std::vector<double> *mrpc_avg_ly;
  std::vector<double> *mrpc_avg_lz;
  std::vector<double> *mrpc_px;
  std::vector<double> *mrpc_py;
  std::vector<double> *mrpc_pz;
  std::vector<double> *mrpc_vx;
  std::vector<double> *mrpc_vy;
  std::vector<double> *mrpc_vz;
  std::vector<double> *mrpc_mvx;
  std::vector<double> *mrpc_mvy;
  std::vector<double> *mrpc_mvz;
  std::vector<double> *mrpc_avg_t;
  std::vector<double> *mrpc_id;
  std::vector<double> *mrpc_hitn;

};

#endif /* __UTLSOLIDBASE_H__ */
