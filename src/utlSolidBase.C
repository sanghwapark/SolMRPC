#include <iostream>

//ROOT headers
#include <TBranch.h>
#include <TMath.h>

#include "utlSolidBase.h"
#include "SolidConst.h"

using namespace std;

utlSolidBase::utlSolidBase(string name)
{

  outfname=name;
  //  ClearContainers();

}

//___________________________________________________________


utlSolidBase::~utlSolidBase()
{
}

//___________________________________________________________

int utlSolidBase::Init()
{

  fout = new TFile(outfname.c_str(), "RECREATE");

  h1_edep = new TH1F("h1_edep","h1_edep", 400, 0, 20);
  h2_edep = new TH2F("h2_edep","h2_edep", 110, 0, 11, 400, 0, 20);

  h_edep_rate = new TH1F("h_edep","h_edep", 400, 0, 20);
  h_edep_pho = new TH1F("h_edep_pho","Edep from photons to SPD", 400, 0, 20);
  h_edep_conv = new TH1F("h_edep_conv","Edep from e+ and e- pairs", 400, 0, 20);
  h_edep_conv0 = new TH1F("h_edep_conv0","Edep from e+ or e- pairs", 400, 0, 20);
  h_edep_conv1 = new TH1F("h_edep_conv1","Edep from one e+ and e- pairs", 400, 0, 20);
  h_edep_conv2 = new TH1F("h_edep_conv2","Edep from two e+ and e- pairs", 400, 0, 20);
  h_edep_conv_solo = new TH1F("h_edep_conv_solo","Edep from single e", 400, 0, 20);

  return 0;
}

//___________________________________________________________

int utlSolidBase::End()
{
  fout->cd();
  h1_edep->Write();
  h2_edep->Write();  
  h_edep_rate->Write();
  h_edep_pho->Write();
  h_edep_conv->Write();
  h_edep_conv0->Write();
  h_edep_conv1->Write();
  h_edep_conv2->Write();
  h_edep_conv_solo->Write();

  fout->Close();
  return 0;

}

//___________________________________________________________
bool utlSolidBase::ChkConversion(std::vector<double> *v_hit, 
		     std::vector<double> *v_id)
{

  bool conv = 0;
  for(int i=0; i<(int)v_hit->size(); i++)
    {
      int pid = (int)v_id->at(i);
      if(abs(pid) == ELECTRON || abs(pid) == PION || abs(pid) == KAON || abs(pid) == PROTON)
	conv = 1;
    }

  return conv;
}

//___________________________________________________________

int utlSolidBase::process_tree()
{

  std::cout << "READ: " << filename.c_str() << std::endl;

  TFile* fin = new TFile(filename.c_str());
  if(!fin->IsOpen())
    {
      cerr << "Input file is not open" << endl;
      return 1;
    }

  InitContainers();

  TTree* T_header = (TTree*)fin->Get("header");
  SetHeaderTreeVars(T_header);
  TTree* T_gen = (TTree*)fin->Get("generated");
  SetGenTreeVars(T_gen);
  TTree* T_flux = (TTree*)fin->Get("flux");
  SetFluxTreeVars(T_flux);
  TTree* T_spd = (TTree*)fin->Get("solid_spd");
  SetSPDTreeVars(T_spd);  

  double r_min = 96;
  double r_max = 210;

  Long64_t nentries = T_gen->GetEntries();
  for( Long64_t ientry=0; ientry<nentries; ientry++)
    {
      if(ientry%100000==0) cout << ientry << " of " << nentries << endl;
      ClearContainers(ALL);

      T_header->GetEntry(ientry);
      T_spd->GetEntry(ientry);
      T_flux->GetEntry(ientry);

      double rate = var8->at(0);
      double edep = 0;
      edep = get_totEdep(spd_hitn, spd_id, spd_totEdep);

      h1_edep->Fill(edep);
      h_edep_rate->Fill(edep, rate);

      for(int i=0; i<(int)flux_hitn->size(); i++)
	{

	  int vp_id = (int)flux_id->at(i);
	  int vp_tid = (int)flux_tid->at(i);
	  int vp_pid = (int)flux_pid->at(i);

	  double x = flux_avg_x->at(i) * 1.e-1;
	  double y = flux_avg_y->at(i) * 1.e-1;
	  double r = sqrt(pow(x,2)+pow(y,2));

	  if(vp_id != 5110000) continue;
	  if(r<r_min || r>r_max) continue;
	  if(flux_pz->at(i) < 0) continue;

	  // if(vp_pid != PION) continue;
	  // if(abs(vp_pid) != PION) continue;
	  if(abs(vp_pid) != ELECTRON) continue;

	  double edep2 = 0;
	  bool basic_cut = false;
	  bool conv = false;
	  bool find_ele = false; bool find_pos = false;
	  int ne=0; int np = 0;
	  int nconv=0;
	  for(int j=0; j<(int)spd_hitn->size(); j++)
	    {
	      int det_id = (int)spd_id->at(j);
	      if(det_id != 5100000) continue; //FASPD

	      int det_pid = (int)spd_pid->at(j);
	      int det_tid = (int)spd_tid->at(j);
	      int det_mpid = (int)spd_mpid->at(j);
	      int det_mtid = (int)spd_mtid->at(j);
	      
	      double det_x = spd_avg_x->at(j) * 1.e-1;
	      double det_y = spd_avg_y->at(j) * 1.e-1;
	      double det_r = sqrt(pow(det_x,2)+pow(det_y,2));

	      if(det_tid!=vp_tid && det_mtid!=vp_tid) continue;
	      if(det_pid!=vp_pid && det_mpid!=vp_pid) continue;
	      if(det_r > r_max || det_r < r_min) continue;
	      basic_cut = true;

	      if(det_pid == 11)
		{
		  find_ele = true;
		  ne++;
		}
	      if(det_pid == -11)
		{
		  find_pos = true;
		  np++;
		}

	      edep2 += spd_totEdep->at(j);
	    }

	  if(basic_cut) h_edep_pho->Fill(edep2);

	  if(find_ele || find_pos) h_edep_conv0->Fill(edep2);
	  if((ne+np)==1) h_edep_conv_solo->Fill(edep2);

	  if(find_ele && find_pos) conv = true;
	  if(ne==1 && np==1) nconv = 1;
	  else if(ne==2 && np==2) nconv = 2;
	  if(conv) h_edep_conv->Fill(edep2);

	  if(nconv==1) h_edep_conv1->Fill(edep2);
	  if(nconv==2) h_edep_conv2->Fill(edep2);

	}

    }

  fin->Close();

  return 0;
}

//___________________________________________________________

double utlSolidBase::get_totEdep(std::vector<double> *v_hit, std::vector<double> *v_id, std::vector<double> *v_totE)
{

  double totEdep = 0;
  for(int i=0; i<(int)v_hit->size(); i++)
    {

      /*
      cout << (int)v_id->at(i) << endl;
      int det_ID = ((int)v_id->at(i))/1000000;
      int subdet_ID = ((int)v_id->at(i)%1000000)/100000;
      int subsubdet_ID = (((int)v_id->at(i)%1000000)%1000000)/10000;
      int comp_ID = (int)v_id->at(i)%10000;
      if(det_ID == 5 && subdet_ID == 1 && subsubdet_ID == 0)
      */
      int det_ID = (int)v_id->at(i);
      if( det_ID == 5100000)
	{
	  totEdep += v_totE->at(i);
	}
    }

  return totEdep;
}

//___________________________________________________________

void utlSolidBase::ClearContainers(int mode)
{
  if(mode == HEADER || mode == ALL)
    {
      evn->clear();
      evn_type->clear();
      beamPol->clear();
      var1->clear();
      var2->clear();
      var3->clear();
      var4->clear();
      var5->clear();
      var6->clear();
      var7->clear();
      var8->clear();
    }
  else if(mode == GENERATED || mode == ALL)
    {
  gen_pid->clear();
  gen_px->clear();
  gen_py->clear();
  gen_pz->clear();
  gen_vx->clear();
  gen_vy->clear();
  gen_vz->clear();
    }
  else if(mode == FLUX || mode == ALL)
    {
  flux_hitn->clear();
  flux_id->clear();
  flux_pid->clear();
  flux_mpid->clear();
  flux_tid->clear();
  flux_otid->clear();
  flux_trackE->clear();
  flux_totEdep->clear();
  flux_avg_x->clear();
  flux_avg_y->clear();
  flux_avg_z->clear();
  flux_avg_lx->clear();
  flux_avg_ly->clear();
  flux_avg_lz->clear();
  flux_px->clear();
  flux_py->clear();
  flux_pz->clear();
  flux_vx->clear();
  flux_vy->clear();
  flux_vz->clear();
  flux_mvx->clear();
  flux_mvy->clear();
  flux_mvz->clear();
  flux_avg_t->clear();
    }  
  else if(mode == SPD || mode == ALL)
    {
  spd_id->clear();
  spd_hitn->clear();  
  spd_pid->clear();
  spd_mpid->clear();
  spd_tid->clear();
  spd_mtid->clear();
  spd_otid->clear();
  spd_trackE->clear();
  spd_totEdep->clear();
  spd_avg_x->clear();
  spd_avg_y->clear();
  spd_avg_z->clear();
  spd_avg_lx->clear();
  spd_avg_ly->clear();
  spd_avg_lz->clear();
  spd_px->clear();
  spd_py->clear();
  spd_pz->clear();
  spd_vx->clear();
  spd_vy->clear();
  spd_vz->clear();
  spd_mvx->clear();
  spd_mvy->clear();
  spd_mvz->clear();
  spd_avg_t->clear();
    }
  else if(mode == MRPC || mode == ALL)
    {
  mrpc_id->clear();
  mrpc_hitn->clear();  
  mrpc_pid->clear();
  mrpc_mpid->clear();
  mrpc_tid->clear();
  mrpc_mtid->clear();
  mrpc_otid->clear();
  mrpc_trackE->clear();
  mrpc_totEdep->clear();
  mrpc_avg_x->clear();
  mrpc_avg_y->clear();
  mrpc_avg_z->clear();
  mrpc_avg_lx->clear();
  mrpc_avg_ly->clear();
  mrpc_avg_lz->clear();
  mrpc_px->clear();
  mrpc_py->clear();
  mrpc_pz->clear();
  mrpc_vx->clear();
  mrpc_vy->clear();
  mrpc_vz->clear();
  mrpc_mvx->clear();
  mrpc_mvy->clear();
  mrpc_mvz->clear();
  mrpc_avg_t->clear();
    }

  return;

}

//___________________________________________________________
void utlSolidBase::InitContainers()
{

  //vectors need to be initialized first
  evn = 0;
  evn_type = 0;
  beamPol = 0;
  var1 = 0;
  var2 = 0;
  var3 = 0;
  var4 = 0;
  var5 = 0;
  var6 = 0;
  var7 = 0;
  var8 = 0;

  gen_pid = 0;
  gen_px = 0;
  gen_py = 0;
  gen_pz = 0;
  gen_vx = 0;
  gen_vy = 0;
  gen_vz = 0;

  flux_hitn = 0;
  flux_id = 0;
  flux_pid = 0;
  flux_mpid = 0;
  flux_tid = 0;
  flux_otid = 0;
  flux_trackE = 0;
  flux_totEdep = 0;
  flux_avg_x = 0;
  flux_avg_y = 0;
  flux_avg_z = 0;
  flux_avg_lx = 0;
  flux_avg_ly = 0;
  flux_avg_lz = 0;
  flux_px = 0;
  flux_py = 0;
  flux_pz = 0;
  flux_vx = 0;
  flux_vy = 0;
  flux_vz = 0;
  flux_mvx = 0;
  flux_mvy = 0;
  flux_mvz = 0;
  flux_avg_t = 0;
  
  spd_id = 0;
  spd_hitn = 0;  
  spd_pid = 0;
  spd_mpid = 0;
  spd_tid = 0;
  spd_mtid = 0;
  spd_otid = 0;
  spd_trackE = 0;
  spd_totEdep = 0;
  spd_avg_x = 0;
  spd_avg_y = 0;
  spd_avg_z = 0;
  spd_avg_lx = 0;
  spd_avg_ly = 0;
  spd_avg_lz = 0;
  spd_px = 0;
  spd_py = 0;
  spd_pz = 0;
  spd_vx = 0;
  spd_vy = 0;
  spd_vz = 0;
  spd_mvx = 0;
  spd_mvy = 0;
  spd_mvz = 0;
  spd_avg_t = 0;

  mrpc_id = 0;
  mrpc_hitn = 0;  
  mrpc_pid = 0;
  mrpc_mpid = 0;
  mrpc_tid = 0;
  mrpc_mtid = 0;
  mrpc_otid = 0;
  mrpc_trackE = 0;
  mrpc_totEdep = 0;
  mrpc_avg_x = 0;
  mrpc_avg_y = 0;
  mrpc_avg_z = 0;
  mrpc_avg_lx = 0;
  mrpc_avg_ly = 0;
  mrpc_avg_lz = 0;
  mrpc_px = 0;
  mrpc_py = 0;
  mrpc_pz = 0;
  mrpc_vx = 0;
  mrpc_vy = 0;
  mrpc_vz = 0;
  mrpc_mvx = 0;
  mrpc_mvy = 0;
  mrpc_mvz = 0;
  mrpc_avg_t = 0;

  return;

}

//___________________________________________________________

void utlSolidBase::SetMRPCTreeVars(TTree* T)
{

  T->SetBranchAddress("pid",&mrpc_pid);
  T->SetBranchAddress("mpid",&mrpc_mpid);
  T->SetBranchAddress("tid",&mrpc_tid);
  T->SetBranchAddress("mtid",&mrpc_mtid);
  T->SetBranchAddress("otid",&mrpc_otid);
  T->SetBranchAddress("trackE",&mrpc_trackE);
  T->SetBranchAddress("totEdep",&mrpc_totEdep);
  T->SetBranchAddress("avg_x",&mrpc_avg_x);
  T->SetBranchAddress("avg_y",&mrpc_avg_y);
  T->SetBranchAddress("avg_z",&mrpc_avg_z);
  T->SetBranchAddress("avg_lx",&mrpc_avg_lx);
  T->SetBranchAddress("avg_ly",&mrpc_avg_ly);
  T->SetBranchAddress("avg_lz",&mrpc_avg_lz);
  T->SetBranchAddress("px",&mrpc_px);
  T->SetBranchAddress("py",&mrpc_py);
  T->SetBranchAddress("pz",&mrpc_pz);
  T->SetBranchAddress("vx",&mrpc_vx);
  T->SetBranchAddress("vy",&mrpc_vy);
  T->SetBranchAddress("vz",&mrpc_vz);
  T->SetBranchAddress("mvx",&mrpc_mvx);
  T->SetBranchAddress("mvy",&mrpc_mvy);
  T->SetBranchAddress("mvz",&mrpc_mvz);
  T->SetBranchAddress("avg_t",&mrpc_avg_t);
  T->SetBranchAddress("id",&mrpc_id);
  T->SetBranchAddress("hitn",&mrpc_hitn);

  return;
}
//___________________________________________________________

void utlSolidBase::SetHeaderTreeVars(TTree* T)
{

  T->SetBranchAddress("evn",&evn);
  T->SetBranchAddress("evn_type",&evn_type);
  T->SetBranchAddress("beamPol",&beamPol);
  T->SetBranchAddress("var1",&var1);
  T->SetBranchAddress("var2",&var2);
  T->SetBranchAddress("var3",&var3);
  T->SetBranchAddress("var4",&var4);
  T->SetBranchAddress("var5",&var5);
  T->SetBranchAddress("var6",&var6);
  T->SetBranchAddress("var7",&var7);
  T->SetBranchAddress("var8",&var8);

 return;
}

//___________________________________________________________

void utlSolidBase::SetGenTreeVars(TTree* T)
{
  T->SetBranchAddress("pid",&gen_pid);
  T->SetBranchAddress("px",&gen_px);
  T->SetBranchAddress("py",&gen_py);
  T->SetBranchAddress("pz",&gen_pz);
  T->SetBranchAddress("vx",&gen_vx);
  T->SetBranchAddress("vy",&gen_vy);
  T->SetBranchAddress("vz",&gen_vz);
 
  return;
}

//___________________________________________________________

void utlSolidBase::SetFluxTreeVars(TTree* T)
{

  T->SetBranchAddress("hitn",&flux_hitn);
  T->SetBranchAddress("id",&flux_id);
  T->SetBranchAddress("pid",&flux_pid);
  T->SetBranchAddress("mpid",&flux_mpid);
  T->SetBranchAddress("tid",&flux_tid);
  T->SetBranchAddress("mtid",&flux_mtid);
  T->SetBranchAddress("otid",&flux_otid);
  T->SetBranchAddress("trackE",&flux_trackE);
  T->SetBranchAddress("totEdep",&flux_totEdep);
  T->SetBranchAddress("avg_x",&flux_avg_x);
  T->SetBranchAddress("avg_y",&flux_avg_y);
  T->SetBranchAddress("avg_z",&flux_avg_z);
  T->SetBranchAddress("avg_lx",&flux_avg_lx);
  T->SetBranchAddress("avg_ly",&flux_avg_ly);
  T->SetBranchAddress("avg_lz",&flux_avg_lz);
  T->SetBranchAddress("px",&flux_px);
  T->SetBranchAddress("py",&flux_py);
  T->SetBranchAddress("pz",&flux_pz);
  T->SetBranchAddress("vx",&flux_vx);
  T->SetBranchAddress("vy",&flux_vy);
  T->SetBranchAddress("vz",&flux_vz);
  T->SetBranchAddress("mvx",&flux_mvx);
  T->SetBranchAddress("mvy",&flux_mvy);
  T->SetBranchAddress("mvz",&flux_mvz);
  T->SetBranchAddress("avg_t",&flux_avg_t);

  return;
}

//___________________________________________________________

void utlSolidBase::SetSPDTreeVars(TTree* T)
{

  T->SetBranchAddress("hitn",&spd_hitn);
  T->SetBranchAddress("id",&spd_id);
  T->SetBranchAddress("pid",&spd_pid);
  T->SetBranchAddress("mpid",&spd_mpid);
  T->SetBranchAddress("tid",&spd_tid);
  T->SetBranchAddress("mtid",&spd_mtid);
  T->SetBranchAddress("otid",&spd_otid);
  T->SetBranchAddress("trackE",&spd_trackE);
  T->SetBranchAddress("totEdep",&spd_totEdep);
  T->SetBranchAddress("avg_x",&spd_avg_x);
  T->SetBranchAddress("avg_y",&spd_avg_y);
  T->SetBranchAddress("avg_z",&spd_avg_z);
  T->SetBranchAddress("avg_lx",&spd_avg_lx);
  T->SetBranchAddress("avg_ly",&spd_avg_ly);
  T->SetBranchAddress("avg_lz",&spd_avg_lz);
  T->SetBranchAddress("px",&spd_px);
  T->SetBranchAddress("py",&spd_py);
  T->SetBranchAddress("pz",&spd_pz);
  T->SetBranchAddress("vx",&spd_vx);
  T->SetBranchAddress("vy",&spd_vy);
  T->SetBranchAddress("vz",&spd_vz);
  T->SetBranchAddress("mvx",&spd_mvx);
  T->SetBranchAddress("mvy",&spd_mvy);
  T->SetBranchAddress("mvz",&spd_mvz);
  T->SetBranchAddress("avg_t",&spd_avg_t);

  return;
}
