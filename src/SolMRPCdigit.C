#include "SolMRPCdigit.h"
#include "SolMRPCStrip.h"
//#include "SolidConst.h"

#include <iostream>

#include "TRandom3.h"
#include "TVector3.h"

using namespace std;

const double DEG=180./3.1415926;   //rad to degree	

const double glass_width = 0.7; //mm
const double gas_width = 0.25; //mm
const double mylar_width = 0.15; //mm
const int ngap = 5;

const double e0 = 1.602177e-19;
const double CtopC = 1.e12;
const double CtofC = 1.e15;

const double RMax = 2100.0; //mm
const double RMin = 960.0; //mm

SolMRPCdigit::SolMRPCdigit(std::string name, string fname) : utlSolidBase(name)
{

  E_weight = calEweight();

  fin = fname;

  _mode = 0;
  _nsteps = 200;
  _stepsize = gas_width/_nsteps;

  // hard copy for now (E = 108 kV/cm for the test beam)
  alpha = 129;
  eta = 5.435;
  drift_v = 0.2012;
  emin = 20.e-6; //BESIII reference 20ev.. (unit of a edep branch is MeV) 

  Qthreshold = 15.e-3; //in pC

}


void SolMRPCdigit::avalanche(int gap, double lz, vector<avlc_pair>& v_sig, double& ind_chg)
{

  ind_chg = 0;
  v_sig.clear();

  double x = 0;
  double n = 1;
  double dx = _stepsize;
  double dt = dx / drift_v;
  double k = eta / alpha;
  double t = 0;

  TRandom3* rnd = new TRandom3();
  TRandom3* rand = new TRandom3();
  rnd->SetSeed(0);
  rand->SetSeed(0);

  double hit_x;
  //count 1-10
  if(gap < 6)
    hit_x = fabs(lz - 0.125);
  else
    hit_x = lz + 0.125;

  //apply effective model when n>200
  bool is_saturated = false;

  double signal = 0;

  for(int i=0; i<_nsteps; i++)
    {

      //within the gap
      if((hit_x + x) >= gas_width) break;

      double tot = 0;

      x = (i+1)*dx;
      t = (i+1)*dt; // time passed after dx

      double nbar = exp((alpha - eta)*dx);

      if(n > 1.5e7) is_saturated = true;

      if(!is_saturated && n > 200)
	{

	  double mean = n * exp((alpha - eta)*dx);
	  double sigma = sqrt(n) * getSigma(dx);

	  double n2 = rand->Gaus(mean, sigma);
	  if(n2 > 1.5e7) n2 = 1.5e7+1;

	  tot = n2;
	}
      else if(!is_saturated && n<=200)
	{
	  for(int j=0; j<n; j++)
	    {
	      double s = rnd->Uniform(1);
	      double det = k* ((nbar -1 ) / (nbar - k));
	      double n1 = 0;
	      if(s < det) n1 = 0;
	      if(s > det) n1 = 1 + (int)(prob(dx,s));

	      tot += n1;
	    }
	}
      else if(is_saturated)
	{
	  tot = 1.5e7+1;
	}

      n = tot;

      double current = n * E_weight * e0 * drift_v * 1.e9; //A
      double fast_charge = current * dt * 1.e-9; //C
      signal = signal + fast_charge;

      avlc_pair av_pair;
      av_pair.istep = i;
      av_pair.nelec = n;
      av_pair.time = t;

      v_sig.push_back(av_pair);

    }

  //  nt.nElectrons = n;
  ind_chg = signal * CtopC;

  delete rnd;
  delete rand;

  return;

}

double SolMRPCdigit::getSigma(double x)
{
  double k = eta / alpha;
  double nbar = exp((alpha - eta)*x);

  double a1 = (1+k)/(1-k);
  double sig2 = sqrt(a1 * nbar * (nbar -1));

  return sig2;
}

double SolMRPCdigit::prob(double x, double s)
{
  double nbar = exp((alpha - eta)*x);
  double k = eta / alpha;

  double a1 = 1 - ((1-k)/(nbar-k));
  double a2 = (nbar-k)*(1-s);
  double a3 = nbar * (1-k);
  double a4 = a2 / a3;

  double f = 1.0/(log(a1)) * log(a4);

  return f;
}

void SolMRPCdigit::SetConstant(int mode)
{

//  switch (mode)
//    {
//    case TEST_BEAM:
//      break;
//    case PHYSICS:
//      break;
//    }
  return;
}

SolMRPCdigit::~SolMRPCdigit()
{
  //destructor
}

void SolMRPCdigit::SetOutTree(TTree* tree)
{

  tree->Branch("EvtNumber",    &nt.EvtNumber,   "EvtNumber/I");
  tree->Branch("rate",         &nt.rate,        "rate/D");
  tree->Branch("pid",          &nt.pid,         "pid/I");
  tree->Branch("mpid",         &nt.mpid,        "mpid/I");
  tree->Branch("tid",          &nt.tid,         "tid/I");
  tree->Branch("mtid",         &nt.mtid,        "mtid/I");

  tree->Branch("gen_px",       &nt.gen_px,      "gen_px/D");
  tree->Branch("gen_py",       &nt.gen_py,      "gen_py/D");
  tree->Branch("gen_pz",       &nt.gen_pz,      "gen_pz/D");
  tree->Branch("gen_vx",       &nt.gen_vx,      "gen_vx/D");
  tree->Branch("gen_vy",       &nt.gen_vy,      "gen_vy/D");
  tree->Branch("gen_vz",       &nt.gen_vz,      "gen_vz/D");

  tree->Branch("timeLeading",  &nt.timeLeading, "timeLeading/D");
  tree->Branch("ToT",          &nt.timeOverThr, "ToT/D");
  tree->Branch("timeRef",      &nt.timeRef,     "timeRef/D");
  tree->Branch("totcharge",    &nt.totcharge,   "totcharge/D");	       
  tree->Branch("fastCharge",   nt.fastCharge,   "fastCharge[200]/D");

  tree->Branch("vp_x",         &nt.vp_x,        "vp_x/D");
  tree->Branch("vp_y",         &nt.vp_y,        "vp_y/D");
  tree->Branch("vp_z",         &nt.vp_z,        "vp_z/D");
  tree->Branch("vp_vx",        &nt.vp_vx,       "vp_vx/D");
  tree->Branch("vp_vy",        &nt.vp_vy,       "vp_vy/D");
  tree->Branch("vp_vz",        &nt.vp_vz,       "vp_vz/D");
  tree->Branch("vp_px",        &nt.vp_px,       "vp_px/D");
  tree->Branch("vp_py",        &nt.vp_py,       "vp_py/D");
  tree->Branch("vp_pz",        &nt.vp_pz,       "vp_pz/D");

  tree->Branch("vp_module",    &nt.vp_module,   "vp_module/I");
  tree->Branch("vp_strip",     &nt.vp_strip,    "vp_strip/I");
  tree->Branch("trackE",       &nt.trackE,      "trackE/D");

  tree->Branch("pid_in_gas",   nt.pid_in_gas,   "pid_in_gas[10]/I");
  tree->Branch("tid_in_gas",   nt.tid_in_gas,   "tid_in_gas[10]/I");
  tree->Branch("mpid_in_gas",  nt.mpid_in_gas,  "mpid_in_gas[10]/I");
  tree->Branch("mtid_in_gas",  nt.mtid_in_gas,  "mtid_in_gas[10]/I");
  tree->Branch("E",            nt.E,            "E[10]/D");
  tree->Branch("edep",         nt.edep,         "edep[10]/D");
  tree->Branch("nion",         nt.nion,         "nion[10]/I");
  tree->Branch("x",            nt.x,            "x[10]/D");
  tree->Branch("y",            nt.y,            "y[10]/D");
  tree->Branch("z",            nt.z,            "z[10]/D");
  tree->Branch("px",           nt.px,           "px[10]/D");
  tree->Branch("py",           nt.py,           "py[10]/D");
  tree->Branch("pz",           nt.pz,           "pz[10]/D");

  tree->Branch("charge",       nt.charge,       "charge[10]/D");
  tree->Branch("nAvalanche",   &nt.nAvalanche,  "nAvalanche/I");
  tree->Branch("module",       nt.module,       "module[10]/I");
  tree->Branch("strip",        nt.strip,        "strip[10]/I");

}

void SolMRPCdigit::Init()
{

  fout = new TFile(outfilename.c_str(),"RECREATE");
  fout->cd();
  T = new TTree("T","T");
  SetOutTree(T);

  h2ChargeTime = new TH2F("h2ChargeTime","Time-Charge correlation", 1500, 0, 15, 6000, -1, 5);
  h2ChargeToT = new TH2F("h2ChargeToT","ToT-Charge correlation", 1500, 0, 15, 6000, -1, 5);
  htotCharge = new TH1F("htotCharge","total induced charge per particle",1500,0,15);
  hRate = new TH1F("hRate", "rate", 20000, 90, 290);
  hT0 = new TH1F("hT0", "time at VP", 30000, 20, 50); //ns
  hToT = new TH1F("hToT", "Time over threshold", 10000, -5, 5); 
  hTimeLeading = new TH1F("hTimeLeading", "Leading edge time", 6000, -1, 5);
  hchannel = new TH1F("hchannel", "hchannel", 2000, 0, 2000);

}

void SolMRPCdigit::End()
{

  fout->cd();
  T->Write();
  htotCharge->Write();
  hRate->Write();
  hT0->Write();
  hToT->Write();
  hTimeLeading->Write();
  h2ChargeTime->Write();
  h2ChargeToT->Write();
  hchannel->Write();

  fout->Close();

  return;

}

void SolMRPCdigit::InitStruct()
{

  nt.rate        = -1.0;
  nt.nAvalanche  = -1.0;
  nt.EvtNumber   = -1.0;
  nt.timeLeading = -1.0;
  nt.timeOverThr = -1.0;
  nt.timeRef     = -1.0;
  nt.totcharge   = -1.0;
  nt.pid         = -999;
  nt.mpid        = -999;
  nt.tid         = -999;
  nt.mtid        = -999;
  nt.trackE      = -999;

  nt.gen_px      = -9.e9;
  nt.gen_py      = -9.e9;
  nt.gen_pz      = -9.e9;
  nt.gen_vx      = -9.e9;
  nt.gen_vy      = -9.e9;
  nt.gen_vz      = -9.e9;

  nt.vp_px        = -9.e9;
  nt.vp_py        = -9.e9;
  nt.vp_pz        = -9.e9;
  nt.vp_vx        = -9.e9;
  nt.vp_vy        = -9.e9;
  nt.vp_vz        = -9.e9;

  nt.vp_x        = -9.e9;
  nt.vp_y        = -9.e9;
  nt.vp_z        = -9.e9;

  nt.vp_module   = -1;
  nt.vp_strip    = -1;

  for(int i=0; i<10; i++)
    {
      nt.pid_in_gas[i]  = -999;
      nt.mpid_in_gas[i] = -999;
      nt.tid_in_gas[i]  = -999;
      nt.mtid_in_gas[i] = -999;

      nt.E[i] = 0;
      nt.edep[i] = 0;
      nt.nion[i] = 0;

      nt.module[i] = -1;
      nt.strip[i] = -1;

      nt.charge[i] = -1;

      nt.x[i] = -9.e9;
      nt.y[i] = -9.e9;
      nt.z[i] = -9.e9;

      nt.px[i] = -9.e9;
      nt.py[i] = -9.e9;
      nt.pz[i] = -9.e9;

    }

  for(int i=0; i<200; i++)
    nt.fastCharge[i] = 0;

}

double SolMRPCdigit::calEweight()
{

  double e_weight = 0;

  //  double epsilon_mylar = 3.1;
  double epsilon_glass = 9; //5-10.. need to check the correct value with Yi

  //  e_weight = 1./(5*0.25 + ((5-1)*0.7 + 2*0.7)/epsilon_glass +2*0.15/epsilon_mylar);
  e_weight = 1./(5*gas_width + ((5-1)*0.7 + 2*0.7)/epsilon_glass);

  return e_weight;

}


int SolMRPCdigit::process(int nevt, int ngen)
{

  SolMRPCStrip* MRPCStrip = new SolMRPCStrip();

  Init();

  double rate = 1.;

  switch (_mode)
    {
    case BEAM:
      rate = 15.e-6/1.6e-19/ngen;
      break;
    case GENERATOR:
      rate = 1.0;
      break;
    case COSMIC:
      rate = 1.0;
      break;
    case SINGLE:
      rate = 1.0;
      break;
    }

  cout << E_weight << "/mm" << endl;

  TFile* file = new TFile(fin.c_str());
  if(!file->IsOpen())
    {
      std::cout << "File is not open" << std::endl;
      return 1;
    }

  file->cd();

  InitContainers();
  TTree* T0 = (TTree*)file->Get("generated");
  SetGenTreeVars(T0);
  TTree* T1 =  (TTree*)file->Get("flux");
  SetFluxTreeVars(T1);
  TTree* T2 = (TTree*)file->Get("solid_mrpc");
  SetMRPCTreeVars(T2);
  TTree* T3 = (TTree*)file->Get("header");
  SetHeaderTreeVars(T3);

  cout << "start loop .." << endl;
  Long64_t nentries = T0->GetEntries();
  int np;
  if(nevt == 0) np = nentries;
  else np = nevt;

  for( Long64_t ientry = 0; ientry<np; ientry++)
    {
      if(ientry%10000==0) cout << ientry << " of " << nentries << endl;
      
      ClearContainers(GENERATED);
      ClearContainers(MRPC);
      ClearContainers(FLUX);
      ClearContainers(HEADER);
      
      T0->GetEntry(ientry); 
      T2->GetEntry(ientry); 
      T1->GetEntry(ientry); //flux tree 
      T3->GetEntry(ientry); 
      
      if(_mode == GENERATOR)
	rate = var8->at(0);

      double evnNumber = evn->at(0);

      if((int)gen_pid->size() > 1)
	cout << "More than one particle generated?" << endl;
      
      double g_px = gen_px->at(0);
      double g_py = gen_py->at(0);
      double g_pz = gen_pz->at(0);
      double g_vx = gen_vx->at(0);
      double g_vy = gen_vy->at(0);
      double g_vz = gen_vz->at(0);

      for(int iflux=0; iflux<(int)flux_hitn->size(); iflux++)
	{
	  
	  InitStruct();

	  int fNAvalanche = 0;

	  int vp_id = (int)flux_id->at(iflux);
	  int vp_pid = (int)flux_pid->at(iflux);
	  int vp_mpid = (int)flux_mpid->at(iflux);
	  int vp_tid = (int)flux_tid->at(iflux);
	  int vp_mtid = (int)flux_mtid->at(iflux);

	  if(vp_id != 4110000) continue;

	  // if(vp_tid != 1) continue; // beam particle only

	  double vp_x = flux_avg_x->at(iflux);
	  double vp_y = flux_avg_y->at(iflux);
	  double vp_z = flux_avg_z->at(iflux);
	  double vp_r = sqrt(pow(vp_x,2) + pow(vp_y,2));

	  if(vp_r > RMax || vp_r < RMin) continue;

	  hRate->Fill(vp_r, rate);

	  double tref = flux_avg_t->at(iflux);
	  hT0->Fill(tref);

	  double flux_E = flux_trackE->at(iflux);
	  nt.trackE = flux_E;

	  //unphysical energy.. something is wrong.
	  if(flux_E > 11000) continue;

	  double vp_vx = flux_vx->at(iflux);
	  double vp_vy = flux_vy->at(iflux);
	  double vp_vz = flux_vz->at(iflux);

	  double vp_px = flux_px->at(iflux);
	  double vp_py = flux_py->at(iflux);
	  double vp_pz = flux_pz->at(iflux);

	  nt.vp_vx = vp_vx;
	  nt.vp_vy = vp_vy;
	  nt.vp_vz = vp_vz;

	  nt.vp_px = vp_px;
	  nt.vp_py = vp_py;
	  nt.vp_pz = vp_pz;

	  nt.EvtNumber = evnNumber;
	  nt.rate = rate;
	  nt.gen_px = g_px;
	  nt.gen_py = g_py;
	  nt.gen_pz = g_pz;
	  nt.gen_vx = g_vx;
	  nt.gen_vy = g_vy;
	  nt.gen_vz = g_vz;

	  nt.pid = vp_pid;
	  nt.tid = vp_tid;
	  nt.mpid = vp_mpid;
	  nt.mtid = vp_mtid;
	  nt.timeRef = tref;
	  nt.vp_x = vp_x;
	  nt.vp_y = vp_y;
	  nt.vp_z = vp_z;
	  
	  int moduleID = -1; 
	  int GlobalStripNum = -1;
	  
	  moduleID = MRPCStrip->FindModuleID(vp_x, vp_y);
	  GlobalStripNum = MRPCStrip->FindStrip(moduleID, vp_x, vp_y);

	  nt.vp_module = moduleID;	
	  nt.vp_strip  = GlobalStripNum;

	  double SumChargeUp[200];
	  double SumChargeLow[200];
	  for(int iarray=0; iarray<200; iarray++)
	    {
	      SumChargeUp[iarray] = 0;
	      SumChargeLow[iarray] = 0;
	    }

	  double tot_charge = 0;
	  double ch_charge[10];
	  for(int i=0; i<10; i++)
	    ch_charge[i] = 0;

	  int sum_nions = 0;
	  for(int i=0; i<(int)mrpc_hitn->size(); i++)
	    {
	      
	      int detector_ID=(int)mrpc_id->at(i)/1000000;
	      int subdetector_ID=((int)mrpc_id->at(i)%1000000)/100000;
	      int subsubdetector_ID=(((int)mrpc_id->at(i)%1000000)%100000)/10000;
	      int component_ID=(int)mrpc_id->at(i)%10000;
	      
	      if (detector_ID==4 && subdetector_ID == 1 && subsubdetector_ID == 0){//in gas     
		
		if(component_ID<1 || component_ID>10){
		  cout<<"MRPC index is wrong "<<mrpc_id->at(i)<<endl;
		  continue;
		}
		
		int pid_in_gas = mrpc_pid->at(i);
		int mpid_in_gas = mrpc_mpid->at(i);
		int tid_in_gas = mrpc_tid->at(i);
		int mtid_in_gas = mrpc_mtid->at(i);

		double x = mrpc_avg_x->at(i);
		double y = mrpc_avg_y->at(i);
		double z = mrpc_avg_z->at(i);

		if( tid_in_gas != vp_tid && mtid_in_gas != vp_tid ) continue;
		if( pid_in_gas != vp_pid && mpid_in_gas != vp_pid ) continue;

		nt.pid_in_gas[component_ID-1]  = pid_in_gas;
		nt.tid_in_gas[component_ID-1]  = tid_in_gas;
		nt.mpid_in_gas[component_ID-1] = mpid_in_gas;
		nt.mtid_in_gas[component_ID-1] = mtid_in_gas;

		nt.x[component_ID-1] = x;
		nt.y[component_ID-1] = y;
		nt.z[component_ID-1] = z;

		//if(abs(vp_pid) == 211)
		// if(mpid_in_gas != 0) continue;

		TVector3 dir_vec;
		dir_vec.SetX(x - vp_x);
		dir_vec.SetY(y - vp_y);
		dir_vec.SetZ(z - vp_z);

		moduleID = MRPCStrip->FindModuleID(x, y);
		GlobalStripNum = MRPCStrip->FindStrip(moduleID, x, y);

		nt.module[component_ID-1] = moduleID;
		nt.strip[component_ID-1] = GlobalStripNum;

		double lz = mrpc_avg_lz->at(i); //z position in local ref system
		
		double px = mrpc_px->at(i);
		double py = mrpc_py->at(i);
		double pz = mrpc_pz->at(i);

		nt.px[component_ID-1] = px;
		nt.py[component_ID-1] = py;
		nt.pz[component_ID-1] = pz;

		double E = mrpc_trackE->at(i);
		nt.E[component_ID-1] = E;

		double edep = mrpc_totEdep->at(i);

		int fNion = (int)(edep / emin);

		if(fNion < 1) continue;

		TRandom3 rand(0);
		int nion = rand.Poisson(fNion);

		if(nion<0) cout << "wrong number of ion .. negative value" << endl;

		nt.edep[component_ID-1] += edep;
		nt.nion[component_ID-1] += nion;

		sum_nions += nion;

		if(nion < 1) continue;

		// if(nion > 200) continue; //too many ions?.. dump the event

		int StripNum = -1;

		for(int j=0; j<nion; j++)
		  {
		    
		    vector<avlc_pair> v_avlc;

		    double induced_charge = 0; //pC
		    
		    double Zrand = rand.Uniform(0,1);

		    //z - lz: center z position of the gas gap
		    //z - lz - (gas_width/2): global z0 position of the gap
		    double z2 = z - lz - (gas_width/2) + (Zrand * gas_width);

		    double ionX = (dir_vec(0)/dir_vec(2)) * (z2 - vp_z) + vp_x;
		    double ionY = (dir_vec(1)/dir_vec(2)) * (z2 - vp_z) + vp_y;
		    StripNum = MRPCStrip->FindStrip(moduleID, ionX, ionY);
		    
		    if(StripNum != -1)
		      hchannel->Fill(StripNum + moduleID*33);
		    else
		      hchannel->Fill(1999);

		    double lz_ion = gas_width*Zrand - (gas_width/2);
		    
		    avalanche(component_ID, lz_ion, v_avlc, induced_charge);

		    if((int)v_avlc.size() < 1) continue;
		    if((int)v_avlc.size() > 200){ cout << "vector size exceeds the limit" << endl; continue;}
		    fNAvalanche++;

		    ch_charge[component_ID-1] = ch_charge[component_ID-1] + induced_charge;
		    tot_charge = tot_charge + induced_charge;

		    for(int k=0; k<(int)(v_avlc.size()); k++)
		      {
			
			double emult = v_avlc[k].nelec;
			double charge = emult * E_weight * e0 * drift_v * (gas_width/200.0/drift_v) * CtopC; //pC

			if(component_ID >5)
			  SumChargeLow[k] = SumChargeLow[k] + charge;
			else
			  SumChargeUp[k] = SumChargeUp[k] + charge;			  
		      }

		  }//nion
	      }//in gas

	    }//mrpc hit loop

	  for(int i=0; i<10; i++)
	    nt.charge[i] = ch_charge[i];

	  htotCharge->Fill(tot_charge, rate);

	  bool IsOverThreshold = false;
	  double time_leading = 0;
	  int tOverThreshold = 0;
	  for(int k=0; k<200; k++)
	    {
	      double SumCharge = SumChargeLow[k] + SumChargeUp[k];
	      nt.fastCharge[k] = SumCharge;

	      if( SumCharge > Qthreshold )
		{
		  if(!IsOverThreshold)
		    {
		      IsOverThreshold = true;
		      time_leading = (gas_width/200.0/drift_v) * (k+1);
		    }
		  tOverThreshold++;
		}
	    }

	  double ToT = (gas_width/200.0/drift_v) * tOverThreshold;

	  nt.nAvalanche = fNAvalanche;
	  nt.totcharge = tot_charge;		
	  nt.timeLeading = time_leading;
	  nt.timeOverThr = ToT;

	  T->Fill();

	  hTimeLeading->Fill(time_leading, rate);
	  hToT->Fill(ToT, rate);
	  h2ChargeTime->Fill(tot_charge, time_leading, rate);
	  h2ChargeToT->Fill(tot_charge, ToT, rate);

	}//flux loop	       

    }// event loop

  delete file;
  delete MRPCStrip;

  return 0;

}


/*
  BESIII references
DO NOT REMOVE LINES BELOW
*/

/*
double SolMRPCdigit::getEta(double E)
{
  
  double alpha[22]=
    {
      383.5/10  ,
      471  /10  ,
      564.5/10  ,
      663.6/10  ,
      777.1/10  ,
      877  /10  ,
      990.8/10  ,
      1106 /10  ,
      1154 /10  ,
      1199 /10  ,
      1253 /10  ,
      1296 /10  ,
      1344 /10  ,
      1396 /10  ,
      1448 /10  ,
      1502 /10  ,
      1545 /10  ,
      1597 /10  ,
      1726 /10  ,
      1858 /10  ,
      1992 /10  ,
      2124 /10  ,
    };

  return alpha[E];
}

double SolMRPCdigit::getAlpha(double E)
{
  
  double alpha[22]=
    {
      383.5/10  ,
      471  /10  ,
      564.5/10  ,
      663.6/10  ,
      777.1/10  ,
      877  /10  ,
      990.8/10  ,
      1106 /10  ,
      1154 /10  ,
      1199 /10  ,
      1253 /10  ,
      1296 /10  ,
      1344 /10  ,
      1396 /10  ,
      1448 /10  ,
      1502 /10  ,
      1545 /10  ,
      1597 /10  ,
      1726 /10  ,
      1858 /10  ,
      1992 /10  ,
      2124 /10  ,
    };

  return alpha[E];
}
*/
