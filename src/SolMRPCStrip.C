#include <iostream> 

#include "TVector2.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"

#include "SolMRPCStrip.h"

#define PI 3.14159265

SolMRPCStrip::SolMRPCStrip()
{

  fStripGapWidth = 3.; //mm
  fStripWidth = 25.; //mm

  InitStripBank(fStripVector);

}

//____________________________________________________

SolMRPCStrip::~SolMRPCStrip()
{

}

//____________________________________________________
void SolMRPCStrip::drawStrips(int moduleID)
{

  bool draw_all = false;
  if(moduleID < 0 || moduleID > 49) draw_all = true;

  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->cd();

  TH2F* h2line = new TH2F("h2line", "Strip map", 4400, -2200, 2200, 4400, -2200, 2200);
  h2line->Draw();

  TLine line;

  for(int imodule = 0; imodule<50; imodule++)
    {
      for(int istrip = 0; istrip<33; istrip++)
	{
	  double r_center = fStripVector[istrip].r_center;
	  double Lstrip = fStripVector[istrip].strip_length;

	  double angle = imodule * 7.2; //degree
	  double angle_in_rad =  angle * PI / 180.;

	  double x1 = r_center; 
	  double x2 = r_center; 
	  double y1 = Lstrip * 0.5;
	  double y2 = Lstrip * 0.5 * -1;

	  TVector2 v1(x1, y1);
	  TVector2 v2(x2, y2);

	  double new_x1 = v1.Rotate(angle_in_rad).X();
	  double new_x2 = v2.Rotate(angle_in_rad).X();
	  double new_y1 = v1.Rotate(angle_in_rad).Y();
	  double new_y2 = v2.Rotate(angle_in_rad).Y();

	  line.SetLineColor(kBlue);
	  if(draw_all)
	    {
	      line.DrawLine(new_x1, new_y1, new_x2, new_y2);
	    }
	  else
	    {
	      if(imodule == moduleID)
		line.DrawLine(new_x1, new_y1, new_x2, new_y2);
	      else
		continue;
	    }
	}
    }

}

//____________________________________________________
int SolMRPCStrip::FindModuleID(double x, double y)
{

  //moduleID 0-49
   int sector = -1;

   double phi = atan2(y,x) * 180 / PI; // -pi to +pi

   double sec_shift= 3.6;  // shift to match electron turning in field

   double hit_phi = (phi<0)?(phi+360):phi;
   double new_phi = 0;


   if(hit_phi+sec_shift > 360)
     {
       new_phi = hit_phi - 360 + sec_shift;
       sector = int(new_phi/7.2);
     }
   else
     {
       new_phi = hit_phi + sec_shift;
       sector = int(new_phi/7.2);
     }

   /*
   //original
   int sec_shift = 0;
   if (hit_phi>=90+sec_shift) sector=int((hit_phi-90-sec_shift)/7.2);
   else sector=int((hit_phi+360-90-sec_shift)/7.2);
   */

   return sector;
}

//____________________________________________________
void SolMRPCStrip::Print()
{

  for( unsigned int it = 0; it < fStripVector.size(); it++)
    {
      std::cout << "Strip: " << it << " " 
		<< fStripVector[it].r_low << " "
		<< fStripVector[it].r_center << " "
		<< fStripVector[it].r_up << std::endl;
    }

}

//____________________________________________________

void SolMRPCStrip::InitStripBank(StripVector &vstrips)
{

  vstrips.clear();

  for(int i=0; i<33; i++)
    {

      //pCDR design: 13-17 cm for the inner module
      double Lstrip = 130. + 4 * i;

      //model #1
      /*
      double r0 = 960;
      double edge = 36.; //mm
      double r_cent = r0 + (edge + 12.5) +  edge * (int(i/11) * 2) + (25 * i) + (3 * i);
      double r_low  = r_cent - 12.5;
      double r_up   = r_cent + 12.5;
      */

      //model #2
      /*
      double r0 = 960;
      double edge = 109.5;
      double r_cent = r0 + edge + 12.5 + (25 * i) + (3 * i);
      double r_low  = r_cent - 12.5;
      double r_up   = r_cent + 12.5;
      */

      //model #3 - 1
      double r0 = 1033.15;
      double r_cent = r0 + 12.5 + (25 * i) + (3 * i);
      double r_low  = r_cent - 12.5;
      double r_up   = r_cent + 12.5;

      //model #3 - 2
      /*
      double r0 = 1050.;
      double r_cent = r0 + 12.5 + (25 * i) + (3 * i);
      double r_low  = r_cent - 12.5;
      double r_up   = r_cent + 12.5;
      */

      /*
      double half_angle = 360./50./2.;
      double Lstrip = 130. + 4 * i;
      double r_cent = (Lstrip/2) / tan(half_angle*PI/180);
      double r_low = r_cent - 12.5;
      double r_up = r_cent + 12.5;
      */

      Strip_t st;
      st.lcStripNum = i;
      st.r_center = r_cent;
      st.r_low = r_low;
      st.r_up = r_up;
      st.strip_length = Lstrip;

      vstrips.push_back(st);

    }

}

//____________________________________________________

int SolMRPCStrip::FindStrip(int moduleID, double x, double y)
{

  int lcNum = -1;

  double r = sqrt(x*x + y*y);

  double center_angle = moduleID * 7.2;

  TVector2 v0(1,0); //start from phi=0, i.e. positive x-axis
  double new_x0 = v0.Rotate(center_angle*PI/180).X();
  double new_y0 = v0.Rotate(center_angle*PI/180).Y();

  v0.Set(new_x0,new_y0);

  TVector2 vr(x,y);

  double dphi = vr.DeltaPhi(v0);

  double r2 = r * cos(fabs(dphi)*PI/180);

  if( r2 - 3 < fStripVector[0].r_low || r2 + 3 > fStripVector[32].r_up ) return -1;

  for( unsigned int it = 0; it < fStripVector.size(); it++)
    {
      if( r2 >= fStripVector[it].r_low && r2 < fStripVector[it].r_up )
	{
	  lcNum = fStripVector[it].lcStripNum;
	  break;
	}
    }

  double dr = 25.;
  if( lcNum == -1 )
    {
      for( unsigned int it = 0; it < fStripVector.size(); it++)
	{
	  double _dr = fabs(r2 - fStripVector[it].r_center);
	  if(_dr < dr) 
	    {
	      dr = _dr;
	      lcNum = fStripVector[it].lcStripNum;
	    }
	}
    }

  return lcNum;

}

