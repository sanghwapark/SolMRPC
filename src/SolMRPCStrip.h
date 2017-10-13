#ifndef __SOLMRPCSTRIP_H__
#define __SOLMRPCSTRIP_H__


struct Strip_t
{
  int lcStripNum;
  double r_center;
  double r_low;
  double r_up;
  double strip_length;

  //center position
  //  double x0;
  //  double y0;
};

typedef std::vector<Strip_t> StripVector;

class SolMRPCStrip
{

 public:

  SolMRPCStrip();
  //  SolMRPCStrip(int ModuleID, int &StripNum);

  ~SolMRPCStrip();
  
  int getStrip(){ return fStripNum; }
  
  //  TVector3 getGlobalPositionBegin(int moudleID, int StripNum);
  //  TVector3 getGlobalPositionEnd(int moduleID, int StripNum);
  int FindModuleID(double x, double y);
  int FindStrip(int moduleID, double x, double y);
  void InitStripBank(StripVector &vstrips);
  void drawStrips(int moduleID);
  void Print();


//  double get_strip_center(int lcStripNum){ return fStrip_dat->r_cent[i]; }
//  double get_strip_low(int lcStripNum){ return fStrip_dat->r_low[i]; }
//  double get_strip_up(int lcStripNum){ return fStrip_dat->r_up[i]; }

 private:

  int fModuleID;
  int fStripNum;

  StripVector fStripVector;

  double fStripGapWidth;
  double fStripWidth;
  double fStripLength;
  double fAngle;
  double fR;
  double fX;
  double fY;
  double fZ;

  //  TVector3 fGlobalPositionBegin;
  //  TVector3 fGlobalPositionEnd;

};

#endif
