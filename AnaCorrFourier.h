#ifndef _AnaCorrFourier_HH_
#define _AnaCorrFourier_HH_

#include "Track.h"
#include "Event.h"

#include <vector>
#include <iostream>
#include <string>

//#ifndef __CINT__
// #include "boost/smart_ptr.hpp"
// typedef boost::shared_ptr<Event> EVENT_PTR;
// #else
typedef Event* EVENT_PTR;
// #endif

using namespace std;

class TChain;
class TH1;
class TH2;

class AnaCorrFourier{

 public:
  AnaCorrFourier();

  string filelist;
  string outfile;
  void run();
  void SetTotal(int tot);
  void SetVerbosity(int v);
  void SetZbinSize(float zsize);

  int from, to;
 private:

  TChain *track;


  void InitHistos();
  void SaveHistos();
  bool Fillforeground( EVENT_PTR event);
  void Fill(Event* event1, Event* event2, int mixtype);
  int index(int c,int z);

  int get_centbin(float cent);
  int get_centPool(float cent);
  int get_zPool(float z);
  int get_ptBin1(float pt);
  int get_ptBin2(float pt);
  int get_centcut(float cent);

  int get_ptBinSingle(float pt);
  int get_centBinSingle(float pt);
  enum {
    NCENT = 20,
    NPT1 = 5,
    NPT2 = 27,
    NCENTSINGLE = 19,
    NPTSINGLE = 4,

    MAX  = 20000,//ntrk in one event
    NZMON = 40,  //5cm binning for Idiff histograms zvtx
    //pool bins:
    nc = 20,
    nz = 200,
    plnbr = nc*nz
  };
  double depth[nc];
  vector< EVENT_PTR > pool[plnbr];
  int nmix;
  int nevents;
  int total;
  int verbosity;

  TH2* acc[NCENTSINGLE][NPTSINGLE][2];

  TH2* fg[NCENT][NPT1][NPT2][2];
  TH2* bg[NCENT][NPT1][NPT2][2];

  TH1* hfgpt[NCENT];
  TH2* hd0p[NCENT];
  TH2* hz0sinp[NCENT];
  TH2* hd0n[NCENT];
  TH2* hz0sinn[NCENT];

  TH2* hd0   [NCENT];
  TH2* hz0sin[NCENT];

  //global information
  TH1* Idiff[NZMON];
  TH1* hEt;
  TH1* hcent;
  TH1* hzvtx;
  TH1* hxvtx;
  TH1* hyvtx;

  float cent_i;
  float zvtx_i;
  float Et;
  float Zbin_Size;
};

#endif
