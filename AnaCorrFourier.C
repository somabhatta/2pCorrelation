#include "AnaCorrFourier.h"
#include "string.h"
#include <cmath>
#include <fstream>
#include <memory>
#include <TRandom.h>
#include <TRandom3.h>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2D.h"

#define USE_OLD1
const double PI = acos(-1.0);
static TFile *tmpf;
static TRandom *ran;
 

using namespace std;

Event::Event(){
  Tracks.clear();
}

Event::~Event()
{

}

void  Event::identify(){
}

void
Event::AddTrack(float pt, float eta, float phi0, int qual)
{
  // TRACK_PTR ph = new Track();                                                                                                                        
  Track *ph=new Track();
  ph->set_pt(pt);
  ph->set_eta(eta);
  ph->set_phi0(phi0);
  ph->set_qual(qual);
  ph->set_qop(-9999);

  Tracks.push_back(ph);
}
void
Event::AddTrack(float pt, float eta, float phi0, float qop, int qual)
{
  //  TRACK_PTR ph = new Track();                                                                                                                       
  Track* ph = new Track();
  ph->set_pt(pt);
  ph->set_eta(eta);
  ph->set_phi0(phi0);
  ph->set_qop(qop);
  ph->set_qual(qual);

  Tracks.push_back(ph);
}


Track*
Event::GetTrack(unsigned int i) const

{
  if(i>=Tracks.size()) return 0;
  return Tracks[i];
  //   return Tracks[i].get();                                                                                                                         

}

/* //Our Old private cuts
static double centcuts[]={
  3.44565,3.30005,3.16895,3.04535,2.92845,2.81735,2.71095,2.60965,2.51385,2.42155,//2.44
  2.33225,2.24765,2.16585,2.08535,2.00995,1.93635,1.86495,1.79595,1.72815,1.66365,//1.67
  1.60145,1.54025,1.48175,1.42445,1.36885,1.31535,1.26375,1.21275,1.16435,1.11655,//1.09
  1.07035,1.02615,0.98245,0.94055,0.90045,0.86125,0.82335,0.78645,0.75065,0.71625,//0.7
  0.68285,0.65075,0.61965,0.58955,0.56045,0.53265,0.50625,0.48035,0.45555,0.43145,//0.41
  0.40825,0.38595,0.36495,0.34435,0.32475,0.30605,0.28805,0.27085,0.25455,0.23885,//0.22
  0.22395,0.20955,0.19615,0.18315,0.17085,0.15925,0.14845,0.13815,0.12845,0.11905,//0.11
  0.11025,0.10215,0.09445,0.08725,0.08055,0.07425,0.06825,0.06275,0.05765,0.05275,//0.05
  0.04825,0.04405,0.04015,0.03655,0.03315,0.02995,0.02705,0.02435,0.02175,0.01935,//0.019
  0.01705,0.01505,0.01305,0.01125,0.00955,0.00795,0.00645,0.00485,0.00315,0
};// */


// /* //New 100% efficiency cuts
static double centcuts[]={
  3.44515,3.29855,3.16725,3.04515,  2.93035  ,2.81985,2.71485,2.61315,2.51595,  2.42335,//2.44
  2.33245,2.24735,2.16545,2.08525,  2.00895  ,1.93425,1.86295,1.79255,1.72565,  1.66125,//1.67
  1.60145,1.54025,1.48175,1.42445,  1.36585  ,1.31535,1.26375,1.21275,1.16435,  1.11605,//1.09
  1.07035,1.02615,0.98245,0.94055,  0.90045  ,0.86125,0.82335,0.78645,0.75065,  0.71645,//0.7
  0.68285,0.65075,0.61965,0.58955,  0.56045  ,0.53265,0.50625,0.48035,0.45555,  0.43035,//0.41
  0.40825,0.38595,0.36495,0.34435,  0.32425  ,0.30605,0.28805,0.27085,0.25455,  0.23865,//0.22
  0.22395,0.20955,0.19615,0.18315,  0.17065  ,0.15925,0.14845,0.13815,0.12845,  0.11885,//0.11
  0.11025,0.10215,0.09445,0.08725,  0.08035  ,0.07425,0.06825,0.06275,0.05765,  0.05265,//0.05
  0.04825,0.04405,0.04015,0.03655,  0.03295  ,0.02995,0.02705,0.02435,0.02175,  0.01915,//0.019
  0.01705,0.01505,0.01305,0.01125,  0.00935  ,0.00795,0.00645,0.00485,0.00315,  0
}; // */


 /* //New 102% efficiency cuts
static double centcuts[]={
  3.44805,3.30415,3.17455,3.05455,  2.94135  ,2.83295,2.72885,2.62855,2.53355,  2.44105,//2.44
  2.35135,2.26725,2.18615,2.10655,  2.03185  ,1.95735,1.88615,1.81675,1.75035,  1.68585,//1.67
  1.60145,1.54025,1.48175,1.42445,  1.39225  ,1.31535,1.26375,1.21275,1.16435,  1.14265,//1.09
  1.07035,1.02615,0.98245,0.94055,  0.92815  ,0.86125,0.82335,0.78645,0.75065,  0.74275,//0.7
  0.68285,0.65075,0.61965,0.58955,  0.58605  ,0.53265,0.50625,0.48035,0.45555,  0.45415,//0.41
  0.40825,0.38595,0.36495,0.346  ,  0.34515  ,0.30605,0.28805,0.27085,0.257  ,  0.25685,//0.22
  0.22395,0.20955,0.19615,0.187  ,  0.18615  ,0.15925,0.14845,0.13815,0.132  ,  0.13175,//0.11
  0.11025,0.10215,0.09445,0.091  ,  0.09035  ,0.07425,0.06825,0.06275,0.061  ,  0.06035,//0.05
  0.04825,0.04405,0.04015,0.039  ,  0.03875  ,0.02995,0.02705,0.02435,0.024  ,  0.02345,//0.019
  0.01705,0.01505,0.01305,0.013  ,  0.01265  ,0.00795,0.00645,0.00485,0.003  ,  0
}; // */
  
  
 /* //New 98% efficiency cuts
static double centcuts[]={
  3.44215,3.29315,3.15945,3.03555,  2.91875  ,2.80735,2.70025,2.59805,2.49885,  2.40495,//2.44
  2.31285,2.22725,2.14325,2.06365,  1.98555  ,1.91105,1.83825,1.76835,1.69975,  1.63505,//1.67
  1.60145,1.54025,1.48175,1.42445,  1.33825  ,1.31535,1.26375,1.21275,1.16435,  1.08765,//1.09
  1.07035,1.02615,0.98245,0.94055,  0.87205  ,0.86125,0.82335,0.78645,0.75065,  0.69855,//0.7
  0.68285,0.65075,0.61965,0.58955,  0.53445  ,0.53265,0.50625,0.48035,0.45555,  0.40695,//0.41
  0.406  ,0.38595,0.36495,0.34435,  0.30345  ,0.302  ,0.28805,0.27085,0.25455,  0.22045,//0.22
  0.220  ,0.20955,0.19615,0.18315,  0.15545  ,0.154  ,0.14845,0.13815,0.12845,  0.10655,//0.11
  0.106  ,0.10215,0.09445,0.08725,  0.07085  ,0.070  ,0.06825,0.06275,0.05765,  0.04555,//0.05
  0.045  ,0.04405,0.04015,0.03655,  0.02755  ,0.0271 ,0.02705,0.02435,0.02175,  0.01515 ,//0.019
  0.0151 ,0.01505,0.01305,0.01125,  0.00625  ,0.0061 ,0.006  ,0.00485,0.00315,  0
}; // */





AnaCorrFourier::AnaCorrFourier():
  filelist(""),
  outfile(""),
  track(0),
  nevents(0),
  total(0),
  verbosity(0){ 
  cout << " mixing code" << endl;
  from = 0; to = 0;
  for(int i=0;i<12;i++){
    depth[i] = 2;
  }
  depth[12] = 4;
  depth[13] = 4;
  depth[14] = 8;
  depth[15] = 8;
  depth[16] = 20;
  depth[17] = 20;
  depth[18] = 20;
  depth[19] = 20;
  ran = new TRandom;
  Zbin_Size=1.0;
}

void
AnaCorrFourier::SetZbinSize(float zsize)
{
  if (zsize<1.0) {
    cout<<"Zbin_Size cant be smaller than 1.0:: It is being set to 1"<<endl;
    Zbin_Size=1.0;
  }
  else {
    cout<<"Zbin_Size is being set to "<<zsize<<endl;
    Zbin_Size=zsize;
  }
}

void
AnaCorrFourier::SetTotal(int tot)
{
  total = tot;
}

void
AnaCorrFourier::SetVerbosity(int v)
{
  verbosity = v;
}

void 
AnaCorrFourier::run()
{

  InitHistos();

  if(total!=0) cout << "will run " << total << " events." << endl;

  cout << "Let us start running! " << endl;

  track = new TChain("HeavyIonD3PD");
  cout << track->GetMaxVirtualSize() << endl;
  char fname[400];
  ifstream lis(filelist.c_str());
  int cnt=0;
  while(!lis.eof()){
    string filename;
    lis >> filename;
    sprintf(fname,"%s",filename.c_str());
    if(cnt>=from&&cnt<to){
      cout << fname << endl;
      if(!filename.empty()) track->Add(fname);
    }
    cnt++;
  }

  cout << " track : " << track << endl;

  int   ntrk,cent;
  float xvtx=-5,yvtx=-5, zvtx;

  float pt[MAX],phi0[MAX],eta[MAX],qop[MAX];
  float d0[MAX],z0[MAX],d0_err[MAX],z0_err[MAX],theta_err[MAX];
  int   qual[MAX];
 
  track->SetBranchStatus("*",0);
  track->SetBranchStatus("trk_n",1);
  track->SetBranchStatus("Centrality",1);
  track->SetBranchStatus("Fcal_Et",1);
  track->SetBranchStatus("vx_x",1);
  track->SetBranchStatus("vx_y",1);
  track->SetBranchStatus("vx_z",1);
  track->SetBranchStatus("trk_pt",1);
  track->SetBranchStatus("trk_phi0_wrt_PV",1);
  track->SetBranchStatus("trk_eta",1);
  track->SetBranchStatus("trk_d0_wrtPV",1);
  track->SetBranchStatus("trk_z0_wrtPV",1);
  track->SetBranchStatus("trk_Quality",1);
  track->SetBranchStatus("trk_qop_wrtPV",1);
  track->SetBranchStatus("trk_err_d0",1);
  track->SetBranchStatus("trk_err_z0",1);
  track->SetBranchStatus("trk_err_theta",1);


  track->SetBranchAddress("trk_n",&ntrk);
  track->SetBranchAddress("Centrality",&cent);
  track->SetBranchAddress("Fcal_Et",&Et);
  track->SetBranchAddress("vx_x",&xvtx);
  track->SetBranchAddress("vx_y",&yvtx);
  track->SetBranchAddress("vx_z",&zvtx);

  track->SetBranchAddress("trk_pt",&pt);
  track->SetBranchAddress("trk_phi0_wrt_PV",&phi0);
  track->SetBranchAddress("trk_eta",&eta);
  track->SetBranchAddress("trk_d0_wrtPV",&d0);
  track->SetBranchAddress("trk_z0_wrtPV",&z0);
  track->SetBranchAddress("trk_Quality",&qual);
  track->SetBranchAddress("trk_qop_wrtPV",&qop);

  track->SetBranchAddress("trk_err_d0",&d0_err);
  track->SetBranchAddress("trk_err_z0",&z0_err);
  track->SetBranchAddress("trk_err_theta",&theta_err);
 

  int nentries = track->GetEntries();
  cout << "nentries = "<<nentries << " out of total  =  "<<total<<endl;

  for(int i=0; i<nentries;i++){
    track->GetEntry(i);
//////////////////////////////////////////////////////
//******************************** IMPORTANT*************************************
Et=Et*1.041; //This is because in the november 2011 run the EM scale in the Fcal was off by 4%
////////////////////////////////////



    cent_i = get_centbin(Et);
    if(i%1000==0) cout<<i<<"    "<<Et<<" "<<cent_i<<endl;    
    hcent->Fill(cent_i);
    hEt->Fill(Et);

    int centcb = get_centcut(cent_i);
    if(centcb<0) continue;
    //if(cent_i>95) continue;

    //if(mbtime_countA<9 || mbtime_countC<9) continue;     // MBTS HITS

    hzvtx->Fill(zvtx);
    int zbin   = get_zPool(zvtx);
    if(zbin<0) continue;
    if(cent_i<80){
      hxvtx->Fill(xvtx);
      hyvtx->Fill(yvtx);
    }

    int centbinsingle = get_centBinSingle(cent_i);

    Event* ev = new Event();
    
    ev->set_id(i);
    ev->set_zvtx(zvtx);
    ev->set_cent(cent);
    int ns=0;
    for(int j=0; j<ntrk; j++){
      pt[j]/=1000;

      if(qual[j]!=1) continue;  
      double z0sin = z0[j]*(2*exp(-eta[j])/(1+exp(-2*eta[j])));
      //if(fabs(d0[j])>1||fabs(z0sin)>1) continue; // already applied by qual[j]==1 cut

      double d0cut = sqrt(0.0738*0.0738+0.637*0.637/pt[j]/pt[j]);
      double z0cut = sqrt(0.19*0.19+0.882*0.882/pt[j]/pt[j]);
      if(fabs(d0[j])<d0cut){
	if(qop[j]>0) hz0sinp[centcb]->Fill(pt[j],z0sin);
	else  hz0sinn[centcb]->Fill(pt[j],z0sin);
      }
      if(fabs(z0sin)<z0cut){
        if(qop[j]>0) hd0p[centcb]->Fill(pt[j],d0[j]);
	else  hd0n[centcb]->Fill(pt[j],d0[j]);
      }


//---------------------------------------------------------------------------
      double sin_theta=2*exp(-eta[j])/(1+exp(-2*eta[j]));
      double cos_theta=sqrt(1.0-sin_theta*sin_theta); 
      double e_z0= pow(z0_err[j]*sin_theta,2) +  pow(theta_err[j]*z0[j]*cos_theta,2);
      e_z0 =(e_z0>0)? sqrt(e_z0):0;
      double e_d0 = d0_err[j]; 
      if(e_d0!=0 && e_z0!=0) {
        if(fabs(d0[j]/e_d0)<3.0) hz0sin[centcb]->Fill(pt[j],z0sin/e_z0);
        if(fabs(z0sin/e_z0)<3.0) hd0   [centcb]->Fill(pt[j],d0[j]/e_d0);
      } 

 //Our Default tracking cuts
//      if(fabs(d0[j])>d0cut  || fabs(z0sin)>z0cut) continue; // */


/* //Atillo's tracking cuts1
      if(fabs(d0[j])>0.5 || fabs(z0sin) >0.5)      continue; // */

/* //Atillo's tracking cuts2
      if(e_d0==0 || e_z0==0 ) continue;
      if(fabs(z0sin/e_z0) >3) continue;  
      if(fabs(d0[j]/e_d0) >3) continue;  // */
//------------------------------------------------------------------------

      int ptbinsingle = get_ptBinSingle(pt[j]);
      if(ptbinsingle>-1){
	if(qop[j]>0) acc[centbinsingle][ptbinsingle][0]->Fill(phi0[j],eta[j]);
	else         acc[centbinsingle][ptbinsingle][1]->Fill(phi0[j],eta[j]);
      }
      ev->AddTrack(pt[j],eta[j],phi0[j],qop[j],1);
      hfgpt[centcb]->Fill(pt[j]);//raw spectra
      ns++;
    }
    ev->set_npart(ns);


    if(nevents%1000==0) cout << "Has run "<< nevents << " events" << endl;
    if(nevents%3000==0&&nevents>0){
      int nz1=nz*1.0/Zbin_Size;
      if (nz1>nz) {nz1=nz; cout<<"ERROR nz1 > nz"<<endl;}
      for(int cbin=0;cbin<17;cbin+=3){
	for(int zbin=0;zbin<nz1;zbin++){
	  int indx = index(cbin,zbin); 
	  int size = pool[indx].size();
	  if(size<depth[0]){
	    cout<<"cbin: "<<cbin<<" "<<"zbin: "<<zbin<<" "<< pool[indx].size()<<endl;
	  }
	}
      }
    }


    if(total!=0 && nevents>=total) break;
    nevents++;
    
    Fillforeground(ev);

  }//entry
  SaveHistos();
}

bool AnaCorrFourier::Fillforeground( Event* event)
{
  float zvtx = event->get_zvtx();
  int zbin   = get_zPool(zvtx);
  int centbin = get_centPool(cent_i);
  int indx = index(centbin,zbin); 
  int dep = depth[centbin];
  if(zbin==-1||centbin==-1) return false;

  // Fill( event.get(), event.get(), 0);//0: forground;
  Fill( event, event, 0);//0: forground; 
    int zbinmon=(zbin*NZMON)/(nz+0.0001);

    // int id0=event.get()->get_id();
    int id0=event->get_id(); 
    int KiddiePool = pool[indx].size();
    nmix = KiddiePool;
    if(KiddiePool>0 && KiddiePool<=dep){
      for(int k=0; k<KiddiePool; k++){
        int id1=(pool[indx][k])->get_id();
	//Fill(event.get(), (pool[indx][k]).get(), 1);//;mixing
	Fill(event, (pool[indx][k]), 1);//;mixing
        Idiff[zbinmon]->Fill(id0-id1);//ID Difference 
      }
    }
    
    if(KiddiePool>=dep){
      pool[indx].erase( (pool[indx]).begin() );
    }
    pool[indx].push_back( event );
  return true;
}


void
AnaCorrFourier::Fill(Event* event1, Event *event2, int mixtype)
{
  
  int icent = get_centcut(cent_i);
  if(icent<0) return;

  int ntrk1 = event1->get_npart();
  int ntrk2 = event2->get_npart();

  float pt1,phi01,eta1,qop1;
  float pt2,phi02,eta2,qop2;

  int sig1 =(int)(ran->Rndm()*2);
  
  float wei = 1.0/nmix;
  
  for( int i=0; i<ntrk1; i++){
     Track* trk1 = event1->GetTrack(i);
    eta1  = (trk1)->get_eta();
    pt1   = (trk1)->get_pt();
    phi01 = (trk1)->get_phi0();
    qop1  = (trk1)->get_qop();
    
     int ptbin1 = get_ptBin1(pt1);//trigger bin
    if(ptbin1<0) continue;
     int j=0;
    while(j<ntrk2)  {	
      if(mixtype==0&&j==i) {j++;continue;}
      Track* trk2 = event2->GetTrack(j);
      j++;
       pt2   = (trk2)->get_pt();	
       phi02 = (trk2)->get_phi0();
       eta2  = (trk2)->get_eta();
       qop2  = (trk2)->get_qop();
       int ptbin2 = get_ptBin2(pt2);
       if(ptbin2<0) continue;
      
       float dphi0 = phi01-phi02;
       float deta  = fabs(eta1-eta2);
       if(sig1) dphi0 = -dphi0;
       if     (dphi0>1.5*PI)  dphi0 -=2*PI;
        else if(dphi0<-0.5*PI) dphi0 +=2*PI;	
       // cout<<"del phi = "<<dphi0<<endl;
       // cout<<"del eta = "<<deta<<endl;
            if(mixtype==0){
		if((qop1*qop2)>0) fg[icent][ptbin1][ptbin2][0]->Fill(dphi0,deta);
		else              fg[icent][ptbin1][ptbin2][1]->Fill(dphi0,deta);
	      }else if(mixtype==1){
		if((qop1*qop2)>0) bg[icent][ptbin1][ptbin2][0]->Fill(dphi0,deta,wei);
		else              bg[icent][ptbin1][ptbin2][1]->Fill(dphi0,deta,wei);
       }
         if(cent_i<1.01){//top 1%
		if(mixtype==0){
	    if((qop1*qop2)>0) fg[19][ptbin1][ptbin2][0]->Fill(dphi0,deta);
	  	  else              fg[19][ptbin1][ptbin2][1]->Fill(dphi0,deta);
      	}else if(mixtype==1){
      	  if((qop1*qop2)>0) bg[19][ptbin1][ptbin2][0]->Fill(dphi0,deta,wei);
      	  else              bg[19][ptbin1][ptbin2][1]->Fill(dphi0,deta,wei);
         }
       }
}
  }//while 
}

int 
AnaCorrFourier::get_ptBinSingle(float lpt){//8 bins in pt
  int bin=-1;
  if(lpt<0.5||lpt>=3) return -1;

  //total 3 bins
  if(lpt<1) bin=0;
  else if(lpt<1.5) bin =1;
  else if(lpt<2.0) bin =2;
  else  bin =3;
  return bin; 
}
int AnaCorrFourier::get_centBinSingle(float cent){//6 bins
  int bin=-1;
  if(cent<0 || cent > 95) return -1;
  bin = cent/5.01;
  return bin;
}

int 
AnaCorrFourier::get_ptBin1(float lpt){
  int bin=-1;
  if(lpt<0.5||lpt>=4) return -1;

  //total 5 bins
  if(lpt<1) bin=0;
  else if(lpt<1.5) bin =1;
  else if(lpt<2.0) bin =2;
  else if(lpt<3.0) bin =3;
  else  bin =4;

  return bin; 
}

int AnaCorrFourier::get_ptBin2(float lpt){
  int bin=-1;
  if(lpt<0.5||lpt>20) return -1;

  //total 26 bins +1 extra bin for the 1.5-1.6 bin
  if(lpt<1.5){//0-17 // now it is 0-5  
    bin = (lpt-0.4)/0.2;
  }else if(lpt<1.6){// the extra 1.5-1.6 bin
    bin = 6;   
  }else if(lpt<4){ // added 1 from here on
    bin = (lpt-0.4)/0.2 +1;
  }else if(lpt<5){//18,19 +1
    bin = (lpt-4)/0.5+19;
  }else if(lpt<8){//20-22 +1
    bin = (lpt-5)+21;
  }else if(lpt<12){//23,24 +1
    bin = (lpt-8)/2+24;
  }else{
    bin = 26; //25 +1
  }
  return bin;

}

int 
AnaCorrFourier::index(int c,int z)
{
  return  z + c*nz;
}


int AnaCorrFourier::get_centbin(float Et){
  int cbin=100;
  for(int i=0;i<100;i++){
    if(Et>centcuts[i]){
      cbin = i+1; break;//centrality bin starts with 1
    }
  }
  return cbin;
}

int AnaCorrFourier::get_centPool(float cent){
  if(cent<0 || cent > 95) return -1;
  int bin=-1;
  bin = cent/5.01;
  return bin;
}
int AnaCorrFourier::get_zPool(float z){
  int bin=-1;
  if(fabs(z)>=100) return -1;
  bin = (z+100.0)/Zbin_Size;  // pools 
  if(bin<0||bin>=nz) return -1;
  return bin;
}

int AnaCorrFourier::get_centcut(float cent){//19 bins here, 1% bin filled seperatly.
  if(cent<0 || cent > 95) return -1;
  int bin=-1;
  bin = cent/5.01;
  return bin;
}

void
AnaCorrFourier::InitHistos()
{
  char histname[100];

  tmpf = new TFile(outfile.c_str(),"recreate");
  double ptbins[]={0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,12,14,16,18,20,25,30,40,50};


  for(int izbin=0;izbin<NZMON;izbin++)
  {
    sprintf(histname,"Idiff%d",izbin);
    Idiff[izbin]=new TH1D(histname,histname,10000,-0.5,199999.5);  //Max difference is 200,000
  }


  for(int icent=0; icent<NCENT; icent++){
    sprintf(histname,"hd0p_cent%d",icent);
    hd0p[icent] = new TH2D(histname,histname,31,ptbins,1000,-0.5,0.5);
    sprintf(histname,"hz0sinp_cent%d",icent);
    hz0sinp[icent] = new TH2D(histname,histname,31,ptbins,1000,-1.0,1.0);

    sprintf(histname,"hd0n_cent%d",icent);
    hd0n[icent] = new TH2D(histname,histname,31,ptbins,1000,-0.5,0.5);
    sprintf(histname,"hz0sinn_cent%d",icent);
    hz0sinn[icent] = new TH2D(histname,histname,31,ptbins,1000,-1.0,1.0);

    sprintf(histname,"hz0sin_cent%d",icent);
    hz0sin[icent]=new TH2D(histname,"hz0sin;pt;z0/z0err",31,ptbins,2000,-10,10);

    sprintf(histname,"hd0_cent%d",icent);
    hd0   [icent]=new TH2D(histname,"hd0;pt;d0/d0err"   ,31,ptbins,2000,-10,10);
  }

  for(int icent=0; icent<NCENT; icent++){
    sprintf(histname,"pt_cent%d",icent);
    hfgpt[icent] = new TH1D(histname,histname,150,0,30);
    hfgpt[icent]->Sumw2();
  }
  for(int icent=0; icent<NCENTSINGLE; icent++){
    for(int ipt=0; ipt<NPTSINGLE; ipt++){
      for(int it=0; it<2; it++){
	sprintf(histname,"acc_cent%d_pt%d_ch%d",icent,ipt,it);
	acc[icent][ipt][it] = new TH2D(histname,histname,100,-PI,PI,100,-2.5,2.5);
	acc[icent][ipt][it]->Sumw2();
      }
    }
  }
  int netabin = 18;
  double etabins[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
  for(int icent=0; icent<NCENT; icent++){
    for(int ipt1=0; ipt1<NPT1; ipt1++){
      for(int ipt2=0; ipt2<NPT2; ipt2++){
	for(int it=0; it<2; it++){
	  sprintf(histname,"fg_cent%d_pta%d_ptb%d_ch%d",icent,ipt1,ipt2,it);
	  fg[icent][ipt1][ipt2][it] = new TH2D(histname,histname,200,-PI/2,1.5*PI,netabin,etabins);
	  fg[icent][ipt1][ipt2][it]->Sumw2();
	  sprintf(histname,"bg_cent%d_pta%d_ptb%d_ch%d",icent,ipt1,ipt2,it);
	  bg[icent][ipt1][ipt2][it] = new TH2D(histname,histname,200,-PI/2,1.5*PI,netabin,etabins);
	  bg[icent][ipt1][ipt2][it]->Sumw2();
	}
      }
    }
  }
  hEt = new TH1F("hEt","hEt",450000,-0.1,4.4);
  hEt->Sumw2();

  hcent = new TH1F("hcent","hcent",101,-0.5,100.5);
  hcent->Sumw2();

  hzvtx = new TH1F("hzvtx","hzvtx",300,-300,300);
  hzvtx->Sumw2();
  hxvtx = new TH1F("hxvtx","hxvtx",20000,-10,10);
  hxvtx->Sumw2();
  hyvtx = new TH1F("hyvtx","hyvtx",20000,-10,10);
  hyvtx->Sumw2();
}

void
AnaCorrFourier::SaveHistos()
{
  tmpf->cd();
  tmpf->Write();
  cout << "Saving Finished" << endl;
  cout << " total events: " << nevents << endl;
}
