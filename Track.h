#ifndef _TRACK_HH_
#define _TRACK_HH_

#include <vector>
#include <string>

class Track{
  
 public:
  Track();
  ~Track() {};
  
  float get_pt()   {return pt;}
  float get_eta()  {return eta;}
  float get_phi0() {return phi0;}
  int   get_qual() {return qual;}
  float get_qop()  {return qop;}
  
  void set_pt   (float apt  )  { pt   = apt;}
  void set_eta  (float aeta )  { eta  = aeta;}
  void set_phi0 (float aphi0)  { phi0 = aphi0;}
  void set_qual (int   aqual)  { qual = aqual;}
  void set_qop  (float aqop )  { qop  = aqop;}

  int isValid(float f);
  int isValid(int i);
  
  Track* clone();

 private:
  float pt;
  float eta;
  float phi0;
  int   qual;  
  float qop;  
};
#endif
