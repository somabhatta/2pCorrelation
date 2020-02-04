#ifndef _EVENT_HH_
#define _EVENT_HH_

#include <iostream>
#include <vector>
#include "Track.C"

//#ifndef __CINT__
// #include "boost/shared_ptr.hpp"
// typedef boost::shared_ptr<Track> TRACK_PTR;
// #else
typedef Track* TRACK_PTR;
// #endif


class Event{

 public:
  Event();
  ~Event();

  void set_npart   (int npart){ntrk    = npart;}
  void set_id      (int var)  {id      = var;}
  void set_zvtx    (float var){zvtx    = var;}
  void set_cent    (float var){cent    = var;}

  int   get_npart() {return ntrk;}
  int   get_id()    {return id;}
  float get_cent()  {return cent;}
  float get_zvtx()  {return zvtx;}

  void   identify();

  typedef std::vector<TRACK_PTR>::iterator TrackIter;

  void AddTrack(float pt, float eta, float phi0, int qual);
  void AddTrack(float pt, float eta, float phi0, float qop, int qual);

  Track* GetTrack(unsigned int i) const;

 private:
  int   id;
  float zvtx;
  float cent;
  int   ntrk;

  std::string ThisName;
  std::vector<TRACK_PTR > Tracks;

};

#endif
