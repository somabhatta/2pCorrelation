#include "Event.h"

#include <iostream>
#include <cmath>
//#include <boost/smart_ptr.hpp>

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
  TRACK_PTR ph(new Track());
  //Track* ph= new Track();
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
  //TRACK_PTR ph(new Track());
  Track* ph= new Track();
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
}
