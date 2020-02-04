#include "Track.h"

using namespace std;

Track::Track(){
  set_pt  (-9999);
  set_eta (-9999);
  set_phi0(-9999);
  set_qop (-9999);
  set_qual(-9999);
}

int
Track::isValid(float f)
{
  if (f == -9999)
    {
      return 0;
    }
  return 1;
}

int
Track::isValid(int i)
{
  if (i == -9999)
    {
      return 0;
    }
  return 1;
}

Track*
Track::clone()
{
  Track* mytrack = new Track();

  mytrack->set_pt(pt);
  mytrack->set_eta(eta);
  mytrack->set_phi0(phi0);
  mytrack->set_qop(qop);
  mytrack->set_qual(qual);

  return mytrack;
}
