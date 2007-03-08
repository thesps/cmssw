#ifndef KinematicParameters_H
#define KinematicParameters_H

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

/**
 * Class to store the 7-vector of
 * particle parameters: (x,y,z,p_x,p_y,p_z,m)
 *
 * Kirill Prokofiev Febrauary 2003
 */


class KinematicParameters{

public:
  KinematicParameters():
              vl(false)
  {}

  KinematicParameters(const AlgebraicVector& pr):
                                par(pr),vl(true)
  {}

/**
 * access methods
 */

  AlgebraicVector vector() const
  {return par;}
  
  GlobalVector momentum() const;
  
  GlobalPoint position() const;
  
  bool isValid() const
  {return vl;}

private:
   AlgebraicVector par;
   bool vl;
};


#endif
