#ifndef _TRACKER_KFUPDATOR_H_
#define _TRACKER_KFUPDATOR_H_

/** \class KFUpdator
 * Update trajectory state by combining predicted state and measurement 
 * as prescribed in the Kalman Filter algorithm 
 * (see R. Fruhwirth, NIM A262 (1987) 444). <BR>
 *
 * x_filtered = x_predicted + K * (measurement - H * x_predicted) <BR> 
 *
 * x_filtered, x_predicted    filtered and predicted state vectors <BR>
 * measurement                measurement vector <BR>
 * H "measurement matrix"     projects state vector onto measurement space <BR>
 * K                          Kalman gain matrix <BR>
 * (formulae for K and error matrix of filtered state not shown) <BR>
 *
 * This implementation works for measurements of all dimensions.
 * It relies on CLHEP double precision vectors and matrices for 
 * matrix calculations. <BR>
 *
 * Arguments: TrajectoryState &   predicted state <BR>
 *            RecHit &            reconstructed hit <BR>
 *
 * Initial author: P.Vanlaer 25.02.1999
 * Ported from ORCA.
 *
 *  \author vanlaer, cerati
 */

#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include <fstream>

class KFUpdator final : public TrajectoryStateUpdator {

public:

  // methods of Updator

  KFUpdator(){
    std::ofstream ofile("KFUpdator_FWSW_comparisons.csv");
    ofile << "sw_x[0], sw_x[1], sw_x[2], sw_x[3], sw_x[4], fw_[0], fw_[1], fw_[2], fw_[3], fw_[4]" << std::endl;
  }

  KFUpdator(const KFUpdator &existing){
  }

  TrajectoryStateOnSurface update(const TrajectoryStateOnSurface&,
                                  const TrackingRecHit&) const;

  template<unsigned int D>
  TrajectoryStateOnSurface dfeUpdate(const TrajectoryStateOnSurface&,
                                     const TrackingRecHit&) const;


  virtual KFUpdator * clone() const {
    return new KFUpdator(*this);
    //return kfUpdator;
  }

};

#endif
