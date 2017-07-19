#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/MeasurementExtractor.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/invertPosDefMatrix.h"
#include "DataFormats/Math/interface/ProjectMatrix.h"
#include "DataFormats/Maxeler/interface/KFUpdatorPacker.h"
#include <iostream>
// MaxCompiler includes
#include "DataFormats/Maxeler/interface/single_updator_sim.h"
#include <MaxSLiCInterface.h>

namespace {

template <unsigned int D>
TrajectoryStateOnSurface lupdate(const TrajectoryStateOnSurface& tsos,
				           const TrackingRecHit& aRecHit) {

  KFUpdatorPacker packer(0, 1, 2, 3);
  int nFields = 5 + D + KFUpdatorPacker::SMatDD_nUnique(5) + KFUpdatorPacker::SMatDD_nUnique(D);
  nFields += (4 - (nFields % 4)); // Add PCIE padding
  int nFieldsRet = 5 + KFUpdatorPacker::SMatDD_nUnique(5);
  nFieldsRet += (4 - (nFieldsRet % 4)); // Add PCIE padding

  float packed[nFields];
  packer.pack(packed, tsos, aRecHit);

  float packed_ret[nFieldsRet];
  single_updator_sim(1, packed, packed_ret);
  return packer.unpack(nFieldsRet, packed_ret, tsos);
}

}
TrajectoryStateOnSurface KFUpdator::update(const TrajectoryStateOnSurface& tsos,
                                           const TrackingRecHit& aRecHit) const {
    switch (aRecHit.dimension()) {
        case 1: return lupdate<1>(tsos,aRecHit);
        case 2: return lupdate<2>(tsos,aRecHit);
        case 3: return lupdate<3>(tsos,aRecHit);
        case 4: return lupdate<4>(tsos,aRecHit);
        case 5: return lupdate<5>(tsos,aRecHit);
    }
    throw cms::Exception("Rec hit of invalid dimension (not 1,2,3,4,5)") <<
         "The value was " << aRecHit.dimension() <<
        ", type is " << typeid(aRecHit).name() << "\n";
}
