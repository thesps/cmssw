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

  typedef typename AlgebraicROOTObject<D,5>::Matrix MatD5;
  typedef typename AlgebraicROOTObject<5,D>::Matrix Mat5D;
  typedef typename AlgebraicROOTObject<D,D>::SymMatrix SMatDD;
  typedef typename AlgebraicROOTObject<D>::Vector VecD;
  using ROOT::Math::SMatrixNoInit;
  double pzSign = tsos.localParameters().pzSign();

  auto && x = tsos.localParameters().vector();
  auto && C = tsos.localError().matrix();

  // projection matrix (assume element of "H" to be just 0 or 1)
  ProjectMatrix<double,5,D>  pf;
 
  // Measurement matrix
  VecD r, rMeas; 
  SMatDD V(SMatrixNoInit{}), VMeas(SMatrixNoInit{});

  KfComponentsHolder holder; 
  holder.template setup<D>(&r, &V, &pf, &rMeas, &VMeas, x, C);
  aRecHit.getKfComponents(holder);

  std::cout << "hit: (" << r[0] <<", "<< r[1] <<")" << std::endl << " x: (";
  for(int i = 0; i < 5; i++)
    std::cout << x[i] << ", ";
  std::cout << ")" << std::endl;
  std::cout << "C:(";
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      std::cout << C(i, j) << ", ";
    }
    std::cout << ")" << std::endl << "(";
  }
  
  r -= rMeas;

  // and covariance matrix of residuals
  SMatDD R = V + VMeas;
  bool ok = invertPosDefMatrix(R);

  // Compute Kalman gain matrix
  AlgebraicMatrix55 M = AlgebraicMatrixID();
  Mat5D K = C*pf.project(R);
  pf.projectAndSubtractFrom(M,K);
 

  // Compute local filtered state vector
  AlgebraicVector5 fsv = x + K * r;
  // Compute covariance matrix of local filtered state vector
  AlgebraicSymMatrix55 fse = ROOT::Math::Similarity(M, C) + ROOT::Math::Similarity(K, V);

  std::cout << "upState: (";
  for(int i = 0; i < 5; i++)
    std::cout << fsv[i] << ", ";
  std::cout << ")" << std::endl;
  std::cout << "(";
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      std::cout << fse(i, j) << ", ";
    }
    std::cout << ")" << std::endl << "(";
  }
  /*
  // expanded similariy
  AlgebraicSymMatrix55 fse; 
  ROOT::Math::AssignSym::Evaluate(fse, (M* C) * ( ROOT::Math::Transpose(M)));
  AlgebraicSymMatrix55 tmp;
  ROOT::Math::AssignSym::Evaluate(tmp, (K*V) * (ROOT::Math::Transpose(K)));
  fse +=  tmp;
  */

  if (ok) {
    return TrajectoryStateOnSurface( LocalTrajectoryParameters(fsv, pzSign),
				     LocalTrajectoryError(fse), tsos.surface(),&(tsos.globalParameters().magneticField()), tsos.surfaceSide() );
  }else {
    edm::LogError("KFUpdator")<<" could not invert martix:\n"<< (V+VMeas);
    return TrajectoryStateOnSurface();
  }
}
}

template<unsigned int D>
TrajectoryStateOnSurface KFUpdator::dfeUpdate(const TrajectoryStateOnSurface& tsos,
                                              const TrackingRecHit& aRecHit) const{
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

TrajectoryStateOnSurface KFUpdator::update(const TrajectoryStateOnSurface& tsos,
                                           const TrackingRecHit& aRecHit) const {
    TrajectoryStateOnSurface updated;
    switch (aRecHit.dimension()) {
        case 1: updated = lupdate<1>(tsos,aRecHit);
                break;
        case 2: updated = lupdate<2>(tsos,aRecHit);
                break;
        case 3: updated = lupdate<3>(tsos,aRecHit);
                break;
        case 4: updated = lupdate<4>(tsos,aRecHit);
                break;
        case 5: updated = lupdate<5>(tsos,aRecHit);
                break;
        default:
        throw cms::Exception("Rec hit of invalid dimension (not 1,2,3,4,5)") <<
         "The value was " << aRecHit.dimension() <<
         ", type is " << typeid(aRecHit).name() << "\n";
    }

    /*std::ofstream f("KFUpdator_FWSW_comparisons.csv", std::fstream::app);
    // Dump the track parameters to a file
    for(int i = 0; i < 5; i++){
      f << updated.localParameters().vector()[i] << ",";
    }

    if(aRecHit.dimension() == 2){
      TrajectoryStateOnSurface dfeUpdated = dfeUpdate<5>(tsos, aRecHit);
      for(int i = 0; i < 5; i++){
        f << dfeUpdated.localParameters().vector()[i] << ", ";
      }
    }
    f << std::endl;
    f.close();*/
    return updated;
}
