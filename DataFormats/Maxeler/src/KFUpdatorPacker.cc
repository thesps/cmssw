#include "DataFormats/Maxeler/interface/KFUpdatorPacker.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "DataFormats/Math/interface/ProjectMatrix.h"
#include <iostream>

/**
 * Pack a vector of trajectory states into a float array
 * Ordering is vector followed by covariance matrix
 */
void KFUpdatorPacker::pack(float* packed, const std::vector<TrajectoryStateOnSurface>& tsoss){
  typedef std::vector< float > vf;

  const int nFields = nFieldsState;
  int i = 0;

  for(TrajectoryStateOnSurface tsos : tsoss){
    //auto && predictedState = tt.predictedState();
    auto && x = tsos.localParameters().vector();
    auto && C = tsos.localError().matrix();

    vf cVector = unrollSMat<5>(&C);
    float* startAddr = packed + i*nFields;
    std::copy(x.Array(), x.Array() + 5, startAddr);
    std::copy(std::begin(cVector), std::end(cVector), startAddr + 5);
    i++;
  }
}

void KFUpdatorPacker::pack(float* packed, const std::vector<TrajectoryMeasurement>& hits){
  using ROOT::Math::SMatrixNoInit;
  typedef std::vector<float> vf;
  typedef typename AlgebraicROOTObject<2,2>::SymMatrix SMat22;
  typedef typename AlgebraicROOTObject<2>::Vector Vec2;

  const int nFields = 2 + SMatDD_nUnique(2);
  const int nPaddingFields = (4 - (nFields % 4)); // PCIE padding
  int i = 0;
  for(auto hitm : hits){
    TrackingRecHit::ConstRecHitPointer hit = hitm.recHit();
    Vec2 rMeas = asSVector<2>(hit->parameters());
    SMat22 VMeas = asSMatrix<2>(hit->parametersError());
    // TODO get VMeas data
    vf VMeasVector = unrollSMat<2>(&VMeas);
    float* startAddr = packed + i*(nFields + nPaddingFields);
    std::copy(rMeas.Array(), rMeas.Array() + 2, startAddr);
    std::copy(std::begin(VMeasVector), std::end(VMeasVector), startAddr+2);
    for(int j = 0; j < nPaddingFields; j++){
      *(startAddr + nFields + j) = 0.;
    }
    i++;
  }
}

void KFUpdatorPacker::pack(float* packed, const TrajectoryStateOnSurface& tsos,
        const TrackingRecHit& hit) const{
  typedef std::vector< float > vf;
  typedef std::pair< int, vf > pivf;
  typedef std::vector< pivf > vpivf;
  typedef typename AlgebraicROOTObject<5,5>::SymMatrix SMat55;
  typedef typename AlgebraicROOTObject<2,2>::SymMatrix SMat22;
  typedef typename AlgebraicROOTObject<5>::Vector Vec5;
  typedef typename AlgebraicROOTObject<2>::Vector Vec2;
  // For now just assume 5D state, 2D hit
  using ROOT::Math::SMatrixNoInit;
  auto && x = tsos.localParameters().vector();
  auto && C = tsos.localError().matrix();

  ProjectMatrix<double,5,2> pf;

  // Measurement matrix
  Vec2 r, rMeas;
  SMat22 V(SMatrixNoInit{});
  SMat22 VMeas(SMatrixNoInit{});

  KfComponentsHolder holder;
  holder.template setup<2>(&r, &V, &pf, &rMeas, &VMeas, x, C);
  hit.getKfComponents(holder);

  /*std::cout << "state: ";
  for(auto xi : x)
      std::cout << xi << ", ";
  std::cout << std::endl << "covariance: ";
  for(auto ci : unrollSMat<5>(&C))
      std::cout << ci << ", ";
  std::cout << std::endl << "hit: ";
  for(auto ri : r)
      std::cout << ri << ", ";
  std::cout << std::endl << "covariance: ";
  for(auto vi : unrollSMat<2>(&VMeas))
      std::cout << vi << ", ";
  std::cout << std::endl;*/

  //Vec2 params = asSVector<2>(hit.parameters());
  float pArray[2];
  std::copy(r.Array(), r.Array() + 2, pArray);
  vf pVector(pArray, pArray + sizeof(pArray) / sizeof(pArray[0]));
  //SMat22 errors = asSMatrix<2>(hit.parametersError());

  //Vec5 x = tsos.localParameters().vector();
  float xArray[5];
  std::copy(x.Array(), x.Array() + 5, xArray);
  vf xVector(xArray, xArray + sizeof(xArray) / sizeof(xArray[0]));
  //SMat55 C = tsos.localError().matrix();

  vpivf fields;
  fields.push_back(pivf(posR_, pVector));
  fields.push_back(pivf(posV_, unrollSMat<2>(&VMeas)));
  fields.push_back(pivf(posX_, xVector));
  fields.push_back(pivf(posC_, unrollSMat<5>(&C)));

  // sort by order of position
  std::sort(fields.begin(), fields.end(), [](const pivf &left, const pivf &right) {
    return left.first < right.first;
   });

  int i = 0;
  for(pivf pair : fields)
    for(float field : pair.second)
      packed[i++] = field;

  int nFields = 2 + 5 + SMatDD_nUnique(2) + SMatDD_nUnique(5);
  nFields += (4 - (nFields % 4)); // PCIE padding
  while(i < nFields)
    packed[i++] = 1.;
}

TrajectoryStateOnSurface KFUpdatorPacker::unpack(const int nFields, const float* packed,
        const TrajectoryStateOnSurface& tsos) const{
  typedef std::vector< float > vf;
  typedef std::pair< int, vf > pivf;
  typedef std::vector< pivf > vpivf;
  typedef typename AlgebraicROOTObject<5,5>::SymMatrix SMat55;
  typedef typename AlgebraicROOTObject<5>::Vector Vec5;
  // Copy with implicit typecast
  double packedD[nFields]; 
  std::copy(packed, packed + nFields, packedD);
    // For now just assume 5D state, 2D hit
  // Unpack the data into a ROOT vector and matrix
  Vec5 fsv(packedD, 5);
  // Construct the symmetric matrix from vector of floats
  // Get the vector of floats from the data array
  const double* startC = packedD + 5;
  AlgebraicROOTObject<15>::Vector Cvector(startC, 15);
  //vf cUnrolled(startC, startC + sizeof(startC) / sizeof(startC[0]));
  SMat55 fse(Cvector, false);
  return TrajectoryStateOnSurface(LocalTrajectoryParameters(fsv, tsos.localParameters().pzSign()),
          LocalTrajectoryError(fse), tsos.surface(), &(tsos.globalParameters().magneticField()),
          tsos.surfaceSide() );
}

template<unsigned int D>
std::vector<float> KFUpdatorPacker::unrollSMat(const typename AlgebraicROOTObject<D,D>::SymMatrix* mat){
    std::vector<float> unrolled;
    for(unsigned int i = 0; i < D; i++)
        for(unsigned int j = 0; j < D; j++)
            if(j >= i)
                unrolled.push_back((float) (*mat)(i, j));
    return unrolled;
}

/**
 * Return number of unique elements in a symmetric matrix DxD
 */
int KFUpdatorPacker::SMatDD_nUnique(const unsigned int D){
  int nUnique = 0;
  for(unsigned int i = 1; i < (D+1); i++)
    nUnique += i;
  return nUnique;
}

