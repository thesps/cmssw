#ifndef _MAXELER_KFUPDATOR_H
#define _MAXELER_KFUPDATOR_H

/** \class KFUpdatorPacker
 * Class to pack a state and hit for Maxeler device
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include <vector>

class KFUpdatorPacker{

  private:
    // Position of each item in packed data
    const int posX_;
    const int posC_;
    const int posR_;
    const int posV_;

  public:

    KFUpdatorPacker(const int posX, const int posC, const int posR, const int posV) : 
        posX_(posX), posC_(posC), posR_(posR), posV_(posV){}

    void pack(float* packed, const TrajectoryStateOnSurface& tsos,
            const TrackingRecHit& hit) const;

    static void pack(float* packed, const std::vector<TrajectoryStateOnSurface>& tsoss) const;
    static void pack(float* packed, const std::vector<TrackingRecHit>& hits) const;

    TrajectoryStateOnSurface unpack(const int nFields, const float* packed, const TrajectoryStateOnSurface& tsos) const;

    template <unsigned int D>
    std::vector<float> unrollSMat(const typename AlgebraicROOTObject<D, D>::SymMatrix* mat) const;

    static int SMatDD_nUnique(const unsigned int D);

    static const int nFieldsHit = 8;
    static const int nFieldsState = 20;
        
};

#endif
