#ifndef __L1Trigger_VertexFinder_RecoVertexWithTP_h__
#define __L1Trigger_VertexFinder_RecoVertexWithTP_h__


#include <set>
#include <vector>

#include "L1Trigger/VertexFinder/interface/L1TrackTruthMatched.h"
#include "L1Trigger/VertexFinder/interface/TP.h"

#include "L1Trigger/VertexFinder/interface/RecoVertex.h"


namespace l1tVertexFinder {

class RecoVertexWithTP {

public:
  // Fill useful info about tracking particle.
  RecoVertexWithTP(const float z0);
  RecoVertexWithTP(RecoVertex& vertex, std::map<const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>>, const L1TrackTruthMatched*> trackAssociationMap);
  ~RecoVertexWithTP() {}

  /// Tracking Particles in vertex
  const std::vector<const L1TrackTruthMatched*>& tracks() const { return tracks_; }
  /// Tracking Particles in vertex
  const std::set<const TP*>& trueTracks() const { return trueTracks_; }
  /// Number of tracks originating from this vertex
  unsigned int numTracks() const { return tracks_.size(); }
  /// Number of true particles assigned to this vertex
  unsigned int numTrueTracks() const { return trueTracks_.size(); }
  /// Assign fitted track to this vertex
  void insert(const L1TrackTruthMatched* fitTrack);
  /// Compute vertex parameters
  void computeParameters(bool weightedmean = false);
  /// Set z0 position
  void setZ(double z) { z0_ = z; }
  /// Sum ot fitted tracks transverse momentum
  double pT() const { return pT_; }
  /// Vertex z0 position
  double z0() const { return z0_; }
  /// Vertex z0 width
  double z0width() const { return z0width_; }
  /// Clear track vector
  void clear()
  {
    tracks_.clear();
    trueTracks_.clear();
  }
  /// True if primary vertex
  bool PrimaryVertex() const { return pv_; }
  /// Set primary vertex tag
  void isPrimary(bool is) { pv_ = is; }
  /// Contain high-pT track?
  bool hasHighPt() const { return highPt_; }
  /// highest track pT in the vertex
  double highestPt() const { return highestPt_; }
  /// Number of high-pT tracks (pT > 10 GeV)
  unsigned int numHighPtTracks() const { return numHighPtTracks_; }
  /// Vertec MET
  double met() const { return met_; }


private:
  double z0_;
  double z0width_;
  double pT_;
  double met_;
  double metX_;
  double metY_;
  double highestPt_;

  std::vector<const L1TrackTruthMatched*> tracks_;
  std::set<const TP*> trueTracks_;
  bool pv_;
  bool highPt_;
  unsigned int numHighPtTracks_;
};

} // end ns l1tVertexFinder

#endif
