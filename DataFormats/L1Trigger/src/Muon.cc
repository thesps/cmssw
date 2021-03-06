#include "DataFormats/L1Trigger/interface/Muon.h"

l1t::Muon::Muon()
  : L1Candidate(math::PtEtaPhiMLorentzVector{0., 0., 0., 0.}, 0., 0., 0., 0, 0),
    hwCharge_(0),
    hwChargeValid_(0),
    tfMuonIndex_(-1),
    hwTag_(0),
    debug_(false),
    hwIsoSum_(0),
    hwDPhiExtra_(0),
    hwDEtaExtra_(0),
    hwRank_(0),
    hwEtaAtVtx_(0),
    hwPhiAtVtx_(0),
    etaAtVtx_(0.),
    phiAtVtx_(0.)
{

}

l1t::Muon::Muon( const LorentzVector& p4,
    int pt,
    int eta,
    int phi,
    int qual,
    int charge,
    int chargeValid,
    int iso,
    int tfMuonIndex,
    int tag,
    bool debug,
    int isoSum,
    int dPhi,
    int dEta,
    int rank,
    int hwEtaAtVtx,
    int hwPhiAtVtx,
    double etaAtVtx,
    double phiAtVtx)
  : L1Candidate(p4, pt, eta, phi, qual, iso),
    hwCharge_(charge),
    hwChargeValid_(chargeValid),
    tfMuonIndex_(tfMuonIndex),
    hwTag_(tag),
    debug_(debug),
    hwIsoSum_(isoSum),
    hwDPhiExtra_(dPhi),
    hwDEtaExtra_(dEta),
    hwRank_(rank),
    hwEtaAtVtx_(hwEtaAtVtx),
    hwPhiAtVtx_(hwPhiAtVtx),
    etaAtVtx_(etaAtVtx),
    phiAtVtx_(phiAtVtx)
{

}

l1t::Muon::Muon( const PolarLorentzVector& p4,
    int pt,
    int eta,
    int phi,
    int qual,
    int charge,
    int chargeValid,
    int iso,
    int tfMuonIndex,
    int tag,
    bool debug,
    int isoSum,
    int dPhi,
    int dEta,
    int rank,
    int hwEtaAtVtx,
    int hwPhiAtVtx,
    double etaAtVtx,
    double phiAtVtx)
  : L1Candidate(p4, pt, eta, phi, qual, iso),
    hwCharge_(charge),
    hwChargeValid_(chargeValid),
    tfMuonIndex_(tfMuonIndex),
    hwTag_(tag),
    debug_(debug),
    hwIsoSum_(isoSum),
    hwDPhiExtra_(dPhi),
    hwDEtaExtra_(dEta),
    hwRank_(rank),
    hwEtaAtVtx_(hwEtaAtVtx),
    hwPhiAtVtx_(hwPhiAtVtx),
    etaAtVtx_(etaAtVtx),
    phiAtVtx_(phiAtVtx)
{

}

l1t::Muon::Muon( const l1t::Muon& muon) :
  l1t::L1Candidate(muon)
{
  hwCharge_ = muon.hwCharge();
  hwChargeValid_ = muon.hwChargeValid();
  tfMuonIndex_ = muon.tfMuonIndex();
  hwTag_ = muon.hwTag();
  hwEtaAtVtx_ = muon.hwEtaAtVtx();
  hwPhiAtVtx_ = muon.hwPhiAtVtx();
  etaAtVtx_ = muon.etaAtVtx();
  phiAtVtx_ = muon.phiAtVtx();
  hwIsoSum_ = muon.hwIsoSum();
  hwDPhiExtra_ = muon.hwDPhiExtra();
  hwDEtaExtra_ = muon.hwDEtaExtra();
  hwRank_ = muon.hwRank();
  debug_ = muon.debug();
}

l1t::Muon::~Muon()
{

}

bool l1t::Muon::operator==(const l1t::Muon& rhs) const
{
  return l1t::L1Candidate::operator==(static_cast<const l1t::L1Candidate &>(rhs))
      && hwCharge_ == rhs.hwCharge()
      && hwChargeValid_ == rhs.hwChargeValid()
      && tfMuonIndex_ == rhs.tfMuonIndex()
      && hwEtaAtVtx_ == rhs.hwEtaAtVtx()
      && hwPhiAtVtx_ == rhs.hwPhiAtVtx();
}

