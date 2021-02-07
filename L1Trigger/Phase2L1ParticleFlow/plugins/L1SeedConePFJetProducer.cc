#include "ap_int.h"
#include "ap_fixed.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1SeedConePFJetProducer.h"

l1t::PFJet L1SeedConePFJetProducer::makeJet_SW(const std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const {
  l1t::PFCandidate seed = *parts.at(0);

  auto sumpt = [](float a, const edm::Ptr<l1t::PFCandidate>& b) { return a + b->pt(); };

  // Sum the pt
  float pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);

  // pt weighted d eta
  std::vector<float> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
    return (part->pt() / pt) * (part->eta() - seed.eta());
  });
  // Accumulate the pt weighted etas. Init to the seed eta, start accumulating at begin()+1 to skip seed
  float eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), seed.eta());

  // pt weighted d phi
  std::vector<float> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
    return (part->pt() / pt) * reco::deltaPhi(part->phi(),seed.phi());
  });
  // Accumulate the pt weighted phis. Init to the seed phi, start accumulating at begin()+1 to skip seed
  float phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), seed.phi());

  l1t::PFJet jet(pt, eta, phi);
  for (auto it = parts.begin(); it != parts.end(); it++) {
    jet.addConstituent(*it);
  }

  return jet;
}

std::vector<l1t::PFJet> L1SeedConePFJetProducer::processEvent_SW(std::vector<edm::Ptr<l1t::PFCandidate>>& work) const{
  // The floating point algorithm simulation
  std::sort(work.begin(), work.end(), [](edm::Ptr<l1t::PFCandidate> i, edm::Ptr<l1t::PFCandidate> j) {
    return (i->pt() > j->pt());
  });
  std::vector<l1t::PFJet> jets;
  jets.reserve(_nJets);
  while (!work.empty() && jets.size() < _nJets) {
    // Take the first (highest pt) candidate as a seed
    edm::Ptr<l1t::PFCandidate> seed = work.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<edm::Ptr<l1t::PFCandidate>> particlesInCone;
    std::copy_if(
        work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const edm::Ptr<l1t::PFCandidate>& part) {
          return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
        });
    jets.push_back(makeJet_SW(particlesInCone));
    // remove the clustered particles
    work.erase(std::remove_if(work.begin(),
                              work.end(),
                              [&](const edm::Ptr<l1t::PFCandidate>& part) {
                                return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
                              }),
               work.end());
  }
  return jets;
}

L1SeedConePFJetProducer::L1SeedConePFJetProducer(const edm::ParameterSet& cfg)
    : _coneSize(cfg.getParameter<double>("coneSize")),
      _nJets(cfg.getParameter<unsigned>("nJets")),
      _HW(cfg.getParameter<bool>("HW")),
      _debug(cfg.getParameter<bool>("debug")),
      _l1PFToken(consumes<std::vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("L1PFObjects"))) {
  produces<l1t::PFJetCollection>();
}

void L1SeedConePFJetProducer::produce(edm::StreamID /*unused*/,
                                      edm::Event& iEvent,
                                      const edm::EventSetup& iSetup) const {
  std::unique_ptr<l1t::PFJetCollection> newPFJetCollection(new l1t::PFJetCollection);

  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(_l1PFToken, l1PFCandidates);

  std::vector<edm::Ptr<l1t::PFCandidate>> particles;
  for (unsigned i = 0; i < (*l1PFCandidates).size(); i++) {
    particles.push_back(edm::Ptr<l1t::PFCandidate>(l1PFCandidates, i));
  }

  std::vector<l1t::PFJet> jets;
  if(_HW){
    jets = processEvent_HW(particles);
  }else{
    jets = processEvent_SW(particles);
  }

  std::sort(jets.begin(), jets.end(), [](l1t::PFJet i, l1t::PFJet j) { return (i.pt() > j.pt()); });
  newPFJetCollection->swap(jets);
  iEvent.put(std::move(newPFJetCollection));
}


namespace L1SCJetEmu{
  // Data types used in the FPGA and FPGA-optimized functions

  //static float etaphi_base = 100./128;
  static float etaphi_base = 100./64;
  typedef ap_ufixed<16,14>  pt_t;      // 1 unit = 0.25 GeV; max = 16 TeV
  typedef ap_fixed<10,4>    etaphi_t;   // 1 unit = 0.01;     max = +/- 5.12
  typedef ap_fixed<12,6>    detaphi_t;  // 1 unit = 0.01;     max = +/- 10.24 
  typedef ap_fixed<18,9>    detaphi2_t;  // 1 unit = 0.01;     max = +/- 10.24 
  typedef ap_fixed<22,16> pt_etaphi_t; // type for product of pt with eta or phi
  typedef ap_uint<5> count_t; // type for multiplicity

  // constants for the axis update
  typedef ap_ufixed<18,-2> inv_pt_t;
  static constexpr int N_table_inv_pt = 1024;
  //
  //static const detaphi_t TWOPI = 3.14159 * 2; // 0.78125 is 100 / 128
  //static const detaphi_t PI = 3.14159; // 0.78125 is 100 / 128
  //static const detaphi_t HALFPI = 3.14159 / 2; // 0.78125 is 100 / 128

  static const detaphi_t TWOPI = 3.14159 * 2. * etaphi_base;  //
  static const detaphi_t PI = 3.14159 * etaphi_base;         // 
  static const detaphi_t HALFPI = 3.14159 / 2 * etaphi_base; // 0.78125 is 100 / 128
  //static const detaphi_t RCONE = 0.4 * 100 / 128;
  //static const detaphi_t R2CONE = RCONE * RCONE;
  //
  static const etaphi_t FIDUCIAL_ETA_PHI = 5.11 * etaphi_base;
  //static const etaphi_t FIDUCIAL_ETA_PHI = 5.11;
  static const pt_t JET_PT_CUT = 5;


  constexpr int ceillog2(int x){
    return (x <= 2) ? 1 : 1 + ceillog2((x+1) / 2);
  }

  constexpr int floorlog2(int x){
    return (x < 2) ? 0 : 1 + floorlog2(x / 2);
  }

  constexpr int pow2(int x){
    return x == 0 ? 1 : 2 * pow2(x - 1);
  }

  template<class data_T, int N>
  inline float real_val_from_idx(unsigned i){
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N); // number of address bits for table
    data_T x(0);
    // The MSB of 1 is implicit in the table
    x[x.width-1] = 1;
    // So we can use the next NB bits for real data
    x(x.width-2, x.width-NB-1) = i;
    return (float) x;
  }

  template<class data_T, int N>
  inline unsigned idx_from_real_val(data_T x){
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N); // number of address bits for table
    // Slice the top-1 NB bits of the value
    // the MSB of '1' is implicit, so only slice below that
    ap_uint<NB> y = x(x.width-2, x.width-NB-1);
    return (unsigned) y(NB-1, 0);
  }

  template<class data_T, class table_T, int N>
  void init_invert_table(table_T table_out[N]){
    // The template data_T is the data type used to address the table
    for(unsigned i = 0; i < N; i++){
      float x = real_val_from_idx<data_T, N>(i);
      table_T inv_x = 1 / x;
      table_out[i] = inv_x;
    }
  }

  template<class in_t, class table_t, int N>
  table_t invert_with_shift(in_t in, bool debug=false){
    table_t inv_table[N];
    init_invert_table<in_t, table_t, N>(inv_table);

    // find the first '1' in the denominator
    int msb = 0;
    for(int b = 0; b < in.width; b++){
      if(in[b]) msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width-msb-1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t,N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width-msb-1);
    if(debug){
      std::cout << "           x " << in << ", msb = " << msb << ", shift = " << (in.width-msb) << ", idx = " << idx << std::endl;
      std::cout << "     pre 1 / " << in_shifted << " = " << inv_in << "(" << 1/(float)in_shifted << ")" << std::endl;
      std::cout << "    post 1 / " << in << " = " << out << "(" << 1/(float)in << ")" << std::endl;
    }
    return out;
  }


detaphi_t deltaPhi(l1t::PFCandidate a, l1t::PFCandidate b){
  etaphi_t aphi = etaphi_t(a.phi() * etaphi_base);
  etaphi_t bphi = etaphi_t(b.phi() * etaphi_base);
  detaphi_t dphi = detaphi_t(aphi) - detaphi_t(bphi);
  // phi wrap
  detaphi_t dphi0 = dphi > PI ? detaphi_t(TWOPI - dphi) : detaphi_t(dphi);
  detaphi_t dphi1 = dphi < -PI ? detaphi_t(TWOPI + dphi) : detaphi_t(dphi);
  detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
  return dphiw;
}

 bool deltaR2(l1t::PFCandidate seed, l1t::PFCandidate part, detaphi_t cone2, bool debug){
    // scale the particle eta, phi to hardware units
    etaphi_t seta = etaphi_t(seed.eta() * etaphi_base);
    etaphi_t peta = etaphi_t(part.eta() * etaphi_base);
    detaphi_t deta = detaphi_t(seta) - detaphi_t(peta);
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = deta*deta + dphi*dphi < cone2;
    //bool ret = r2 < cone2;
    if(debug){
      detaphi2_t r2 = detaphi2_t(deta)*detaphi2_t(deta) + detaphi2_t(dphi)*detaphi2_t(dphi);
      std::cout << "  part eta, seed eta: " << etaphi_t(part.eta()*etaphi_base) << ", " << etaphi_t(seed.eta()*etaphi_base) << std::endl;
      std::cout << "  part phi, seed phi: " << etaphi_t(part.phi()*etaphi_base) << ", " << etaphi_t(seed.phi()*etaphi_base) << std::endl;
      std::cout << "  pt, deta, dphi, r2, lt: " << part.pt() << ", " << deta << ", " << dphi << ", " << r2 << ", " << ret << std::endl;
    }
    return ret;
  }
};

l1t::PFJet L1SeedConePFJetProducer::makeJet_HW(const std::vector<l1t::PFCandidate>& parts) const {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  using namespace L1SCJetEmu;
  l1t::PFCandidate seed = parts.at(0);

  // Fine unless we start using saturation, in which case order matters
  auto sumpt = [](pt_t(a), const l1t::PFCandidate& b) { return a + (pt_t) b.pt(); };

  // Sum the pt
  pt_t pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);
  inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(pt, false);

  // pt weighted d eta
  std::vector<pt_etaphi_t> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const l1t::PFCandidate& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    // Scale the eta, phi coordinate to the HW system
    return pt_etaphi_t(pt_t(part.pt()) * detaphi_t(etaphi_t(part.eta()*etaphi_base) - etaphi_t(seed.eta()*etaphi_base)));
  });
  // Accumulate the pt-weighted etas. Init to the 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), pt_etaphi_t(0));
  etaphi_t eta = etaphi_t(seed.eta()*etaphi_base) + etaphi_t(sum_pt_eta * inv_pt);

  // pt weighted d phi
  std::vector<pt_etaphi_t> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const l1t::PFCandidate& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    // Scale the eta, phi coordinate to the HW system
    return pt_etaphi_t(pt_t(part.pt()) * detaphi_t(deltaPhi(part,seed)));
  });
    // Accumulate the pt-weighted etas. Init to the 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), pt_etaphi_t(0));
  etaphi_t phi = etaphi_t(seed.phi()*etaphi_base) + etaphi_t(sum_pt_phi * inv_pt);
  if(_debug){
      std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x){std::cout << "pt_dphi: " << x << std::endl;});
      std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x){std::cout << "pt_deta: " << x << std::endl;});
      std::cout << " sum_pt_eta: " << sum_pt_eta << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
      std::cout << " sum_pt_phi: " << sum_pt_phi << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
      std::cout << " uncorr eta: " << etaphi_t(seed.eta()*etaphi_base) << ", phi: " << etaphi_t(seed.phi()*etaphi_base) << std::endl;
      std::cout << "   corr eta: " << eta << ", phi: " << phi << std::endl;
      std::cout << "         pt: " << pt << std::endl;
  }

  l1t::PFJet jet(pt, float(eta)/etaphi_base, float(phi)/etaphi_base);
  //for (auto it = parts.begin(); it != parts.end(); it++) {
  //  jet.addConstituent(*it);
  //}

  return jet;
}

std::vector<l1t::PFJet> L1SeedConePFJetProducer::processEvent_HW(std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const {
  // The fixed point algorithm emulation
  using namespace L1SCJetEmu;
  std::vector<l1t::PFCandidate> work;
  work.resize(parts.size());
  std::transform(parts.begin(), parts.end(), work.begin(), [](const edm::Ptr<l1t::PFCandidate>& part){ return *part; });
  std::sort(work.begin(), work.end(), [](l1t::PFCandidate i, l1t::PFCandidate j) {
    return (i.pt() > j.pt());
  });
  // It would be nice to transform the inputs to the etaphi_base of the FW here, as in the line below
  // However the phi may wrap around if the etaphi_base > 1, so don't do it...
  //std::for_each(work.begin(), work.end(), [](l1t::PFCandidate& x){x.setP4(math::PtEtaPhiMLorentzVector(pt_t(x.pt()), etaphi_t(x.eta()*etaphi_base), etaphi_t(x.phi()*etaphi_base), x.mass()));});

  detaphi_t rCone2 = detaphi_t(_coneSize * _coneSize * etaphi_base * etaphi_base);
  
  std::vector<l1t::PFJet> jets;
  jets.reserve(_nJets);
  while (!work.empty() && jets.size() < _nJets) {
    // Take the first (highest pt) candidate as a seed
    l1t::PFCandidate seed = work.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<l1t::PFCandidate> particlesInCone;
    std::copy_if(
      work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const l1t::PFCandidate& part) {
        return deltaR2(seed, part, rCone2, false);
      });
    if(_debug){
        std::cout << " Seed: " << pt_t(seed.pt()) << ", " << etaphi_t(seed.eta()*etaphi_base) << ", " << etaphi_t(seed.phi()*etaphi_base) << std::endl;
        for(unsigned i = 0; i < particlesInCone.size(); i++){
            l1t::PFCandidate part = particlesInCone.at(i);
            std::cout << "  Part: " << pt_t(part.pt()) << ", " << etaphi_t(part.eta()*etaphi_base) << ", " << etaphi_t(part.phi()*etaphi_base) << std::endl;
            deltaR2(seed, part, rCone2, true);
           
        }
    }
    jets.push_back(makeJet_HW(particlesInCone));
    // remove the clustered particles
    work.erase(std::remove_if(work.begin(),
                              work.end(),
                              [&](const l1t::PFCandidate& part) {
                                //return deltaR2(*seed, *part, false) <= (detaphi_t) (_coneSize * _coneSize);
                                return deltaR2(seed, part, rCone2, false);
                              }),
               work.end());
  }
  return jets;
}
