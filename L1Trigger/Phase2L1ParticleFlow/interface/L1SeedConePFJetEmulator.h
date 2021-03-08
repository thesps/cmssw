#ifndef L1SEEDCONEPFJETEMULATOR_H
#define L1SEEDCONEPFJETEMULATOR_H

#include "ap_int.h"
#include "ap_fixed.h"
#include <iostream>
#include <vector>

namespace L1SCJetEmu {
  // Data types and constants used in the FPGA and FPGA-optimized functions
  // This header file is also for use in the standalone FPGA-tools simulation
  // and thus contains no CMSSW/EDM specific content

  //etaphi_base maps physical eta phi units onto bits
  //This way, the least significant bit of etaphi_t is exactly 0.01
  //Even though 0.01 is not a power of 2
  static float etaphi_base = 100. / 64;
  typedef ap_ufixed<16, 14> pt_t;        // 1 unit = 0.25 GeV;
  typedef ap_fixed<10, 4> etaphi_t;      // 1 unit = 0.01;
  typedef ap_fixed<12, 6> detaphi_t;     // type for the difference between etas or phis
  typedef ap_fixed<18, 9> detaphi2_t;    // type for detaphi_t squared
  typedef ap_fixed<22, 16> pt_etaphi_t;  // type for product of pt with deta or phi
  typedef ap_uint<5> count_t;            // type for multiplicity

  class Particle {
  public:
    pt_t hwPt;
    etaphi_t hwEta;
    etaphi_t hwPhi;

    Particle() : hwPt(0), hwEta(0), hwPhi(0) {}
    Particle(pt_t pt, etaphi_t eta, etaphi_t phi) : hwPt(pt), hwEta(eta), hwPhi(phi) {}

    bool operator >= (const Particle &b){
      return hwPt >= b.hwPt;
    }

    pt_t pt() const{
        return hwPt;
    }

    etaphi_t phi() const{
        return hwPhi;
    }

    etaphi_t eta() const{
        return hwEta;
    }
  };

  class Jet : public Particle {
  public:
    Jet(pt_t pt, etaphi_t eta, etaphi_t phi) : Particle(pt, eta, phi) {}
    ap_uint<5> iSeed;
    count_t nCand;
    std::vector<Particle> constituents;

    void addConstituent(const Particle& part){
      constituents.push_back(part);
    }
  };

  struct Config {
  public:
    bool debug;
    float coneSize;
    unsigned nJets;
    Config(bool debug, float coneSize, unsigned nJets) : debug(debug), coneSize(coneSize), nJets(nJets) {}
  };

  // constants for the axis update
  typedef ap_ufixed<18, -2> inv_pt_t;
  static constexpr int N_table_inv_pt = 1024;

  static const detaphi_t TWOPI = 3.14159 * 2. * etaphi_base;
  static const detaphi_t PI = 3.14159 * etaphi_base;
  static const detaphi_t HALFPI = 3.14159 / 2 * etaphi_base;
  //static const detaphi_t RCONE = 0.4 * 100 / 128;
  //static const detaphi_t R2CONE = RCONE * RCONE;
  //
  static const etaphi_t FIDUCIAL_ETA_PHI = 5.11 * etaphi_base;
  //static const etaphi_t FIDUCIAL_ETA_PHI = 5.11;
  static const pt_t JET_PT_CUT = 5;

  constexpr int ceillog2(int x) { return (x <= 2) ? 1 : 1 + ceillog2((x + 1) / 2); }

  constexpr int floorlog2(int x) { return (x < 2) ? 0 : 1 + floorlog2(x / 2); }

  constexpr int pow2(int x) { return x == 0 ? 1 : 2 * pow2(x - 1); }

  template <class data_T, int N>
  inline float real_val_from_idx(unsigned i) {
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    data_T x(0);
    // The MSB of 1 is implicit in the table
    x[x.width - 1] = 1;
    // So we can use the next NB bits for real data
    x(x.width - 2, x.width - NB - 1) = i;
    return (float)x;
  }

  template <class data_T, int N>
  inline unsigned idx_from_real_val(data_T x) {
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    // Slice the top-1 NB bits of the value
    // the MSB of '1' is implicit, so only slice below that
    ap_uint<NB> y = x(x.width - 2, x.width - NB - 1);
    return (unsigned)y(NB - 1, 0);
  }

  template <class data_T, class table_T, int N>
  void init_invert_table(table_T table_out[N]) {
    // The template data_T is the data type used to address the table
    for (unsigned i = 0; i < N; i++) {
      float x = real_val_from_idx<data_T, N>(i);
      table_T inv_x = 1 / x;
      table_out[i] = inv_x;
    }
  }

  template <class in_t, class table_t, int N>
  table_t invert_with_shift(in_t in, bool debug = false) {
    table_t inv_table[N];
    init_invert_table<in_t, table_t, N>(inv_table);

    // find the first '1' in the denominator
    int msb = 0;
    for (int b = 0; b < in.width; b++) {
      if (in[b])
        msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width - msb - 1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t, N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width - msb - 1);
    if (debug) {
      std::cout << "           x " << in << ", msb = " << msb << ", shift = " << (in.width - msb) << ", idx = " << idx
                << std::endl;
      std::cout << "     pre 1 / " << in_shifted << " = " << inv_in << "(" << 1 / (float)in_shifted << ")" << std::endl;
      std::cout << "    post 1 / " << in << " = " << out << "(" << 1 / (float)in << ")" << std::endl;
    }
    return out;
  }

  detaphi_t deltaPhi(Particle a, Particle b) {
    detaphi_t dphi = detaphi_t(a.phi()) - detaphi_t(b.phi());
    // phi wrap
    detaphi_t dphi0 = dphi > PI ? detaphi_t(TWOPI - dphi) : detaphi_t(dphi);
    detaphi_t dphi1 = dphi < -PI ? detaphi_t(TWOPI + dphi) : detaphi_t(dphi);
    detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
    return dphiw;
  }

  bool inCone(Particle seed, Particle part, detaphi_t cone2, bool debug) {
    // scale the particle eta, phi to hardware units
    detaphi_t deta = detaphi_t(seed.eta()) - detaphi_t(part.eta());
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = deta * deta + dphi * dphi < cone2;
    //bool ret = r2 < cone2;
    if (debug) {
      detaphi2_t r2 = detaphi2_t(deta) * detaphi2_t(deta) + detaphi2_t(dphi) * detaphi2_t(dphi);
      std::cout << "  part eta, seed eta: " << part.eta() << ", " << seed.eta() << std::endl;
      std::cout << "  part phi, seed phi: " << part.phi() << ", " << seed.phi() << std::endl;
      std::cout << "  pt, deta, dphi, r2, lt: " << part.pt() << ", " << deta << ", " << dphi << ", " << r2 << ", "
                << ret << std::endl;
    }
    return ret;
  }

Jet makeJet_HW(const std::vector<Particle>& parts, Config config) {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  using namespace L1SCJetEmu;
  Particle seed = parts.at(0);

  // Fine unless we start using saturation, in which case order matters
  auto sumpt = [](pt_t(a), const Particle& b) { return a + b.pt(); };

  // Sum the pt
  pt_t pt = std::accumulate(parts.begin(), parts.end(), pt_t(0), sumpt);
  inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(pt, false);

  // pt weighted d eta
  std::vector<pt_etaphi_t> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    return pt_etaphi_t(part.pt() * detaphi_t(part.eta() - seed.eta()));
  });
  // Accumulate the pt-weighted etas. Init to 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), pt_etaphi_t(0));
  etaphi_t eta = seed.eta() + etaphi_t(sum_pt_eta * inv_pt);

  // pt weighted d phi
  std::vector<pt_etaphi_t> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted dphi
    return pt_etaphi_t(part.pt() * deltaPhi(part, seed));
  });
  // Accumulate the pt-weighted etas. Init to 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), pt_etaphi_t(0));
  etaphi_t phi = seed.phi() + etaphi_t(sum_pt_phi * inv_pt);

  Jet jet(pt, eta, phi);
  for (auto it = parts.begin(); it != parts.end(); it++) {
    jet.addConstituent(*it);
  }

  if (config.debug) {
    std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x) { std::cout << "pt_dphi: " << x << std::endl; });
    std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x) { std::cout << "pt_deta: " << x << std::endl; });
    std::cout << " sum_pt_eta: " << sum_pt_eta << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    std::cout << " sum_pt_phi: " << sum_pt_phi << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
    std::cout << " uncorr eta: " << seed.eta() << ", phi: " << seed.phi() << std::endl;
    std::cout << "   corr eta: " << eta << ", phi: " << phi << std::endl;
    std::cout << "         pt: " << pt << std::endl;
  }

  return jet;
}

std::vector<Jet> emulateEvent(std::vector<Particle>& parts, Config config) {
  // The fixed point algorithm emulation
  using namespace L1SCJetEmu;
  std::vector<Particle> work;
  work.resize(parts.size());
  std::transform(parts.begin(), parts.end(), work.begin(), [](const Particle& part) { return part; });
  std::sort(work.begin(), work.end(), [](Particle i, Particle j) {
    return (i.pt() > j.pt());
  });
  detaphi_t rCone2 = detaphi_t(config.coneSize * config.coneSize * etaphi_base * etaphi_base);

  std::vector<Jet> jets;
  jets.reserve(config.nJets);
  while (!work.empty() && jets.size() < config.nJets) {
    // Take the first (highest pt) candidate as a seed
    Particle seed = work.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<Particle> particlesInCone;
    std::copy_if(
        work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const Particle& part) {
          return inCone(seed, part, rCone2, false);
        });
    if (config.debug) {
      std::cout << "Seed: " << seed.pt() << ", " << seed.eta() << ", " << seed.phi() << std::endl;
      std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part: " << part.pt() << ", " << part.eta() << ", " << part.phi() << std::endl;
        inCone(seed, part, rCone2, true);
      });
    }
    jets.push_back(makeJet_HW(particlesInCone, config));
    // remove the clustered particles
    work.erase(
        std::remove_if(work.begin(),
                       work.end(),
                       [&](const Particle& part) { return inCone(seed, part, rCone2, false); }),
        work.end());
  }
  return jets;
}

};  // namespace L1SCJetEmu

#endif
