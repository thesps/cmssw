#ifndef L1Trigger_Phase2L1ParticleFlow_L1PFJetCorrectionEmulator_h
#define L1Trigger_Phase2L1ParticleFlow_L1PFJetCorrectionEmulator_h
#include <cmath>
#include "../dataformats/jets.h"

namespace l1ct{

  class L1PFJetCorrectionEmulator{

    typedef ap_ufixed<9, 2, AP_RND_CONV, AP_SAT> sf_t; // 2^9 = 512 scale factor values; 1 unit = 1/2^(9-2)=1/128=0.0078125
    typedef ap_ufixed<pt_t::width, pt_t::iwidth, AP_RND_CONV, AP_SAT> pt_rnd_t;

    inline sf_t floatToSF(float sf) { return sf; };
    inline float SFtoFloat(sf_t sf) { return sf.to_float(); };

    static const int NETABINS = 11;
    static const int NPTBINS = 186;
    static const int NCORR = 2046;

    pt_rnd_t pt_bin_upper_edges[NPTBINS];
    glbeta_t eta_bin_upper_edges[NETABINS];
    sf_t corrections[NCORR];

    public:
    L1PFJetCorrectionEmulator(const float corrections_float[NCORR],
                              const float pt_bin_upper_edges_float[NPTBINS],
                              const float eta_bin_upper_edges_float[NETABINS]){
      for (int i = 0; i < NCORR; i++) { 
        corrections[i] = floatToSF(corrections_float[i]);
      }; 
      
      for (int i = 0; i < NPTBINS; i++) { 
        pt_bin_upper_edges[i] = Scales::makePtFromFloat(pt_bin_upper_edges_float[i]);
      };

      for (int i = 0; i < NETABINS; i++) { 
        eta_bin_upper_edges[i] = Scales::makeGlbEta(eta_bin_upper_edges_float[i]);
      };
    }
   
    Jet jet_energy_correction(Jet jet_in) {

      // Determine idx_eta
      glbeta_t abseta = std::abs(jet_in.hwEta);
      unsigned short int idx_eta = NETABINS-1;
      for(int i = 0; i < NETABINS; i++){
        if(abseta < eta_bin_upper_edges[i]){
          idx_eta = i;
          break;
        }
      }
   
      // Determine idx_pt
      pt_t pt = jet_in.hwPt;
      unsigned short int idx_pt = NPTBINS-1;
      for(int i = 0; i < NPTBINS; i++){
        if(pt < pt_bin_upper_edges[i]){
          idx_pt = i;
          break;
        }
      }

      // Apply the correction
      unsigned int idx = NPTBINS*idx_eta+idx_pt;
      pt_t pt_in = jet_in.hwPt;
      sf_t sf = corrections[idx];
      pt_rnd_t pt_out = sf * pt_in;

      Jet jet_out = jet_in;
      jet_out.hwPt = pt_out;

      return jet_out;
    }
  };

}

#endif
