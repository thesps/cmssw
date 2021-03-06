#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "L1Trigger/L1CaloTrigger/interface/L1EGammaEECalibrator.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// we sort the clusters in pt
bool compare_cluster_pt(const l1t::HGCalMulticluster *cl1 ,
                        const l1t::HGCalMulticluster *cl2);

bool compare_cluster_pt(const l1t::HGCalMulticluster *cl1 ,
                        const l1t::HGCalMulticluster *cl2) {
                          return cl1->pt() > cl2->pt();
                        }

int get_eta_bin(const l1t::HGCalMulticluster *cl) {
  float eta_min = 1.;
  float eta_max = 4.;
  unsigned n_eta_bins = 150;
  int eta_bin = floor((fabs(cl->eta()) - eta_min) / ((eta_max - eta_min)/n_eta_bins));
  if (cl->eta() < 0) return -1*eta_bin; // bin 0 doesn't exist
  return eta_bin;
}


int get_phi_bin(const l1t::HGCalMulticluster *cl) {
  float phi_min = -M_PI;
  float phi_max = M_PI;
  unsigned n_phi_bins = 63;
  return floor(fabs(reco::deltaPhi(cl->phi(), phi_min)) / ((phi_max - phi_min)/n_phi_bins));
}

pair<int, int> get_eta_phi_bin(const l1t::HGCalMulticluster *cl) {
  return std::make_pair(get_eta_bin(cl), get_phi_bin(cl));
}

class L1EGammaEEProducer : public edm::EDProducer {
   public:
      explicit L1EGammaEEProducer(const edm::ParameterSet&);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      edm::EDGetToken multiclusters_token_;
      L1EGammaEECalibrator calibrator_;
};

L1EGammaEEProducer::L1EGammaEEProducer(const edm::ParameterSet& iConfig) :
  multiclusters_token_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("Multiclusters"))),
  calibrator_(iConfig.getParameter<edm::ParameterSet>("calibrationConfig")) {
    produces< BXVector<l1t::EGamma> >("L1EGammaCollectionBXVWithCuts");
}


void L1EGammaEEProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  float minEt_ = 0;

  std::unique_ptr< BXVector<l1t::EGamma> > l1EgammaBxCollection(new l1t::EGammaBxCollection);


  // retrieve clusters 3D
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters_h;
  iEvent.getByToken(multiclusters_token_, multiclusters_h);
  const l1t::HGCalMulticlusterBxCollection& multiclusters = *multiclusters_h;

  std::vector<const l1t::HGCalMulticluster *> selected_multiclusters;
  std::map<std::pair<int, int>, std::vector<const l1t::HGCalMulticluster *>> etaphi_bins;

  // here we loop on the TPGs
  for(auto cl3d = multiclusters.begin(0); cl3d != multiclusters.end(0); cl3d++) {
     // std::cout << "-- CL3D is EG: " <<  cl3d->hwQual() << "   "<< cl3d->eta() <<"   "<< cl3d->pt()<<std::endl;
     if(cl3d->hwQual()) {
       if(cl3d->et() > minEt_) {
         int hw_quality = 1; // baseline EG ID passed
         if(fabs(cl3d->eta()) >= 1.52) {
           hw_quality = 2; // baseline EG ID passed + cleanup of transition region
         }

         float calib_factor = calibrator_.calibrationFactor(cl3d->pt(), cl3d->eta());
         l1t::EGamma eg=l1t::EGamma(reco::Candidate::PolarLorentzVector(cl3d->pt()/calib_factor, cl3d->eta(), cl3d->phi(), 0.));
         eg.setHwQual(hw_quality);
         eg.setHwIso(1);
         eg.setIsoEt(-1); // just temporarily as a dummy value
         l1EgammaBxCollection->push_back(0,eg);
         if (hw_quality == 2) {
           selected_multiclusters.push_back(&(*cl3d));
           auto eta_phi_bin = get_eta_phi_bin(&(*cl3d));
           auto bucket = etaphi_bins.find(eta_phi_bin);
           if(bucket == etaphi_bins.end()) {
             std::vector<const l1t::HGCalMulticluster *> vec;
             vec.push_back(&(*cl3d));
             etaphi_bins[eta_phi_bin] = vec;
           } else {
             bucket->second.push_back(&(*cl3d));
           }

         }
      }
    }
  }

  // std::cout << "# of selected HGCAL multiclusters: " << selected_multiclusters.size() << std::endl;

  std::sort(selected_multiclusters.begin(), selected_multiclusters.end(), compare_cluster_pt);
  std::set<const l1t::HGCalMulticluster *> used_clusters;
  for(auto cl3d : selected_multiclusters) {
    if(used_clusters.find(cl3d) == used_clusters.end()) {
      // reco::Candidate::LorentzVector mom = cl3d->p4();
      float pt = cl3d->pt();
      // we drop the Had component of the energy
      if(cl3d->hOverE() != -1) pt = cl3d->pt()/(1+cl3d->hOverE());
      reco::Candidate::PolarLorentzVector mom(pt, cl3d->eta(), cl3d->phi(), 0.);
      // cl3d->pt()/(1+cl3d->hOverE())

      // this is not yet used
      used_clusters.insert(cl3d);
      auto eta_phi_bin = get_eta_phi_bin(cl3d);
      // std::cout << " cl pt: " << cl3d->pt()
      //           << " eta: " << cl3d->eta()
      //           << " phi: " << cl3d->phi()
      //           << " eta, phi bin: " << eta_phi_bin.first << "," << eta_phi_bin.second
      //           << std::endl;

      for(int eta_bin : {eta_phi_bin.first-1, eta_phi_bin.first, eta_phi_bin.first+1}) {
        for(int phi_bin: {eta_phi_bin.second-1, eta_phi_bin.second, eta_phi_bin.second+1}) {
          auto bucket = etaphi_bins.find(std::make_pair(eta_bin, phi_bin));
          if(bucket != etaphi_bins.end()) {
            // this bucket is not empty
            for(auto other_cl_ptr: bucket->second) {
              if(used_clusters.find(other_cl_ptr) == used_clusters.end()) {
                if(fabs(other_cl_ptr->eta() - cl3d->eta()) < 0.02) {
                  if(fabs(reco::deltaPhi(other_cl_ptr->phi(), cl3d->phi())) < 0.1) {

                    // std::cout << "MERGE with cl pt: " << other_cl_ptr->pt()
                    //           << " eta: " << other_cl_ptr->eta()
                    //           << " phi: " << other_cl_ptr->phi()
                    //           <<  std::endl;
                    float pt_other = other_cl_ptr->pt();
                    if(other_cl_ptr->hOverE() != -1) pt_other = other_cl_ptr->pt()/(1+other_cl_ptr->hOverE());
                    mom+=reco::Candidate::PolarLorentzVector(pt_other, other_cl_ptr->eta(), other_cl_ptr->phi(), 0.);
                    used_clusters.insert(other_cl_ptr);
                  }
                }
              }
            }
          }
        }
      }
      float calib_factor = calibrator_.calibrationFactor(mom.pt(), mom.eta());
      l1t::EGamma eg=l1t::EGamma(reco::Candidate::PolarLorentzVector(mom.pt()/calib_factor, mom.eta(), mom.phi(), 0.));
      // l1t::EGamma eg(mom);
      // FIXME: full duplication with HWqual 2
      eg.setHwQual(3);
      eg.setHwIso(1);
      l1EgammaBxCollection->push_back(0,eg);

    }


  }

  iEvent.put(std::move(l1EgammaBxCollection),"L1EGammaCollectionBXVWithCuts");
}



DEFINE_FWK_MODULE(L1EGammaEEProducer);
