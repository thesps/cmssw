#include "L1Trigger/Phase2L1Taus/interface/L1PFJetProducer.hh"

L1PFJetProducer::L1PFJetProducer(const edm::ParameterSet& cfg) : 
    _coneSize( cfg.getParameter<double>("coneSize")),
    _nJets( cfg.getParameter<int>("nJets")),
    _l1PFToken( consumes<vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("L1PFObjects")))
{
  produces< L1PFJetCollection >( "L1PFJets" ).setBranchAlias("L1PFJets");
}

void L1PFJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<L1PFJetCollection> newL1PFJetCollection(new L1PFJetCollection);

  edm::Handle<  l1t::PFCandidateCollection > l1PFCandidates;
  iEvent.getByToken( _l1PFToken, l1PFCandidates);
  l1t::PFCandidateCollection work;

  std::copy(l1PFCandidates.begin(), l1PFCandidates.end(), std::back_inserter(work));

  std::sort(work.begin(), work.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

  std::vector<l1t::PFJet> jets;
  jets.resize(_nJets);
 
    for(int nJet = 0; nJet < _nJets; nJet++){
        // Take the first (highest pt) candidate as a seed
        l1t::PFCandidate seed = *work.begin();
        // Get the particles with a _coneSize of the seed
        auto particlesInCone = std::find_if(work.begin(), work.end(), [&seed](const l1t::PFCandidate &part){
                return reco::deltaR<l1t::PFCandidate,l1t::PFCandidate>(seed, part) <= _coneSize;});
        jets.push_back(makeJet(particlesInCone));
    } 
    std::sort(jets.begin(), jets.end(), [](l1t::PFJet i,l1t::PFJet j){return(i.pt() > j.pt());}); 
    std::copy(jets.begin(), jets.end(), std::back_inserter(newL1PFJetCollection));
    iEvent.put( std::move(newL1PFJetCollection) , "L1PFJets" );
}

L1PFJet makeJet(std::vector<l1t::PFCandidate> parts){

    l1t::PFCandidate seed = parts.at(0);

    auto sumpt = [](float a, const l1t::PFCandidate& b){
        return a + b.pt();
    }
    // Sum the pt
    float pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);

    // pt weighted d eta
    std::vector<float> pt_deta;
    pt_deta.resize(parts.size());
    std::transform(parts.begin(), parts.end(), pt_deta.begin(), [](const l1t::PFCandidate &part){
            return part.pt() * (part.eta() - seed.eta());});
    // Accumulate the pt weighted etas. Init to the seed eta, start accumulating at begin()+1 to skip seed
    float eta = std::accumulate(pt_deta.begin()+1, pt_deta.end(), seed.eta()){

    // pt weighted d phi
    std::vector<float> pt_dphi;
    pt_dphi.resize(parts.size());
    std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [](const l1t::PFCandidate &part){
            return part.pt() * (part.phi() - seed.phi());});
    // Accumulate the pt weighted phis. Init to the seed phi, start accumulating at begin()+1 to skip seed
    float phi = std::accumulate(pt_dphi.begin()+1, pt_dphi.end(), seed.phi()){

    PFJet jet(pt, eta, phi);

    return jet;
}
  
/////////////
// DESTRUCTOR
L1PFJetProducer::~L1PFJetProducer()
{
}  


//////////
// END JOB
void L1PFJetProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
}

////////////
// BEGIN JOB
void L1PFJetProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup )
{
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFJetProducer);
