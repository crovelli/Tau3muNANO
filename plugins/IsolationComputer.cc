#include "../interface/IsolationComputer.h"

IsolationComputer::IsolationComputer(
    edm::Handle<PackedCandidatesCollection>& inputPFcandiadtes,
    const double isoRadius,
    const double isoRadiusForHLT,
    const double MaxDZForHLT,
    const double dZpv,
    const double dBetaCone,
    const double dBetaValue, // optimised for Run2
    const double pT_treshold
){
    PFcandCollection_ = inputPFcandiadtes;
    isoRadius_ = isoRadius;
    isoRadiusForHLT_ = isoRadiusForHLT;
    MaxDZForHLT_ = MaxDZForHLT;
    dZpv_ = dZpv;
    dBetaCone_ = dBetaCone;
    dBetaValue_ = dBetaValue;
    pT_treshold_ = pT_treshold;

    muonsToVeto_.clear();
    tracksToVeto_.clear();

    if(debug) std::cout << " [IsoComputer]> pT threshold " << pT_treshold << std::endl;

    // divide PF candidates in neutral / charged / PU-charged
    for(const pat::PackedCandidate &p : *PFcandCollection_){
       if(p.charge() == 0){
          neutral_.push_back(&p);
       }else if (fabs(p.dz()) < dZpv){
          charged_.push_back(&p);
       }else{
          pileup_.push_back(&p);
       }
    }// loop on PF candidates from miniAOD

    // all charged candidates for HLT emulation
    for(const pat::PackedCandidate &p : *PFcandCollection_){
       if(p.charge() != 0){
          if (!p.trackHighPurity()) continue;
          chargedforhlt_.push_back(&p);
       }
    }// loop on PF candidates from miniAOD

    if(debug) std::cout << "Found " << neutral_.size() << " neutral ptls " << std::endl;
    if(debug) std::cout << "Found " << charged_.size() << " charged ptls w dz<0.2" << std::endl;
    if(debug) std::cout << "Found " << pileup_.size()  << " charged-PU ptls in the iso-Cone" << std::endl;  
}// IsolationComputer()

void IsolationComputer::addMuonsToVeto(const std::vector<edm::Ptr<pat::Muon>> inMuToVeto){
    for (size_t it_mu = 0; it_mu < inMuToVeto.size(); ++it_mu){
        if(fabs(inMuToVeto[it_mu]->pdgId()) != 13) continue;
        muonsToVeto_.push_back(inMuToVeto[it_mu]);
    }
    
}// addMuonsToVeto()

void IsolationComputer::addTracksToVeto(const std::vector<edm::Ptr<pat::CompositeCandidate>> inTrkToVeto, const int& pdgID_trk){
    for (size_t it_tk = 0; it_tk < inTrkToVeto.size(); ++it_tk){
        if(fabs(inTrkToVeto[it_tk]->pdgId()) != pdgID_trk) continue;
        tracksToVeto_.push_back(inTrkToVeto[it_tk]);
    }
}// addTracksToVeto()

double IsolationComputer::pTcharged_iso(const reco::Candidate& tau_cand) const {
    double sum_pTcharge = 0;
    bool to_veto = false;
    for (IT ichargedIso = charged_.begin(); ichargedIso != charged_.end(); ++ichargedIso){
        if( (*ichargedIso)->pt() < pT_treshold_ || reco::deltaR(**ichargedIso, tau_cand) > isoRadius_ ) continue;
        if( fabs((*ichargedIso)->pdgId())  == 13){
           std::vector<edm::Ptr<pat::Muon>>::const_iterator it_mu;
           for(it_mu = muonsToVeto_.begin(); it_mu != muonsToVeto_.end(); ++it_mu){
              if(reco::deltaR(**ichargedIso, **it_mu) < DELTA_R_TOMATCH) {
                 to_veto = true;
                 if(debug) std::cout << " xxx veto muon with pT " << (*it_mu)->pt() << " / track pT "<< (*ichargedIso)->pt() << std::endl;
              }
           }
        }else{
           std::vector<edm::Ptr<pat::CompositeCandidate>>::const_iterator it_tk;
           for(it_tk = tracksToVeto_.begin(); it_tk != tracksToVeto_.end(); ++it_tk){
              if(reco::deltaR(**ichargedIso, **it_tk) < DELTA_R_TOMATCH) {
                 to_veto = true;
                 if(debug) std::cout << " xxx veto track with pT " << (*it_tk)->pt() << " / track pT "<< (*ichargedIso)->pt() << std::endl;
              }
           }
        }
        if (!to_veto) sum_pTcharge += (*ichargedIso)->pt();
        to_veto = false;
    }
    // remove muons used to build the tau cand
    if(debug){
      std::cout << " Tau cand has " << tau_cand.numberOfDaughters() << " source ptls" << std::endl;
      std::cout << " = charged pT = " << sum_pTcharge << std::endl;
    }
    return sum_pTcharge;
} //pTcharged_iso()

double IsolationComputer::pTneutral(const reco::Candidate& tau_cand) const {
   double sum_pTneutral = 0;
   bool to_veto = false;
   for (IT ineutralIso = neutral_.begin(); ineutralIso != neutral_.end(); ++ineutralIso){
      if( (*ineutralIso)->pt() < pT_treshold_ || reco::deltaR(**ineutralIso, tau_cand) > isoRadius_ || fabs((*ineutralIso)->pdgId())  != NeutralPDGid) continue;
      // track matchink one to veto?
      std::vector<edm::Ptr<pat::CompositeCandidate>>::const_iterator it_tk;
      //for(it_tk = tracksToVeto_.begin(); it_tk != tracksToVeto_.end(); ++it_tk){
      for(const edm::Ptr<pat::CompositeCandidate> &it_tk : tracksToVeto_){
         if(reco::deltaR(**ineutralIso, *it_tk) < DELTA_R_TOMATCH) {
             if(debug) std::cout << " xxx veto track with pT " << it_tk->pt() << " / track pT "<< (*ineutralIso)->pt() << std::endl;
             continue;
         }
      }//loop on trk to veto
      sum_pTneutral+= (*ineutralIso)->pt();
   }// loop on neutral track

   return sum_pTneutral;
}//pTneutral()

double IsolationComputer::pTcharged_PU(const reco::Candidate& tau_cand) const{
    double sum_pTcharge_PU = 0;
    for (IT ichargedPU = pileup_.begin(); ichargedPU != pileup_.end(); ++ichargedPU){
        if( (*ichargedPU)->pt() < pT_treshold_ || reco::deltaR(**ichargedPU, tau_cand) > dBetaCone_) continue;
        sum_pTcharge_PU += (*ichargedPU)->pt();
    }
    //std::cout << " = charged pT(PU) = " << sum_pTcharge_PU << std::endl;
    return sum_pTcharge_PU;
}// pTcharged_PU()

double IsolationComputer::pTphoton(const reco::Candidate& tau_cand) const{
    double sum_pTphoton = 0;
    for (IT iGamma = neutral_.begin(); iGamma != neutral_.end(); ++iGamma){
        if( (*iGamma)->pdgId() != 22 || (*iGamma)->pt() < pT_treshold_ || reco::deltaR(**iGamma, tau_cand) > isoRadius_) continue;
        sum_pTphoton += (*iGamma)->pt();
    }
    //std::cout << " = gamma pT = " << sum_pTphoton << std::endl;
    return sum_pTphoton;
}// pTcharged_PU()

double IsolationComputer::pTchargedforhlt_iso(const reco::Candidate& tau_cand, float tau_vz) const {
  double sum_pTcharge = 0;
  for (IT ichargedforhltIso = chargedforhlt_.begin(); ichargedforhltIso != chargedforhlt_.end(); ++ichargedforhltIso) {
    if (!(*ichargedforhltIso)) continue;
    if( reco::deltaR( (*ichargedforhltIso)->p4(), tau_cand.p4() ) > isoRadiusForHLT_ ) continue;
    if( fabs((*ichargedforhltIso)->vz() - tau_vz) > MaxDZForHLT_ ) continue;
    sum_pTcharge += (*ichargedforhltIso)->pt();
  }
  // remove tau pT
  sum_pTcharge -=tau_cand.pt();
  if (sum_pTcharge<0) sum_pTcharge=0;
  return sum_pTcharge;
} 

