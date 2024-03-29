#ifndef IsolationComputer_h
#define IsolationComputer_h

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
typedef std::vector<pat::PackedCandidate> PackedCandidatesCollection;
typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"


class IsolationComputer{

    public:
    IsolationComputer(
        edm::Handle<PackedCandidatesCollection>& inputPFcandiadtes, 
        const double isoRadius,
        const double isoRadiusForHLT,
        const double MaxDZForHLT,
        const double dZpv,
        const double dBetaCone,
        const double dBetaValue = 0.2,
        const double pT_treshold = 0.5
        );
    ~IsolationComputer(){}

    // isolation functions
    double pTcharged_iso(const reco::Candidate& tau_cand) const;
    double pTcharged_PU(const reco::Candidate& tau_cand) const;
    double pTphoton(const reco::Candidate& tau_cand) const;
    double pTneutral(const reco::Candidate& tau_cand) const;
    double pTchargedforhlt_iso(const reco::Candidate& tau_cand, float tau_vz) const;

    // veto muons in the candidate
    void addMuonsToVeto(const std::vector<edm::Ptr<pat::Muon>> inMuToVeto);
    void addTracksToVeto(const std::vector<edm::Ptr<pat::CompositeCandidate>> inMuToVeto, const int& pdgID_trk = 211);

    private:

    const bool debug = false;

    // PF candidates
    edm::Handle<PackedCandidatesCollection> PFcandCollection_ ;
    std::vector<const pat::PackedCandidate *> charged_, neutral_, pileup_, chargedforhlt_;

    // isolation parameters
    double isoRadius_;
    double isoRadiusForHLT_;
    double MaxDZForHLT_;
    double dZpv_;
    double dBetaCone_;
    double dBetaValue_;
    double pT_treshold_;

    // to veto
    std::vector<edm::Ptr<pat::Muon>> muonsToVeto_;
    std::vector<edm::Ptr<pat::CompositeCandidate>> tracksToVeto_;
    const double DELTA_R_TOMATCH = 0.0001;
    const int    ElePDGid = 11; 
    const int    MuonPDGid = 13; 
    const int    GammaPDGid = 22; 
    const int    ChargedPDGid = 211; 
    const int    NeutralPDGid = 130; 

};
#endif
