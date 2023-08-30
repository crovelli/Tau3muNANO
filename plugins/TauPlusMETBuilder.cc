#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <vector>
#include <string>
#include <cmath>
#include "TLorentzVector.h"

#include "helper.h"


class TauPlusMETBuilder : public edm::global::EDProducer<>{

public:

    typedef std::vector<pat::CompositeCandidate> TauCandCollection;

explicit TauPlusMETBuilder(const edm::ParameterSet& cfg):
    src_{consumes<TauCandCollection>(cfg.getParameter<edm::InputTag>("src"))},
    met_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("met") )},
    PuppiMet_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("PuppiMet") )}
    //DeepMet_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("DeepMet") )}
    {
        produces<pat::CompositeCandidateCollection>("builtWbosons");
    }   
    ~TauPlusMETBuilder() override {}

    void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:

    const edm::EDGetTokenT<TauCandCollection>  src_;
    const edm::EDGetTokenT<pat::METCollection> met_;
    const edm::EDGetTokenT<pat::METCollection> PuppiMet_;
    //const edm::EDGetTokenT<pat::METCollection> DeepMet_;

    bool debug = true;

    std::pair<double, double> longMETsolutions( TLorentzVector&,  TLorentzVector &) const;

};


void TauPlusMETBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const{

    // [INPUT]
    edm::Handle<TauCandCollection> tau_cands;
    evt.getByToken(src_, tau_cands);

    edm::Handle<pat::METCollection> Met;
    evt.getByToken(met_, Met);
    const pat::MET &met = Met->front();// get PF Type-1 corrected MET directly available from miniAOD

    edm::Handle<pat::METCollection> P_Met;
    evt.getByToken(PuppiMet_, P_Met);
    const pat::MET &PuppiMet = P_Met->front();

    pat::MET::METCorrectionLevel DeepMETcorr;
    DeepMETcorr = pat::MET::RawDeepResolutionTune;
    //edm::Handle<pat::METCollection> D_Met;
    //evt.getByToken(DeepMet_, D_Met);
    //const pat::MET &DeepMet = D_Met->front();

    // [OUTPUT]
    std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

    for(const pat::CompositeCandidate &tau : *tau_cands){
        // pick the Tau candidate
        TLorentzVector tauCandP4;
        if ( tau.hasUserFloat("fitted_vc_pt")&&tau.hasUserFloat("fitted_vc_eta")&&tau.hasUserFloat("fitted_vc_phi")&&tau.hasUserFloat("fitted_vc_mass")  ){
            if (debug) std::cout << " tau has all the kimetics .. OK" << std::endl;
            tauCandP4.SetPtEtaPhiM( tau.userFloat("fitted_vc_pt"),
                                    tau.userFloat("fitted_vc_eta"),
                                    tau.userFloat("fitted_vc_phi"),
                                    tau.userFloat("fitted_vc_mass") );
        }else{
            std::cout << " [ERROR] tau candidate is fault..." << std::endl;
        }

        // W boson in the transverse plane
        TLorentzVector MetP4, WcandP4, PuppiMetP4, WcandPuppiP4, DeepMetP4, WcandDeepP4;
        
        MetP4.SetPtEtaPhiM(met.pt(), 0.0, met.phi(), 0);
        PuppiMetP4.SetPtEtaPhiM(PuppiMet.pt(), 0.0, PuppiMet.phi(), 0);
        DeepMetP4.SetPtEtaPhiM(met.corPt(DeepMETcorr), 0.0, met.corPhi(DeepMETcorr), 0);

        WcandP4.SetPtEtaPhiM((MetP4 + tauCandP4).Perp(), 0.0, (MetP4 + tauCandP4).Phi(), W_MASS);
        WcandPuppiP4.SetPtEtaPhiM((PuppiMetP4 + tauCandP4).Perp(), 0.0, (PuppiMetP4 + tauCandP4).Phi(), W_MASS);
        WcandDeepP4.SetPtEtaPhiM((DeepMetP4 + tauCandP4).Perp(), 0.0, (DeepMetP4 + tauCandP4).Phi(), W_MASS);
        pat::CompositeCandidate TauPlusMET;
        TauPlusMET.setCharge(tau.charge());
        
        // missing longitudinal momentum
        std::pair<double,double> MET_missPz(longMETsolutions(MetP4,tauCandP4));
        std::pair<double,double> PuppiMET_missPz(longMETsolutions(PuppiMetP4,tauCandP4));
        std::pair<double,double> DeepMET_missPz(longMETsolutions(DeepMetP4,tauCandP4));

        // save variables
        TauPlusMET.addUserInt("charge", TauPlusMET.charge());
        // PF MET type1 correction
        TauPlusMET.addUserFloat("MET_pt", met.pt()),
        TauPlusMET.addUserFloat("METminPz", MET_missPz.first);
        TauPlusMET.addUserFloat("METmaxPz", MET_missPz.second);
        // Puppi correction
        TauPlusMET.addUserFloat("PuppiMET_pt", PuppiMet.pt()),
        TauPlusMET.addUserFloat("PuppiMETminPz", PuppiMET_missPz.first);
        TauPlusMET.addUserFloat("PuppiMETmaxPz", PuppiMET_missPz.second);
        // DeepMET correction
        TauPlusMET.addUserFloat("DeepMET_pt", met.corPt(DeepMETcorr)),
        TauPlusMET.addUserFloat("DeepMETminPz", DeepMET_missPz.first);
        TauPlusMET.addUserFloat("DeepMETmaxPz", DeepMET_missPz.second);
        // Tau + MET ~ W candidate
        TauPlusMET.addUserFloat("pt", WcandP4.Pt());
        TauPlusMET.addUserFloat("Puppi_pt", WcandPuppiP4.Pt());
        TauPlusMET.addUserFloat("Deep_pt", WcandDeepP4.Pt());
        TauPlusMET.addUserFloat("eta", WcandP4.Eta());
        TauPlusMET.addUserFloat("Puppi_eta", WcandPuppiP4.Eta());
        TauPlusMET.addUserFloat("Deep_eta", WcandDeepP4.Eta());
        TauPlusMET.addUserFloat("phi", WcandP4.Phi());
        TauPlusMET.addUserFloat("Puppi_phi", WcandDeepP4.Phi());
        TauPlusMET.addUserFloat("Deep_phi", WcandDeepP4.Phi());
        TauPlusMET.addUserFloat("mass", WcandP4.M());

        // push in the event
        ret_value->push_back(TauPlusMET);

    }// loop on tau candidates


    evt.put(std::move(ret_value),  "builtWbosons");
}// produce()


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauPlusMETBuilder);

std::pair<double, double> TauPlusMETBuilder::longMETsolutions(TLorentzVector& metP4, TLorentzVector & TriMuP4) const{

    bool verbose = false;
    double min_missPz = -999 , max_missPz = -999;
    double A = 0.5 * ( W_MASS*W_MASS - TriMuP4.M()*TriMuP4.M() ) + TriMuP4.Px()*metP4.Px() + TriMuP4.Py()*metP4.Py();

    double denom = TriMuP4.E()*TriMuP4.E() - TriMuP4.Pz()*TriMuP4.Pz();
    double a = 1.;
    double b = A*TriMuP4.Pz()/denom;
    double c = (TriMuP4.E()*TriMuP4.E()*metP4.Pt()*metP4.Pt() - A*A)/denom;

    double delta = b*b - a*c; // already removed factor 2
    if (delta > 0){
        delta = std::sqrt(delta);
        min_missPz = fabs((b - delta)/a) < fabs((b + delta)/a) ? (b - delta)/a : (b + delta)/a;
        max_missPz = fabs((b - delta)/a) > fabs((b + delta)/a) ? (b - delta)/a : (b + delta)/a;
    }
    
    if(verbose){
       std::cout << " --- solve long missing energy ----" << std::endl;
       std::cout << "  Delta \t b \t c" << std::endl;
       std::cout << "  " << delta << "\t" <<  b << "\t" <<  c << std::endl;
       std::cout << " ----------------------------------" << std::endl;
    }
    return  std::make_pair(min_missPz, max_missPz);

}//longMETsolutions()
