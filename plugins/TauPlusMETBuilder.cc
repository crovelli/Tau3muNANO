#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

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
    PuppiMet_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("PuppiMet") )},
    //DeepMet_{consumes<pat::METCollection>( cfg.getParameter<edm::InputTag>("DeepMet") )}
    filterBits_{consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("filter_bits"))},
    Filters_{cfg.getParameter<std::vector<std::string>>("filters")}
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
    const edm::EDGetTokenT<edm::TriggerResults> filterBits_;
    //const edm::EDGetTokenT<pat::METCollection> DeepMet_;
  
    std::vector<std::string> Filters_;
    bool debug = false;
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
    
    edm::Handle<edm::TriggerResults> filterBits;
    evt.getByToken(filterBits_, filterBits);
    const edm::TriggerNames &filterNames = evt.triggerNames(*filterBits); 


    // [OUTPUT]
    std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

    for(const pat::CompositeCandidate &tau : *tau_cands){
        // pick the Tau candidate
        TLorentzVector tauCandP4;
        if ( tau.hasUserFloat("fitted_pt")&&tau.hasUserFloat("fitted_eta")&&tau.hasUserFloat("fitted_phi")&&tau.hasUserFloat("fitted_mass")  ){
            if (debug) std::cout << " tau has all the kimetics .. OK" << std::endl;
            tauCandP4.SetPtEtaPhiM( tau.userFloat("fitted_pt"),
                                    tau.userFloat("fitted_eta"),
                                    tau.userFloat("fitted_phi"),
                                    tau.userFloat("fitted_mass") );
        }else{
            std::cout << " [ERROR] tau candidate is fault..." << std::endl;
        }

        // MET P4 in the transverse plane
        TLorentzVector MetP4, PuppiMetP4, DeepMetP4;
        
        MetP4.SetPtEtaPhiM(met.pt(), 0.0, met.phi(), 0);
        PuppiMetP4.SetPtEtaPhiM(PuppiMet.pt(), 0.0, PuppiMet.phi(), 0);
        DeepMetP4.SetPtEtaPhiM(met.corPt(DeepMETcorr), 0.0, met.corPhi(DeepMETcorr), 0);

        // 2 sol. for missing longitudinal momentum
        std::pair<double,double> MET_missPz(longMETsolutions(MetP4,tauCandP4));
        std::pair<double,double> PuppiMET_missPz(longMETsolutions(PuppiMetP4,tauCandP4));
        std::pair<double,double> DeepMET_missPz(longMETsolutions(DeepMetP4,tauCandP4));

        // Tau candidate transverse mass
        float Tau_mT = std::sqrt(2. * tauCandP4.Perp()* MetP4.Pt() * (1 - std::cos(tauCandP4.Phi()-MetP4.Phi())));
        float Tau_Puppi_mT = std::sqrt(2. * tauCandP4.Perp()* PuppiMetP4.Pt() * (1 - std::cos(tauCandP4.Phi()-PuppiMetP4.Phi())));
        float Tau_Deep_mT = std::sqrt(2. * tauCandP4.Perp()* DeepMetP4.Pt() * (1 - std::cos(tauCandP4.Phi()-DeepMetP4.Phi())));

        // W P4 with 2 sol. for missing Pz 
        TLorentzVector WcandP4_min, WcandP4_max, WcandPuppiP4_min, WcandPuppiP4_max, WcandDeepP4_min, WcandDeepP4_max;
        TLorentzVector MetP4_min, MetP4_max, PuppiMetP4_min, PuppiMetP4_max, DeepMetP4_min, DeepMetP4_max;  
        // PF-type1
        MetP4_min.SetPxPyPzE(MetP4.Px(), MetP4.Py(), MET_missPz.first,  std::sqrt(MetP4.Px()*MetP4.Px()+MetP4.Py()*MetP4.Py()+MET_missPz.first*MET_missPz.first) );
        MetP4_max.SetPxPyPzE(MetP4.Px(), MetP4.Py(), MET_missPz.second, std::sqrt(MetP4.Px()*MetP4.Px()+MetP4.Py()*MetP4.Py()+MET_missPz.second*MET_missPz.second) );
        WcandP4_min = MetP4_min + tauCandP4;
        WcandP4_max = MetP4_max + tauCandP4;
        // Puppi MET
        PuppiMetP4_min.SetPxPyPzE(PuppiMetP4.Px(), PuppiMetP4.Py(), PuppiMET_missPz.first, std::sqrt(PuppiMetP4.Px()*PuppiMetP4.Px()+PuppiMetP4.Py()*PuppiMetP4.Py()+PuppiMET_missPz.first*PuppiMET_missPz.first) );
        PuppiMetP4_max.SetPxPyPzE(PuppiMetP4.Px(), PuppiMetP4.Py(), PuppiMET_missPz.second, std::sqrt(PuppiMetP4.Px()*PuppiMetP4.Px()+PuppiMetP4.Py()*PuppiMetP4.Py()+PuppiMET_missPz.second*PuppiMET_missPz.second) );
        WcandPuppiP4_min = PuppiMetP4_min + tauCandP4;
        WcandPuppiP4_max = PuppiMetP4_max + tauCandP4;
        if (debug) std::cout << Form(" [PuppiMET] W-min pT %.3f \t eta %.3f \t phi %.3f E %.3f", WcandPuppiP4_min.Pt(), WcandPuppiP4_min.Eta(), WcandPuppiP4_min.Phi(), WcandPuppiP4_min.E() ) << std::endl;
        if (debug) std::cout << Form(" [PuppiMET] W-max pT %.3f \t eta %.3f \t phi %.3f E %.3f", WcandPuppiP4_max.Pt(), WcandPuppiP4_max.Eta(), WcandPuppiP4_max.Phi(), WcandPuppiP4_max.E() ) << std::endl;
        // Deep MET
        DeepMetP4_min.SetPxPyPzE(DeepMetP4.Px(), DeepMetP4.Py(), DeepMET_missPz.first, std::sqrt(DeepMetP4.Px()*DeepMetP4.Px()+DeepMetP4.Py()*DeepMetP4.Py()+DeepMET_missPz.first*DeepMET_missPz.first) );
        DeepMetP4_max.SetPxPyPzE(DeepMetP4.Px(), DeepMetP4.Py(), DeepMET_missPz.second, std::sqrt(DeepMetP4.Px()*DeepMetP4.Px()+DeepMetP4.Py()*DeepMetP4.Py()+DeepMET_missPz.second*DeepMET_missPz.second) );
        WcandDeepP4_min = DeepMetP4_min + tauCandP4;
        WcandDeepP4_max = DeepMetP4_max + tauCandP4;
        if (debug) std::cout << Form(" [DeepMET] W-min pT %.3f \t eta %.3f \t phi %.3f E %.3f", WcandDeepP4_min.Pt(), WcandDeepP4_min.Eta(), WcandDeepP4_min.Phi(), WcandDeepP4_min.E() ) << std::endl;
        if (debug) std::cout << Form(" [DeepMET] W-max pT %.3f \t eta %.3f \t phi %.3f E %.3f", WcandDeepP4_max.Pt(), WcandDeepP4_max.Eta(), WcandDeepP4_max.Phi(), WcandDeepP4_max.E() ) << std::endl;

        pat::CompositeCandidate TauPlusMET;
        TauPlusMET.setCharge(tau.charge());
       
        // MET filters 
        for (const std::string& filter: Filters_){

           if(debug) std::cout << " ... checking filter " << filter << std::endl;
           bool filterFound = false;
           unsigned int index = filterNames.triggerIndex(filter);
           if(index == filterBits->size()){
              std::cout << " WARNING filter " << filter << " NOT found in TriggerResults" << std::endl;
           }else{ 
              if (filterNames.triggerName(index) == filter) filterFound = true;
              if(debug) std::cout << " this is filter "<< filterNames.triggerName(index) << " i need " << filter << std::endl; 
              if(debug) std::cout << " extracted value is " << filterBits->accept(index) << std::endl;
           }
           //int filter_val = (filterFound ? filterBits->accept(index) : false); 
           int filter_val = -1;
           if(filterFound){
              filter_val = (filterBits->accept(index) ? 1 : 0);
           }
           //tab->addColumnValue<int>(filter, filter_val, filter);
           if(debug) std::cout << " save value " << filter_val << std::endl;
           TauPlusMET.addUserInt(filterNames.triggerName(index), filter_val);

        }// loop on MET-filters 


        // save variables
        TauPlusMET.addUserInt("charge", TauPlusMET.charge());
        // PF MET type1 correction
        TauPlusMET.addUserFloat("MET_pt", met.pt()),
        TauPlusMET.addUserFloat("Tau_mT", Tau_mT),
        TauPlusMET.addUserFloat("METminPz", MET_missPz.first);
        TauPlusMET.addUserFloat("METmaxPz", MET_missPz.second);
        // Puppi correction
        TauPlusMET.addUserFloat("PuppiMET_pt", PuppiMet.pt()),
        TauPlusMET.addUserFloat("Tau_Puppi_mT", Tau_Puppi_mT),
        TauPlusMET.addUserFloat("PuppiMETminPz", PuppiMET_missPz.first);
        TauPlusMET.addUserFloat("PuppiMETmaxPz", PuppiMET_missPz.second);
        // DeepMET correction
        TauPlusMET.addUserFloat("DeepMET_pt", met.corPt(DeepMETcorr)),
        TauPlusMET.addUserFloat("Tau_Deep_mT", Tau_Deep_mT),
        TauPlusMET.addUserFloat("DeepMETminPz", DeepMET_missPz.first);
        TauPlusMET.addUserFloat("DeepMETmaxPz", DeepMET_missPz.second);
        // Tau + MET ~ W candidate
        TauPlusMET.addUserFloat("pt", (WcandP4_min.Pt() == WcandP4_max.Pt() ? WcandP4_min.Pt() : -1.0 ));
        TauPlusMET.addUserFloat("eta_min", WcandP4_min.Eta());
        TauPlusMET.addUserFloat("eta_max", WcandP4_max.Eta());
        TauPlusMET.addUserFloat("phi", WcandP4_min.Phi());
        TauPlusMET.addUserFloat("mass_min", WcandP4_min.M());
        TauPlusMET.addUserFloat("mass_max", WcandP4_max.M());
        // Puppi correction
        TauPlusMET.addUserFloat("Puppi_pt", (WcandPuppiP4_min.Pt() == WcandPuppiP4_max.Pt() ? WcandPuppiP4_min.Pt() : -1.0 ));
        TauPlusMET.addUserFloat("Puppi_eta_min", WcandPuppiP4_min.Eta());
        TauPlusMET.addUserFloat("Puppi_eta_max", WcandPuppiP4_max.Eta());
        TauPlusMET.addUserFloat("Puppi_phi", WcandPuppiP4_min.Phi());
        TauPlusMET.addUserFloat("Puppi_mass_min", WcandPuppiP4_min.M());
        TauPlusMET.addUserFloat("Puppi_mass_max", WcandPuppiP4_max.M());
        TauPlusMET.addUserFloat("Deep_pt", (WcandDeepP4_min.Pt() == WcandDeepP4_max.Pt()? WcandDeepP4_min.Pt() : -1.0 ));
        TauPlusMET.addUserFloat("Deep_eta_min", WcandDeepP4_min.Eta());
        TauPlusMET.addUserFloat("Deep_eta_max", WcandDeepP4_max.Eta());
        TauPlusMET.addUserFloat("Deep_phi", WcandDeepP4_min.Phi());
        TauPlusMET.addUserFloat("Deep_mass_min", WcandDeepP4_min.M());
        TauPlusMET.addUserFloat("Deep_mass_max", WcandDeepP4_max.M());
        TauPlusMET.addUserFloat("mass_nominal", W_MASS);

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
   
    // equation Pw = missP + P
    //
    // (E** - Pz**) missPz** + 2(A*Pz) * missPz + (E** - A**) = 0
    // a = E** - Pz**
    //   missPz** + 2b  * missPz + c = 0
    // b = (A*Pz)/a
    // c = (E** - A**)/a

    double A = 0.5 * ( W_MASS*W_MASS - TriMuP4.M()*TriMuP4.M() ) + TriMuP4.Px()*metP4.Px() + TriMuP4.Py()*metP4.Py();

    double denom = TriMuP4.E()*TriMuP4.E() - TriMuP4.Pz()*TriMuP4.Pz();
    double a = 1.;
    double b = A*TriMuP4.Pz()/denom;
    double c = (TriMuP4.E()*TriMuP4.E()*metP4.Pt()*metP4.Pt() - A*A)/denom;

    double delta = b*b - a*c; // already removed factor 2
    delta = delta > 0 ? std::sqrt(delta) : 0;
    min_missPz = fabs((b - delta)/a) < fabs((b + delta)/a) ? (b - delta)/a : (b + delta)/a;
    max_missPz = fabs((b - delta)/a) > fabs((b + delta)/a) ? (b - delta)/a : (b + delta)/a;
    
    if(verbose){
       std::cout << " --- solve long missing energy ----" << std::endl;
       std::cout << "  Delta \t b \t c" << std::endl;
       std::cout << "  " << delta << "\t" <<  b << "\t" <<  c << std::endl;
       std::cout << " ----------------------------------" << std::endl;
    }
    return  std::make_pair(min_missPz, max_missPz);

}//longMETsolutions()
