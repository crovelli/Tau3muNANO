//
//class to aplly MET filters
// follow the Run3 recommendation here (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations)
//
#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "helper.h"

using namespace std;

constexpr bool debug = false;

class METFilters : public edm::global::EDProducer <> {
  
public:
    
  explicit METFilters(const edm::ParameterSet &iConfig);
    
  ~METFilters() override {};
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
 
   std::string name_; 
   const edm::EDGetTokenT<edm::TriggerResults> filterBits_;
  
   // for trigger match
   const std::vector<std::string> Filters_;
  
};

METFilters::METFilters(const edm::ParameterSet &iConfig):
   name_(iConfig.getParameter<std::string>("tab_name")),
   filterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   Filters_(iConfig.getParameter<std::vector<std::string>>("filters"))
{
   produces<nanoaod::FlatTable>(); 
}

void METFilters::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const{

   edm::Handle<edm::TriggerResults> filterBits;
   iEvent.getByToken(filterBits_, filterBits);
   const edm::TriggerNames &filterNames = iEvent.triggerNames(*filterBits); 

   // Trigger debug
   if(debug) {
      const edm::TriggerNames &names = iEvent.triggerNames(*filterBits);
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << std::endl;
      std::cout << "\n == FILTERS PATHS= " << std::endl;
      for (unsigned int i = 0, n = filterBits->size(); i < n; ++i) {
         if (filterBits->accept(i)) 
            std::cout << "Event = " << (iEvent.id()).event() << ", Filter " << names.triggerName(i) 
               << ": Pass = " << (filterBits->accept(i)) 
               << ", Was Run = " << (filterBits->wasrun(i))
               << std::endl;
      }

   }

   auto tab  = std::make_unique<nanoaod::FlatTable>(1, name_, true);

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
      tab->addColumnValue<int>(filter, filter_val, filter);
      if(debug) std::cout << " save value " << filter_val << std::endl;

   }// loop on MET-filters 

   if (debug) std::cout << "[+] adding Flat Table with " << tab->nColumns() << " columns" << std::endl;
   iEvent.put(std::move(tab));

}


DEFINE_FWK_MODULE(METFilters);
