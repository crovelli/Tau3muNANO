// Merges the PFPackedCandidates and Lost tracks
// Check compatibility with 3muons trigger

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"
#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks

using namespace std;

constexpr bool debug = false;

class TrackTriggerSelector : public edm::global::EDProducer <> {
  
public:
    
  explicit TrackTriggerSelector(const edm::ParameterSet &iConfig);
    
  ~TrackTriggerSelector() override {};

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const ;
  
private:
  

  reco::Track fix_track(const reco::Track *tk, double delta) const;  

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;  
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  
  // for trigger match
  std::vector<std::string> HLTPaths_;
  const double drForTriggerMatch_;
  
  // Offline selection  
  const double trkPtCut_;
  const double trkEtaCut_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
};

TrackTriggerSelector::TrackTriggerSelector(const edm::ParameterSet &iConfig):
  //bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  bFieldToken_(esConsumes()),
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  tracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  lostTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  vertexToken_(consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertices" ))), 
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),  
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))), 
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),      
  drForTriggerMatch_(iConfig.getParameter<double>("drForTriggerMatch")),
  trkPtCut_(iConfig.getParameter<double>("trkPtCut")),
  trkEtaCut_(iConfig.getParameter<double>("trkEtaCut")),
  trkNormChiMin_(iConfig.getParameter<int>("trkNormChiMin")),
  trkNormChiMax_(iConfig.getParameter<int>("trkNormChiMax")) 
{
  // produce the SelectedTracks collection (all tracks passing the preselection)
  produces<pat::CompositeCandidateCollection>("SelectedTracks");  
  produces<TransientTrackCollection>("SelectedTransientTracks");  
}

void TrackTriggerSelector::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {

  // Inputs
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from Event" ;
  } 
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  const auto& bField = iSetup.getData(bFieldToken_);
 
  edm::Handle<pat::PackedCandidateCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  iEvent.getByToken(lostTracksToken_, lostTracks);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();       

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &trigNames = iEvent.triggerNames(*triggerBits);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  // for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();

  // Outputs
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out(new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection> trans_tracks_out(new TransientTrackCollection);

  std::vector< std::pair<pat::CompositeCandidate,reco::TransientTrack> > vectrk_ttrk; 
  
  // vectors for Run2 trigger
  std::vector<std::vector<int>> fires;
  std::vector<std::vector<float>> matcher; 
  std::vector<std::vector<float>> DR;
  

  // Loop over tracks and apply preselection  
  std::vector<pat::PackedCandidate> preselTracks;

  for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];

    // arranging cuts for speed 
    if (!trk.hasTrackDetails())  continue;    
    if (abs(trk.pdgId()) != 211) continue; 
    if (trk.pt() < trkPtCut_ )   continue;    
    if (fabs(trk.eta()) > trkEtaCut_) continue;
    if ( (trk.bestTrack()->normalizedChi2() < trkNormChiMin_ && trkNormChiMin_>=0 ) || (trk.bestTrack()->normalizedChi2() > trkNormChiMax_ && trkNormChiMax_>0) ) continue; 

    // high purity requirement applied only in packedCands
    if( iTrk < nTracks && !trk.trackHighPurity()) continue;

    // Selected tracks
    preselTracks.push_back(trk);  
  }
  unsigned int numPresTracks = preselTracks.size();

  if (debug) std::cout << "numPresTracks = " << numPresTracks << std::endl;
  
  // First do trigger match, only for selected tracks
  for( unsigned int iTrk=0; iTrk<numPresTracks; ++iTrk ) {
    pat::PackedCandidate trk = preselTracks[iTrk];
    
    // These vectors have one entry per HLT path
    std::vector<int> frs(HLTPaths_.size(),0);              
    std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
    std::vector<float> temp_DR(HLTPaths_.size(),1000.);

    // Loop over trigger paths
    int ipath=-1;
    for (const std::string &path: HLTPaths_){
      
      if(debug) std::cout << "ipath = " << ipath << ", path = " << path << std::endl;
      if(debug) std::cout << std::endl;
      ipath++;
      
      // Here we loop over trigger objects
      float minDr = 1000.;
      float minPt = 1000.;

      if(debug) std::cout << std::endl;
      if(debug) std::cout << "Now start loop over trigger objects" << std::endl;

      // Loop over trigger objects matching the reference path
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
         if(debug) std::cout << "New object" << std::endl;

         // consider only objects which match the ref path    
         obj.unpackPathNames(trigNames);
         obj.unpackFilterLabels(iEvent, *triggerBits);
         std::vector<std::string> pathNamesAll  = obj.pathNames(false);
         bool isPathExist = false;
         for (unsigned int h = 0, n = pathNamesAll.size(); h < n; ++h) {
            // remove the path version
            string pathNameStart;
            pathNameStart = pathNamesAll[h].substr(0,pathNamesAll[h].find("_v")-0);
            if(debug) std::cout << "In loop over trigger object: this is ipath = " << pathNamesAll[h] << ", I need path = " << path << std::endl;	  
            if(pathNameStart==path) isPathExist = true;     
         }
         if(!isPathExist) continue;
      
         if(debug) std::cout << "Path found -> " << path << std::endl;	  

         // muObjNumberTM: hlt candidate firing the request for tracker muons
         int muObjNumberTM = -1;
         for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	

            if(debug) std::cout << "Event: Filter " << hh << " => " << obj.filterLabels()[hh] << ": pt =  " << obj.pt() << ", eta = " << obj.eta() << ", phi = " << obj.phi() << std::endl;
            if(debug) std::cout << "" << std::endl;	  
	  
            if (path=="HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1" || path=="HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15") {
               if(obj.filterLabels()[hh].find("hltTau3MuTriMuon1filter") != std::string::npos) {  
                  muObjNumberTM = hh;
               }
            }
         }// loop on filters
         
         if(debug && muObjNumberTM>=0) std::cout << "Loop over filters done, my TM filter found" << std::endl;
      
         // here HLT obj vs reco track candidates
         TVector3 trkTV3, objTV3;
         trkTV3.SetPtEtaPhi( trk.pt(), trk.eta(), trk.phi() );
         objTV3.SetPtEtaPhi( obj.pt(), obj.eta(), obj.phi() );
         Float_t deltaR = fabs(trkTV3.DeltaR(objTV3));
    
         // here HLT-trk candidates
         if (muObjNumberTM>=0) {  
            if(debug) std::cout << "DeltaR = " << deltaR << std::endl;
            if(debug) std::cout << "This is a tracker muon HLT candidate" << endl;
            if(deltaR < drForTriggerMatch_){    
               frs[ipath]=1;    
               if (deltaR < minDr){
                  minDr = deltaR;
                  minPt = obj.pt(); 
               }
               if(debug) std::cout << "This object is matched with track: minDr = " << minDr << std::endl;
               if(debug) std::cout << "Offline: " << trk.pt() << " " << trk.eta() << " " << trk.phi() << std::endl;
               if(debug) std::cout << "HLT: " << obj.pt() << " " << obj.eta() << " " << obj.phi() << std::endl;
               if(debug) std::cout << muObjNumberTM << std::endl;
            }
         }

      } // Loop over trigger object
      
      // Here we store the minimum between reco track and all its matched HLT objects for this HLT path
      temp_DR[ipath]=minDr;
      temp_matched_to[ipath]=minPt;
      
    } // Loop over HLT paths

    // One vector per track. Each vector : one element per path, corresponding to the closest HLT object
    fires.push_back(frs);                 // This is used in order to see if a reco trk fired a Trigger (1) or not (0).
    matcher.push_back(temp_matched_to);   // PT of the HLT object matched to the reco trk
    DR.push_back(temp_DR);                // Minimum dR between online and offline

    if(debug) {
      std::cout << std::endl;
      std::cout << "Summary for this reco trk: " << std::endl;
      int size1 = frs.size();
      int size2 = temp_matched_to.size();
      int size3 = temp_DR.size();
      if (size1!=size2 || size1!=size3) 
         std::cout << "problem with size: " << size1 << " " << size2 << " " << size3 << std::endl;
      else {
         std::cout << "size ok: " << size1 << std::endl;	
         for (int jj=0; jj<size1; jj++) std::cout << "fired = " << frs[jj] 
            << ", pT HLT obj = " << temp_matched_to[jj] 
               << ", pT RECO obj = " << trk.pt()
               << ", dR = " << temp_DR[jj] << std::endl;
      }
    } // debug

  } // Loop over reco tracks

  if(debug) std::cout << std::endl;

  // Now check for different reco tracks that are matched to the same HLTObject.
  for(unsigned int path=0; path<HLTPaths_.size(); path++){
     for( unsigned int iTrk=0; iTrk<numPresTracks; iTrk++ ) {
        for(unsigned int it=(iTrk+1); it<numPresTracks; it++){
           if(matcher[iTrk][path]!=1000. && matcher[iTrk][path]==matcher[it][path]){ 
              if(DR[iTrk][path]<DR[it][path]){ 
                 fires[it][path]=0;
                 matcher[it][path]=1000.;
                 DR[it][path]=1000.;                       
              }
              else{
                 fires[iTrk][path]=0;
                 matcher[iTrk][path]=1000.;
                 DR[iTrk][path]=1000.;                       
              }
           }              
        }// loop on trk2
     }// loop on trk1
  }// loop on paths

  
  // Loop over tracks and save all tracks passing the selection
  for( unsigned int iTrk=0; iTrk<numPresTracks; ++iTrk ) {

     pat::PackedCandidate trk = preselTracks[iTrk];

     //const reco::TransientTrack trackTT( (*trk.bestTrack()) , &(*bFieldHandle));
     const reco::TransientTrack trackTT( fix_track( &(*trk.bestTrack()), 1e-8 ), &bField ); 
     if(!trackTT.isValid()) continue; 

     // clean tracks wrt muons passing mediumID 
     int matchedToMediumMuon = 0;
     for (const pat::Muon &imutmp : *muons) {
        for (unsigned int i = 0; i < imutmp.numberOfSourceCandidatePtrs(); ++i) {
           if (! ((imutmp.sourceCandidatePtr(i)).isNonnull() && 
                    (imutmp.sourceCandidatePtr(i)).isAvailable()))   continue;

           const edm::Ptr<reco::Candidate> & source = imutmp.sourceCandidatePtr(i);
           if (source.id() == tracks.id() && source.key() == iTrk){
              if (imutmp.isMediumMuon())  matchedToMediumMuon  = 1;
              break;
           }
        } // on muon sources 
     }// on muons
     if (matchedToMediumMuon) continue;
     if(debug) std::cout << " track " << iTrk << " not TM, good to save pT " << trk.pt() << std::endl;
     // track candidate 
     pat::CompositeCandidate pcand;
     pcand.setP4(trk.p4());
     pcand.setCharge(trk.charge());
     pcand.setVertex(trk.vertex());
     pcand.setPdgId(trk.pdgId());
     pcand.addUserInt("charge", trk.charge());
     pcand.addUserInt("nValidHits", trk.bestTrack()->found());
     pcand.addUserFloat("dZpv", trk.dz(PV.position()));
     pcand.addUserFloat("err_dZpv", trk.dzError());    
     pcand.addUserInt("trackQuality", trk.trackHighPurity()); 

     // compatibility with BS, applied at HLT level
     float trkdr = fabs( (- (trk.vx()-beamSpot.x0()) * trk.py() + (trk.vy()-beamSpot.y0()) * trk.px() ) / trk.pt() );
     pcand.addUserFloat("drForHLT", trkdr);

     // trigger match
     for(unsigned int i=0; i<HLTPaths_.size(); i++){
        pcand.addUserInt(HLTPaths_[i],fires[iTrk][i]);
        std::string namedr = HLTPaths_[i]+"_dr";
        pcand.addUserFloat(namedr,DR[iTrk][i]);  
     }

    // in order to avoid revoking the expensive ttrack builder many times and still have everything sorted, we add them to vector of pairs
    vectrk_ttrk.emplace_back( std::make_pair(pcand,trackTT ) );   
  }// loop on presel tracks
  
  // sort to be uniform with leptons
  std::sort( vectrk_ttrk.begin(), vectrk_ttrk.end(), 
        [] ( auto & trk1, auto & trk2) -> 
        bool {return (trk1.first).pt() > (trk2.first).pt();} 
        );

  // finaly save ttrks and trks to the correct _out vectors
  for ( auto & trk: vectrk_ttrk){
     if (debug) std::cout << "[=] save track with pT " << trk.first.pt() << std::endl;
     tracks_out -> emplace_back( trk.first);
     trans_tracks_out -> emplace_back(trk.second);
  }
  if(debug) std::cout << "[==] Number of saved tracks " << (*tracks_out).size() << std::endl;
  if(debug) std::cout << "[==] Number of saved transient tracks " << (*trans_tracks_out).size() << std::endl;
  iEvent.put(std::move(tracks_out),       "SelectedTracks");
  iEvent.put(std::move(trans_tracks_out), "SelectedTransientTracks");
}

// O. Cerri's code to deal with not positive definite covariance matrices
// https://github.com/ocerri/BPH_RDntuplizer/blob/master/plugins/VtxUtils.cc
/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */

reco::Track TrackTriggerSelector::fix_track(const reco::Track *tk, double delta) const {

  unsigned int i, j;
  double min_eig = 1;

  // Get the original covariance matrix. 
  reco::TrackBase::CovarianceMatrix cov = tk->covariance();

  // Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. 
  TMatrixDSym new_cov(cov.kRows);
  for (i = 0; i < cov.kRows; i++) {
    for (j = 0; j < cov.kRows; j++) {
    // Need to check for nan or inf, because for some reason these
    // cause a segfault when calling Eigenvectors().
    //
    // No idea what to do here or why this happens. 
    if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
	cov(i,j) = 1e-6;
      new_cov(i,j) = cov(i,j);
    }
  }

  // Get the eigenvalues. 
  TVectorD eig(cov.kRows);
  new_cov.EigenVectors(eig);
  for (i = 0; i < cov.kRows; i++)
    if (eig(i) < min_eig)
      min_eig = eig(i);

  // If the minimum eigenvalue is less than zero, then subtract it from the
  // diagonal and add `delta`. 
  if (min_eig < 0) {
    for (i = 0; i < cov.kRows; i++)
      cov(i,i) -= min_eig - delta;
  }

  return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}

DEFINE_FWK_MODULE(TrackTriggerSelector);
