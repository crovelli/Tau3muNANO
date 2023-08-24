import FWCore.ParameterSet.Config as cms
from PhysicsTools.Tau3muNANO.common_cff import *
from PhysicsTools.NanoAOD.met_cff import *

########## inputs preparation ################

# Tau -> 3mu
muonTripletForTau3Mu = cms.EDProducer('TriMuonBuilder',
    src = cms.InputTag('triMuonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('triMuonTrgSelector', 'SelectedTransientMuons'),
    packedCandidatesSrc = cms.InputTag('packedPFCandidates'),
    met = cms.InputTag('slimmedMETs'),
    PuppiMet = cms.InputTag('slimmedMETsPuppi'),
    beamSpot   = cms.InputTag("offlineBeamSpot"),
    # selection definition
    lep1Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep2Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep3Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    preVtxSelection = cms.string('mass() < 3 && abs(charge()) == 1'), # selection for tau candidates pre-fit
    postVtxSelection =  cms.string('userInt("vtx_isValid")'),
    # isolation parameters
    isoRadius = cms.double(0.4), # dR of the isolation cone
    dBetaCone = cms.double(0.8),
    dBetaValue = cms.double(0.2), # optimised for Run2... check validity for Run3
)

# W -> Tau + MET

METforWnuTau3Mu = cms.EDProducer('TauPlusMETBuilder',
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriMuons'),
    met = cms.InputTag('slimmedMETs'),
    PuppiMet = cms.InputTag('slimmedMETsPuppi'),
    #DeepMet = cms.InputTag('deepMetResolutionTuneTable','DeepMETResolutionTune'),
)


################################### Tables #####################################

Tau3MuTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriMuons'),
    cut = cms.string(""),
    name = cms.string("TauTo3Mu"),
    doc = cms.string("Tau Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        mu1_idx = uint('l1_idx'),
        mu2_idx = uint('l2_idx'),
        mu3_idx = uint('l3_idx'),
        charge  = uint('charge'),
        mu1_charge = uint("mu1_charge"),
        mu2_charge = uint("mu2_charge"),
        mu3_charge = uint("mu3_charge"),

        vtx_prob = ufloat("vtx_prob"),
        vtx_chi2 = ufloat("vtx_chi2"),
        vtx_Ndof = ufloat("vtx_Ndof"),
        vtx_isValid = uint("vtx_isValid"),

        fitted_wovc_pt   = ufloat('fitted_wovc_pt'),
        fitted_wovc_eta  = ufloat('fitted_wovc_eta'),
        fitted_wovc_phi  = ufloat('fitted_wovc_phi'),
        fitted_wovc_mass = ufloat("fitted_wovc_mass"),

        fitted_vc_pt   = ufloat('fitted_vc_pt'),
        fitted_vc_eta  = ufloat('fitted_vc_eta'),
        fitted_vc_phi  = ufloat('fitted_vc_phi'),
        fitted_vc_mass = ufloat("fitted_vc_mass"),
        
        MET_pt = ufloat('MET_pt'),
        mT = ufloat("mT"), 
        MET_isPf = uint('MET_isPf'),
        PuppiMET_pt = ufloat('PuppiMET_pt'),
        Puppi_mT = ufloat("Puppi_mT"), 
        PuppiMET_isPf = uint('PuppiMET_isPf'),

        diMuVtxFit_bestProb = ufloat("diMuVtxFit_bestProb"),
        diMuVtxFit_bestMass = ufloat("diMuVtxFit_bestMass"),
        diMuVtxFit_toVeto   = uint("diMuVtxFit_toVeto"),

        iso_ptChargedFromPV = ufloat("iso_ptChargedFromPV"),
        iso_ptChargedFromPU = ufloat("iso_ptChargedFromPU"),
        iso_ptPhotons       = ufloat("iso_ptPhotons"),
        absIsolation        = ufloat("absIsolation"),

        dZmu12 = ufloat('dZmu12'),
        dZmu13 = ufloat('dZmu13'),
        dZmu23 = ufloat('dZmu23'),
        sigLxy_3muVtxBS = ufloat('sigLxy_3muVtxBS'),
        CosAlpha2D_LxyP3mu = ufloat('Cos2D_LxyP3mu'),

        mu1_pt = ufloat("mu1_pt"),
        mu1_eta = ufloat("mu1_eta"),
        mu1_phi = ufloat("mu1_phi"),
        mu1_dr = ufloat("mu1_dr"),
        mu1_trackQuality = uint("mu1_trackQuality"),
        mu2_pt = ufloat("mu2_pt"),
        mu2_eta = ufloat("mu2_eta"),
        mu2_phi = ufloat("mu2_phi"),
        mu2_dr = ufloat("mu2_dr"),
        mu2_trackQuality = uint("mu2_trackQuality"),
        mu3_pt = ufloat("mu3_pt"),
        mu3_eta = ufloat("mu3_eta"),
        mu3_phi = ufloat("mu3_phi"),
        mu3_dr = ufloat("mu3_dr"),
        mu3_trackQuality = uint("mu3_trackQuality"),
        mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        mu1_fired_DoubleMu4_3_LowMass = uint("mu1_fired_DoubleMu4_3_LowMass"),


        mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        mu2_fired_DoubleMu4_3_LowMass = uint("mu2_fired_DoubleMu4_3_LowMass"),

        mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        mu3_fired_DoubleMu4_3_LowMass = uint("mu3_fired_DoubleMu4_3_LowMass"),
    )
)

TauPlusMetTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('METforWnuTau3Mu', 'builtWbosons'),
    cut = cms.string(""),
    name = cms.string("TauPlusMET"),
    doc = cms.string("Tau+MET Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        charge = uint("charge"),
        mass   = ufloat("mass"),
        # PF MET type1 correction
        MET_pt   = ufloat('MET_pt'),
        pt       = ufloat("pt"),
        eta      = ufloat("eta"),
        phi      = ufloat("phi"),
        METminPz = ufloat('METminPz'),
        METmaxPz = ufloat('METmaxPz'),
        # Puppi correction to MET
        PuppiMET_pt  = ufloat('PuppiMET_pt'),
        Puppi_pt     = ufloat("Puppi_pt"),
        Puppi_eta    = ufloat("Puppi_eta"),
        Puppi_phi    = ufloat("Puppi_phi"),
        PuppiMETminPz= ufloat('PuppiMETminPz'),
        PuppiMETmaxPz= ufloat('PuppiMETmaxPz'),
        # Deep MET correction
        DeepMET_pt  = ufloat('DeepMET_pt'),
        Deep_pt     = ufloat("Deep_pt"),
        Deep_eta    = ufloat("Deep_eta"),
        Deep_phi    = ufloat("Deep_phi"),
        DeepMETminPz= ufloat('DeepMETminPz'),
        DeepMETmaxPz= ufloat('DeepMETmaxPz')
    
    )
)

CountMuonTriplets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriMuons'),
)    

########################### Sequencies  ############################

Tau3MuSequence = cms.Sequence(
    (muonTripletForTau3Mu * CountMuonTriplets)
)

Tau3MuTableSequence = cms.Sequence( Tau3MuTable )

TauPlusMetSequence = cms.Sequence( METforWnuTau3Mu )
TauPlusMetTableSequence = cms.Sequence( TauPlusMetTable )