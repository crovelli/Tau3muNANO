import FWCore.ParameterSet.Config as cms
from PhysicsTools.Tau3muNANO.common_cff import *
from PhysicsTools.NanoAOD.met_cff import *

########## inputs preparation ################

Path2022=["HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1","HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15","HLT_DoubleMu4_3_LowMass"]
Path=Path2022

# Tau -> 3mu
muonTripletForTau3Mu = cms.EDProducer('TriMuonBuilder',
    src = cms.InputTag('triMuonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('triMuonTrgSelector', 'SelectedTransientMuons'),
    packedCandidatesSrc = cms.InputTag('packedPFCandidates'),
    #met = cms.InputTag('slimmedMETs'),
    #PuppiMet = cms.InputTag('slimmedMETsPuppi'),
    beamSpot   = cms.InputTag("offlineBeamSpot"),
    bits = cms.InputTag("TriggerResults","","HLT"), 
    objects = cms.InputTag("slimmedPatTrigger"),                                         
    # selection definition
    lep1Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep2Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep3Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    preVtxSelection = cms.string('mass() < 3 && abs(charge()) == 1'), # selection for tau candidates pre-fit
    postVtxSelection =  cms.string('userInt("vtx_isValid")'),
    # trigger
    HLTPaths=cms.vstring(Path),                                                                        
    drForTriggerMatch = cms.double(0.1), 
    # isolation parameters
    isoRadius = cms.double(0.4), # dR of the isolation cone
    isoRadiusForHLT = cms.double(0.8), # dR of the isolation cone applied at HLT level
    MaxDZForHLT = cms.double(0.3), # dZ tau-track for isolation
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

        fitted_pt   = ufloat('fitted_pt'),
        fitted_eta  = ufloat('fitted_eta'),
        fitted_phi  = ufloat('fitted_phi'),
        fitted_mass = ufloat("fitted_mass"),
        
        #MET_pt = ufloat('MET_pt'),
        #mT = ufloat("mT"), 
        #MET_isPf = uint('MET_isPf'),
        #PuppiMET_pt = ufloat('PuppiMET_pt'),
        #Puppi_mT = ufloat("Puppi_mT"), 
        #PuppiMET_isPf = uint('PuppiMET_isPf'),

        diMuVtxFit_bestProb = ufloat("diMuVtxFit_bestProb"),
        diMuVtxFit_bestMass = ufloat("diMuVtxFit_bestMass"),
        diMuVtxFit_toVeto   = uint("diMuVtxFit_toVeto"),

        iso_ptChargedFromPV = ufloat("iso_ptChargedFromPV"),
        iso_ptChargedFromPU = ufloat("iso_ptChargedFromPU"),
        iso_ptPhotons       = ufloat("iso_ptPhotons"),
        iso_ptChargedForHLT = ufloat("iso_ptChargedForHLT"),
        absIsolation        = ufloat("absIsolation"),

        dZmu12 = ufloat('dZmu12'),
        dZmu13 = ufloat('dZmu13'),
        dZmu23 = ufloat('dZmu23'),
        dZEasy_mu12 = ufloat('dZEasy_mu12'),
        dZEasy_mu13 = ufloat('dZEasy_mu13'),
        dZEasy_mu23 = ufloat('dZEasy_mu23'),

        DCAmu12 = ufloat('DCAmu12'),
        DCAmu13 = ufloat('DCAmu13'),
        DCAmu23 = ufloat('DCAmu23'),

        vtxFitProbMu12 = ufloat('vtxFitProbMu12'),
        vtxFitProbMu13 = ufloat('vtxFitProbMu13'),
        vtxFitProbMu23 = ufloat('vtxFitProbMu23'),

        sigLxy_3muVtxBS = ufloat('sigLxy_3muVtxBS'),
        CosAlpha2D_LxyP3mu = ufloat('Cos2D_LxyP3mu'),

        fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_dr = ufloat("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_dr"),
        fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_dr = ufloat("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_dr"),

        mu1_pt = ufloat("mu1_pt"),
        mu1_eta = ufloat("mu1_eta"),
        mu1_phi = ufloat("mu1_phi"),
        mu1_dr = ufloat("mu1_drForHLT"),
        mu1_trackQuality = uint("mu1_trackQuality"),
        mu2_pt = ufloat("mu2_pt"),
        mu2_eta = ufloat("mu2_eta"),
        mu2_phi = ufloat("mu2_phi"),
        mu2_dr = ufloat("mu2_drForHLT"),
        mu2_trackQuality = uint("mu2_trackQuality"),
        mu3_pt = ufloat("mu3_pt"),
        mu3_eta = ufloat("mu3_eta"),
        mu3_phi = ufloat("mu3_phi"),
        mu3_dr = ufloat("mu3_drForHLT"),
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
        Tau_mT   = ufloat('Tau_mT'),
        pt       = ufloat("pt"),
        eta      = ufloat("eta"),
        phi      = ufloat("phi"),
        METminPz = ufloat('METminPz'),
        METmaxPz = ufloat('METmaxPz'),
        # Puppi correction to MET
        PuppiMET_pt  = ufloat('PuppiMET_pt'),
        Tau_Puppi_mT = ufloat('Tau_Puppi_mT'),
        Puppi_pt     = ufloat("Puppi_pt"),
        Puppi_eta    = ufloat("Puppi_eta"),
        Puppi_phi    = ufloat("Puppi_phi"),
        PuppiMETminPz= ufloat('PuppiMETminPz'),
        PuppiMETmaxPz= ufloat('PuppiMETmaxPz'),
        # Deep MET correction
        DeepMET_pt  = ufloat('DeepMET_pt'),
        Tau_Deep_mT = ufloat('Tau_Deep_mT'),
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
