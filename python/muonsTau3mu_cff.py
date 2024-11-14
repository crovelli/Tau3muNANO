import FWCore.ParameterSet.Config as cms
from PhysicsTools.Tau3muNANO.common_cff import *
from PhysicsTools.Tau3muNANO.HLTpathsT3m_cff import Path_Tau3Mu2022
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer

Path2022=["HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1","HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15","HLT_DoubleMu4_3_LowMass"]
Path= Path_Tau3Mu2022

triMuonTrgSelector = cms.EDProducer("TriMuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), 
                                 beamSpot   = cms.InputTag("offlineBeamSpot"),
                                 primaryVtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 
                                 ## trigger match
                                 drForTriggerMatch = cms.double(0.1),     

                                 ## for the output selected collection 
                                 ptMin = cms.double(1.0),                            
                                 absEtaMax = cms.double(2.5),

                                 HLTPaths=cms.vstring(Path)
                             )


# we need at least 3 triggering muons
countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("triMuonTrgSelector", "trgMuons")
)

# muons selection
muonT3mTable = simplePATMuonFlatTableProducer.clone(
    src = cms.InputTag("triMuonTrgSelector:SelectedMuons"),
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for Tau3mu analysis after basic selection"),
    variables = cms.PSet(CandVars, 
        isPFcand = Var("userInt('isPFcand')",bool,doc="muon is global muon"),
        isGlobal = Var("userInt('isGlobal')",bool,doc="muon is global muon"),
        isLoose  = Var("userInt('isLoose')",bool,doc="muon is loose muon"),
        isMedium = Var("userInt('isMedium')",bool,doc="muon is medium muon"),
        isSoft   = Var("userInt('isSoft')",bool,doc="soft cut-based ID"), 
        isSoft_BS = uint('isSoft_BS'),
        isTight  = Var("userInt('isTight')",bool,doc="muon is tight muon"),
        isTight_BS = uint('isTight_BS'),
        isTracker= Var("userInt('isTracker')",bool,doc="muon is tight muon"),
        trackQuality = Var("userInt('trackQuality')",int,doc="trackQuality"),
        z = Var("userFloat('z')",float,doc="track z coordinate"),
        dZpv = Var("userFloat('dZpv')",float,doc="long distance from PV"),
        err_dZpv = Var("userFloat('err_dZpv')",float,doc="long error from PV"),
    ),
)

muonsT3mMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonT3mTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesT3m"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

muonT3mMCTable = cms.EDProducer("CandMCMatchTableProducerT3m",
    src     = muonT3mTable.src,
    mcMap   = cms.InputTag("muonsT3mMCMatchForTable"),
    objName = muonT3mTable.name,
    objType = muonT3mTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('triMuonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsT3mMCMatchForTable')
)


muonT3mSequence = cms.Sequence(triMuonTrgSelector * countTrgMuons)
muonT3mMC = cms.Sequence(muonT3mSequence + muonsT3mMCMatchForTable + selectedMuonsMCMatchEmbedded + muonT3mMCTable)
muonT3mTables = cms.Sequence(muonT3mTable)
