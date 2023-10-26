from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.Tau3muNANO.trgbits_cff import * # modified

## for gen and trigger muon
from PhysicsTools.Tau3muNANO.genparticlesT3m_cff import * # define new
from PhysicsTools.Tau3muNANO.particlelevelT3m_cff import * # define new
from PhysicsTools.Tau3muNANO.triggerObjectsTau3Mu_cff import * #define new
from PhysicsTools.Tau3muNANO.muonsTau3mu_cff import * # define new 
from PhysicsTools.Tau3muNANO.tracksTau3mu_cff import * # define new 


## W collections
from PhysicsTools.Tau3muNANO.Wnu_Tau3Mu import * #define new


nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectTau3MuTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            cms.Sequence(vertexTask) + cms.Sequence(metTablesTask) +           
                            cms.Sequence(globalTablesTask) + cms.Sequence(vertexTablesTask) +
                            triggerObjectTau3MuTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelT3mSequence + genParticleT3mSequence + cms.Sequence(metMCTask) + 
                              cms.Sequence(globalTablesMCTask) + cms.Sequence(genWeightsTableTask) + genParticleT3mTables + lheInfoTable)

def nanoAOD_customizeMuonTriggerTau3Mu(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonT3mSequence + muonT3mTables )
    return process

def nanoAOD_customizeTrackTau3Mu(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trackT3mSequence + trackT3mTables)
    return process

def nanoAOD_customizeTriggerBitsTau3Mu(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeWnuTau3Mu(process):
    process.nanoWnuTau3MuSequence = cms.Sequence( Tau3MuSequence + Tau3MuTableSequence + TauPlusMetSequence + TauPlusMetTableSequence)
    return process


from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        massSearchReplaceAnyInputTag(path, 'trackTrgSelector:SelectedTracks', 'trackT3mMCMatchEmbedded')

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonT3mSequence, process.muonT3mMC)
        path.replace(process.trackT3mSequence, process.trackT3mMC)
