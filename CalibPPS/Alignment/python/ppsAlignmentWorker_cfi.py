import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

process = cms.Process('Alignment')

ppsAlignmentWorker = DQMEDAnalyzer("PPSAlignmentWorker",
	tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer")
)