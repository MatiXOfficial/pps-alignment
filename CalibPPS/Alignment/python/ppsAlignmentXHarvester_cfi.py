import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

ppsAlignmentXHarvester = DQMEDHarvester("PPSAlignmentXHarvester",
	folder = cms.string("CalibPPS/Common"),
	debug = cms.bool(True)
)