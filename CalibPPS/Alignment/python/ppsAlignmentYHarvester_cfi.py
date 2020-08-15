import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

ppsAlignmentYHarvester = DQMEDHarvester("PPSAlignmentYHarvester",
	folder = cms.string("CalibPPS/Common"),
	debug = cms.bool(True)
)