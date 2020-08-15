import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

ppsAlignmentXRelativeHarvester = DQMEDHarvester("PPSAlignmentXRelativeHarvester",
	folder = cms.string("CalibPPS/Common"),
	debug = cms.bool(True)
)