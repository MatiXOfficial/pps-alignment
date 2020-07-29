import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

ppsAlignmentWorker = DQMEDHarvester("PPSAlignmentHarvester",
	folder = cms.string("CalibPPS/Common")
)