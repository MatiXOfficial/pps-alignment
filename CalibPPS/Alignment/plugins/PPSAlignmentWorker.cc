/****************************************************************************
 *
 *  CalibPPS/Alignment.plugins/PPSAlignmentWorker.cc
 *
 *  Description : PPS Alignment DQM worker
 *
 * Authors:
 *  - Jan Kašpar
 *  - Mateusz Kocot
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMOneEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "CondFormats/PPSObjects/interface/PPSAlignmentConfig.h"

#include <vector>
#include <string>

//----------------------------------------------------------------------------------------------------

class PPSAlignmentWorker : public DQMOneEDAnalyzer<>
{
public:
    PPSAlignmentWorker(const edm::ParameterSet &)
    ~PPSAlignmentWorker() override {};

private:
    void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;

    void analyze(const edm::Event &, const edm::EventSetup &) override;
}